import requests
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from pyteomics import parser
from Bio.SeqUtils import seq3, seq1
import re

########################## functions 
def flatten_dict_hbvar(d, parent_key='', entries=None):
    
    #This function is to flatten out the dictionary to a dataframe
    
    if entries is None:
        entries = []

    for key, value in d.items():
        new_key = f"{parent_key}_{key}" if parent_key else key
        if isinstance(value, dict):
            flatten_dict_hbvar(value, new_key, entries)
        elif isinstance(value, list):
            for i, item in enumerate(value):
                if isinstance(item, dict):
                    flatten_dict_hbvar(item, f"{new_key}_{i}", entries)
                else:
                    entries.append((f"{new_key}_{i}", item))
        else:
            entries.append((new_key, value))

    return entries


def get_ithanet_info(id, info_list):
    #id = '429'
    #info_list = ['Functionality', 'Pathogenicity', 'HGVS Name', 'Hb Name', 'Protein sequence', 'Protein Info']
    url = f'https://www.ithanet.eu/db/ithagenes?ithaID={id}'
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')

        info_out = {}
        for info in info_list if isinstance(info_list, list) else [info_list]:
            #info = 'Protein sequence'
            try:
                if info != 'Protein sequence':
                    info_span = [i for i in soup.find_all('span') if info in i.text][0]
                    info_th = info_span.find_parent('th') if info_span else None
                    info_td = info_th.find_next_sibling('td') if info_th else None
                    info_text = info_td.get_text(strip=True) if info_td else None
                    info_out[info] = info_text
                    
                else:
                    info_span = [i for i in soup.find_all('span') if info in i.text][0]
                    info_p = info_span.find_parent('p') if info_span else None
                    info_br = info_p.find('br') if info_p else None
                    info_text = info_br.text
                    info_out[info] = info_text

            except Exception as e:
                print(f"{info} info is not there for id {id}")
            
        info_df = pd.DataFrame([info_out])
        info_df['id'] = id

        if not info_df.empty:
            return info_df



def mark_mutation(sequence, position):
    # Convert the amino acid at the mutation position to lowercase
    return sequence[:position] + sequence[position].lower() + sequence[position + 1:]

def cleave_sequence_insensitive(sequence, enzyme='trypsin', missed_cleavages=0):
    # Retrieve the original trypsin rule and modify it for case insensitivity
    original_rule = parser.expasy_rules[enzyme]
    # Modified rule to include lowercase 
    modified_rule = '([KkRr](?=[^Pp]))|((?<=Ww)[Kk](?=Pp))|((?<=Mm)[Rr](?=Pp))'
    return list(parser.cleave(sequence, modified_rule, missed_cleavages=missed_cleavages))

def find_mutated_peptide(sequence, mutation_position, enzyme='trypsin', missed_cleavages=0, regex=True):
    # Mark the mutation by converting the mutated residue to lowercase
    
    marked_sequence = mark_mutation(sequence, mutation_position)
        
    # Digest the protein sequence
    peptides = cleave_sequence_insensitive(marked_sequence, enzyme, missed_cleavages)
    
    # Find the peptide containing the lowercase mutation
    peptides_mutation = []
    for peptide in peptides:
        if any(char.islower() for char in peptide):
            peptides_mutation.append(peptide)
    return peptides_mutation


def find_matching_peptides(modified_peptides):
    # this function check if there is all lowercase peptide (means entire deletion of a pepetide), if not, return uppcase peptide (wildetype peptide) and lowercase removed peptide (deletion peptide), 
    
    # there are 2 occasions, (1) when deletion leads to partial missing after trypsin digestion, then the function tries to 
    # to identify most similar but different peptides, (2) when deletion leads to a missing of an entire digested piptide, this 
    # function will just identify the missed fragment peptide. 

    def remove_lowercase(seq):
        return ''.join([char for char in seq if not char.islower()])

    # check completely missing peptides
    # if there is a peptide all are small letters, then it will be completely missing 
    for pep in modified_peptides:
        if pep.islower():            
            return pep.upper(), '-' + pep
        if not pep.isupper():
            wt_pep = pep.upper()
            mt_pep = remove_lowercase(pep)
            return wt_pep, mt_pep
        

def find_mutated_peptide_deletion(sequence, deletions, missed_cleavages=0, enzyme='trypsin'):
    # sequence is a string 
    # deletions is a list of tuples (starting position, deletion amino acids), so far only accept one position deletion 
    # deletion cannot contain K or R

    # a function to count how many cleavage site
    def count_matches(matches):
        count = 0
        for i in matches:
            for j in i:
                if j:
                    count += 1
        return count 

    # check K and R
    for idx, item in enumerate(deletions):
        del_aa = list(item)[1]
        if any([i == 'K' or i == 'R' for i in del_aa]):
            # there is K-R in deletion
            return None, None 

    modified_sequence = sequence
    deletion_points = []

    for position, amino_acids in sorted(deletions, reverse=True):
        #print(position, amino_acids)
        # position = 6
        # amino_acids = 'D'
        start = position - 1  # Adjust from 1-based to 0-based indexing if necessary
        if sequence[start:start + len(amino_acids)] != amino_acids:
            raise ValueError(f"Deletion sequence mismatch: expected {amino_acids} at position {position}, but found {sequence[start:start + len(amino_acids)]}.")

        end = start + len(amino_acids)
        deletion_points.append((start, end))
        modified_sequence = modified_sequence[:start] + modified_sequence[start:end].lower() + modified_sequence[end:]


    original_peptides_cleavage0 = cleave_sequence_insensitive(sequence, enzyme, missed_cleavages=0)
    modified_peptides_cleavage0 = cleave_sequence_insensitive(modified_sequence, enzyme, missed_cleavages=0)
    
    
    wt_pep_cleavage0, mt_pep_cleavage0 = find_matching_peptides(modified_peptides=modified_peptides_cleavage0)

    if missed_cleavages == 0:
        return [wt_pep_cleavage0], [mt_pep_cleavage0]
    else:
        original_peptides = cleave_sequence_insensitive(sequence, enzyme, missed_cleavages)
        modified_peptides = cleave_sequence_insensitive(modified_sequence, enzyme, missed_cleavages)
        modified_rule = '([KkRr](?=[^Pp]))|((?<=Ww)[Kk](?=Pp))|((?<=Mm)[Rr](?=Pp))'
        modified_peptides_subset = [] # this subset conains exactly missed_cleavages number of matches
        for pep in modified_peptides:
            #pep = 'GHGKL'
            matches = re.findall(modified_rule, pep)
            length = count_matches(matches)
            if length == missed_cleavages:
                modified_peptides_subset.append(pep)
        
        wt_pep_cleavage, mt_pep_cleavage = find_matching_peptides(modified_peptides=modified_peptides_subset)
        return [wt_pep_cleavage], [mt_pep_cleavage]


def get_position(string, pattern):
    # string = 'β 51(D2) Pro>Ser AND β 52(D3) Asp>Asn'
    # string = 'β 51(D2) Pro>Ser'
    # pattern=r'\b(\d+)\('
    if isinstance(string, str):
        match = re.findall(pattern, string)
        if match:
            return match
        else:
            return None 
        

def get_wildtype_mutated_aa(string):
    import re
    from Bio.SeqUtils import seq3, seq1
    '''
    Parses mutation descriptions to extract wild type and mutated amino acids.
    Supports patterns like A>B, A->0, A-B->0, and potentially 0>A.

    Arguments:
    string -- str: A string containing the mutation descriptions.

    Returns:
    tuple: Two lists containing the wild type and mutated amino acids respectively.
    '''
    if not isinstance(string, str):
        return [], []  # If input is not a string, return empty lists

    # Define regex patterns for the mutations
    pattern1 = r'[A-Za-z]+(?:-[A-Za-z]+)*->0'  # Deletions, e.g., 'Glu-Val-Gly-Gly->0'
    pattern2 = r'\b[A-Z][a-z]{2}>[A-Z][a-z]{2}\b'  # Point mutation, e.g., 'Ser>Cys'
    
    # Lists to hold the wild type and mutated amino acids
    wildtype_aa = []
    mutated_aa = []

    # Check for point mutations
    matches = re.findall(pattern2, string)
    for match in matches:
        parts = match.split('>')
        wildtype_aa.append(seq1(parts[0]))
        mutated_aa.append(seq1(parts[1]))

    # Check for deletions
    matches = re.findall(pattern1, string)
    for match in matches:
        match_split = match.split('->')
        wildtype_aa.append(seq1(re.sub('-','',match_split[0])))
        mutated_aa.append('del')  # Indicate deletion

    # Final condition to handle cases where no mutations were found
    if not wildtype_aa and not mutated_aa:
        wildtype_aa = ['X']  # or use [''] if empty strings are preferred
        mutated_aa = ['X']  # or use [''] if empty strings are preferred

    return wildtype_aa, mutated_aa


# select peptide based on both mutation numbers and length
def select_peptide(peptides):
    # This function counts the mutations in a peptide
    def count_mutations(peptide):
        return sum(1 for char in peptide if char.islower())
    
    # This list will store the count of mutations and the peptide length
    peptide_info = [(count_mutations(peptide), len(peptide), peptide) for peptide in peptides]
    
    # Filter to find all peptides with more than one mutation
    multiple_mutations = [info for info in peptide_info if info[0] > 1]
    
    if multiple_mutations:
        # If there are multiple, select the one with the most mutations, or the longest if tied
        return max(multiple_mutations, key=lambda x: (x[0], x[1]))[2]
    else:
        # If no peptides have multiple mutations, select the longest peptide
        return max(peptide_info, key=lambda x: x[1])[2]
    

def choose_mutation_position(A, B):
    # here A is more trustable than B. A can be common name. B can be Protein Info
    # regarding position common name is more correct.
    if A is None and B is None:
        return None  
    elif A is None:
        return B
    elif B is None:
        return A
    else:
        if len(A) == len(B):
            if A != B:
                return A
            else:
                return A  # Assuming A by default if A and B are the same
        else:
            return A if len(A) > len(B) else B

def get_sequence_fasta(file, id):
    from Bio import SeqIO

    for line in SeqIO.parse(file, 'fasta'):
        
        if id in line.id:
            return str(line.seq)
        
    return None

# Function to delete amino acids or sequences from the protein sequence
def delete_amino_acids(sequence, deletions):
    # deletion is a tuble (list of starting positions, list of deletion sequences)
    # position is one based, including leading M
    sequence_list = list(sequence)
    adjustments = []
    for position, amino_acid_sequence in sorted(deletions, reverse=True):
        
        start_index = position - 1
        end_index = start_index + len(amino_acid_sequence)
        if sequence_list[start_index:end_index] == list(amino_acid_sequence):
            adjustments.append((start_index, end_index))
        else:
            raise KeyError('position and amino acids not consistent')
    for start_index, end_index in adjustments:
        del sequence_list[start_index:end_index]
    return ''.join(sequence_list)

# def process_info_trypsin(row, hbvar_id_seq_hgvs_hb, hemoglobin_genes_fasta):
#     # outdated, not working 

#     '''
#     improvements:
#     1, when there are multiple pepeitde containing mutated aa, for now we choose the one has the most most mutated aa and the longest, we could also keep all these peptides. 
#     '''
#     #row = itha_prot.iloc[1,:]
#     id = row['IthaID']
#     check_protein_info = row['Protein Info']
#     check_sequence = row['Protein sequence']
#     check_common_name = row['Common Name']
#     check_hgvs_name = row['HGVS Name']
#     check_hb_name = row['Hb Name']

#     sequence = []
#     mutated_peptides = []
#     wildtype_peptides = []
#     warnings = []

#     ############################################## first make mutation position #################################################################
#     # make mutation info, use Common Name first, then Protein Info, (as I found they can diverge e.g. ithanet 3651, CD position is more often correct maybe)
#     if pd.notna(check_protein_info) and pd.isna(check_common_name):
#         mutation_info = row['Protein Info'] 
#         wildtype_mutated_aa = get_wildtype_mutated_aa(string=row['Protein Info']) # aa mutation info in Protein Info is better formatted
#         mutation_position = get_position(mutation_info, pattern=r'\b(?:\D*)(\d+)\(')
        
#     elif pd.notna(check_common_name) and pd.isna(check_protein_info):
#         mutation_info = row['Common Name'] # make a variable to facilitate following steps
#         wildtype_mutated_aa = get_wildtype_mutated_aa(string=row['Common Name'])
#         mutation_position = get_position(row['Common Name'], pattern=r'CD\s+(\d+)\s+') 

#         # when both info is present, mutated aa is more correct in protein info 
#     elif pd.notna(check_protein_info) and pd.notna(check_common_name):
#         mutation_info = row['Protein Info'] 
#         wildtype_mutated_aa = get_wildtype_mutated_aa(string=row['Protein Info']) # aa mutation info in Protein Info is better formatted
#         mutation_position_protein_info = get_position(mutation_info, pattern=r'\b(?:\D*)(\d+)\(')
#         mutation_position_common_name = get_position(row['Common Name'], pattern=r'[Cc][Dd]\s+(\d+)\s+') 
#         # when position info conficts, choose the correct info from common name
#         mutation_position = choose_mutation_position(mutation_position_protein_info, mutation_position_common_name)

#     # if no mutation info, then skip 
#     if pd.isna(check_protein_info) and pd.isna(check_common_name):
#         warnings.append('No mutation info')
#         return wildtype_peptides, mutated_peptides, warnings  # Skip the rest of the loop for this iteration
    
#     # if we see mutation_position is None, then we skip this iteration    
#     if not isinstance(mutation_position, list):
#         warnings.append('No mutation position info')
#         return wildtype_peptides, mutated_peptides, warnings
    
#     # if mutated aa is not point mutation, not multiple mutations, not deletion, but has other category as 'X', then skip
#     if any([i == 'X' for i in wildtype_mutated_aa[1]]):
#         warnings.append('Exceptional mutation case')
#         return wildtype_peptides, mutated_peptides, warnings
        

#     # if mutation position is multiple, and it is deletion, then skip 
#     if len(mutation_position) > 1 and any([i == 'del' for i in wildtype_mutated_aa[1]]):
#         warnings.append('Complex deletions')
#         return wildtype_peptides, mutated_peptides, warnings

#     # Convert mutation positions to integers 
#     mutation_position = [int(pos) - 1 for pos in mutation_position] # minus position by 1, since we will remove the leading M from the sequence 

#     # # test mutation condition 
#     # 1, point mutation (accept)
#     # 2, multiple point mutation (accept) 
#     # 3, single deletion (accept)
#     # 4, multiple deletion (no)
#     # 5, mixture (no)
    
#     # if multple position and mutations,and it is not just point mutations, then skip 
#     mutated_aa = get_wildtype_mutated_aa(string=mutation_info)[1]
#     if len(mutation_position) > 1 and len(mutated_aa) >1:
#         if not all([i.isalpha() and len(i) == 1 and i != 'X' for i in mutated_aa]):
#             warnings.append('Exceptional mutation case')
#             return wildtype_peptides, mutated_peptides, warnings
            
#     ############################################## then make the sequence ################################################################# 
#     # if no sequence info from ithanet, then get it from hbvar based on first hgsv name and then hb name 
    
#     if pd.notna(check_protein_info):
#         check_deletion = re.findall(r'\b[\w-]+>0$', row['Protein Info'])
#     elif pd.notna(check_common_name):
#         check_deletion = re.findall(r'\b[\w-]+>0$', row['Common Name'])
#     else:
#         check_deletion = None
    
#     # when no sequence info and is not deletion 
#     if pd.isna(check_sequence) and not check_deletion:
        
#         # we first get sequence based on hgvs name 
#         if pd.notna(check_hgvs_name):
#             #sometimes the hgvs name are multiple e.g. HBA1:c.19_21delGAC | HBA2:c.19_21delGAC, so split and loop
#             row_name = row['HGVS Name'].split('|')
#             for name_idx in range(len(row_name)):
#                 #name_idx = 0
#                 row_name_idx = row_name[name_idx].strip()
#                 matched_seq = hbvar_id_seq_hgvs_hb[hbvar_id_seq_hgvs_hb['hgvs'] == row_name_idx]['sequence'].values
#                 if matched_seq.size > 0:
#                     sequence.extend(matched_seq)
                    
#         sequence = list(set(sequence)) #unique sequences in case, multiple hgsv name refer to the same sequence

#         # if we cannot find anything from hgvs name, then we use hb name
#         if not sequence and pd.notna(check_hb_name):
#             row_name = row['Hb Name'].strip()
#             matched_seq = hbvar_id_seq_hgvs_hb[hbvar_id_seq_hgvs_hb['hb_name'] == row_name]['sequence'].values
#             matched_seq = np.unique(matched_seq)
#             if matched_seq.size > 0:
#                 sequence.extend(matched_seq)
        
#     # means it is deletion, we just take the wildtype sequences
#     elif check_deletion:
#         # get gene from Protein Info
#         gene = re.search(r'^[αδβAγGγ]+\d?', row['Protein Info']).group()
#         #convert greek letters to english 
#         conversion = {'α1': 'HBA1', 'α2': 'HBA2', 'β': 'HBB', 'δ': 'HBD', 'Aγ': 'HBG1', 'Gγ': 'HBG2'}
#         try:
#             gene_english = conversion[gene]
#             sequence = [get_sequence_fasta(hemoglobin_genes_fasta, gene_english)] #list sequence to be consistent
#         except:
#             warnings.append('Gene name not in the list')
#             return wildtype_peptides, mutated_peptides, warnings
    
#     # if there is protein sequence info already, take it from ithanet         
#     else:
#         sequence = [row['Protein sequence'] ]
    

#     # make decision based on sequence
#     if len(sequence) > 1:
#         warnings.append('None unique sequence, likely 2 mutations on the same site')
#         return wildtype_peptides, mutated_peptides, warnings
#     elif len(sequence) == 1:
#         sequence = sequence[0]
#         # remove the leading M from the protein sequence
#         if sequence[0] == 'M':
#             sequence = sequence[1:]
#     elif not sequence:
#         warnings.append('No sequence')
#         return wildtype_peptides, mutated_peptides, warnings

#     ############################################## to get wildtype_peptides, mutated_peptides, warnings ################################################################# 
#     try:

#         # get mutated aa and compare if info in Protein Info/Common Name and Protein sequences are consistent, if not we use Protein Info/Common Name as standard
#         # when it is not deletion 
#         if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
#             for pos_idx in range(len(wildtype_mutated_aa[1])):
#                 #pos_idx = 0
#                 mutated_aa = get_wildtype_mutated_aa(string=mutation_info)[1][pos_idx] # this mutation aa is obtained from protein info
#                 if mutated_aa == 'K' or mutated_aa == 'R':
#                     warnings.append(f'mutated amino acid is K/R')
#                 if sequence[mutation_position[pos_idx]] != mutated_aa:
#                     #print(pos_idx)
#                     # when there are multiple mutations, ithanet normally only show one mutation in the sequence, so you always find inconsistency here
#                     warnings.append(f'mutation info not consistent between Protein Info and Protein sequence')
#                     # however, if only one mutation, but the mutation is not consistent, then it is an issue 
#                     if len(mutation_position) == 1:
#                         warnings.append(f'Attention!!!, mutation info not consistent between Protein Info and Protein sequence')
#                     sequence = sequence[:mutation_position[pos_idx]] + mutated_aa + sequence[mutation_position[pos_idx]+1:]
        
#         # Now the sequence is ready to be digested. Ideally the sequences contains all mutations, but does not work for multiple mutation at the same locus
#         # when it is not deletion
#         # get mutated_peptides
#         if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
#             # wildtype_mutated_aa[1] is mutated_aa, can be >1, so loop it 
#             for mut_idx in range(len(wildtype_mutated_aa[1])):    
#                 #mut_idx = 0
#                 miss_cleavage = 0
#                 tmp_mutated_peptide = find_mutated_peptide(sequence, mutation_position[mut_idx], enzyme='trypsin', missed_cleavages=miss_cleavage, regex=True)[0]
#                 while tmp_mutated_peptide and len(str(tmp_mutated_peptide)) < 7:
#                     miss_cleavage += 1
#                     tmp_mutated_peptide = find_mutated_peptide(sequence, mutation_position[mut_idx], missed_cleavages=miss_cleavage)
#                     # can be multiple peptides all containing mutated aa
#                     if isinstance(tmp_mutated_peptide, list) and len(tmp_mutated_peptide) > 1:
#                         tmp_mutated_peptide = select_peptide(peptides=tmp_mutated_peptide) # select peptide having the most mutation number and longest
#                 mutated_peptides.append(tmp_mutated_peptide) # can be multiple peptides that each has a mutation, or one peptide has multiple mutations
#             mutated_peptides = list(dict.fromkeys(mutated_peptides)) # unique lists if 2 peptides both having 2 mutations, so selected multiple times 
            
#         # else means to handle the deletion 
#         else:
#             deletions = []
#             for index, cont in enumerate(mutation_position):
#                 deletions.append((cont+1, wildtype_mutated_aa[0][index]))

#             miss_cleavage = 0
#             peptide_pair = find_mutated_peptide_deletion(sequence, deletions, missed_cleavages=miss_cleavage) # returns None, if deletion is K/R
            
#             # if we see one deletion is an entire peptide, then we skip 
#             if pd.notna(list(peptide_pair)[1]):
#                 if not list(peptide_pair)[1].startswith('-'):
#                     if any([pd.notna(i) for i in peptide_pair]):
#                         while any([len(i) <7 for i in peptide_pair]):
#                             miss_cleavage = miss_cleavage + 1
#                             print(miss_cleavage)
#                             peptide_pair = find_mutated_peptide_deletion(sequence, deletions, miss_cleavage) 

#             wildtype_peptides, mutated_peptides = peptide_pair
#             return [wildtype_peptides], [mutated_peptides], warnings
            

#         # check if mutated position is in same length as mutated_peptides
#         if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
#             if len(mutated_peptides) != len(mutation_position):
#                 warnings.append(f'mutated peptides number is different from mutation positions') # when it happens?
#                 return wildtype_peptides, mutated_peptides, warnings
            
#         # get wildtype_peptides
#         if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
#             for pep_idx in range(len(wildtype_mutated_aa[0])):
#                 #pep_idx = 0    
                
#                 peptide = mutated_peptides[pep_idx]
#                 mut_aa = wildtype_mutated_aa[1][pep_idx].lower()
#                 wild_aa = wildtype_mutated_aa[0][pep_idx]
#                 tmp_wildtype_peptide = re.sub(mut_aa, wild_aa, peptide)
#                 wildtype_peptides.append(tmp_wildtype_peptide)
#                 # this is to check if the mut_aa is inside the sequences and correctly replaced. if not it can be due to (1) mutation position is 1, but mutations is 2
#                 check = any([i.islower() for i in tmp_wildtype_peptide]) 
#                 if check:
#                     warnings.append(f'mutation aa is not in the sequence')

#                 if wild_aa == 'K' or wild_aa == 'R':
#                     warnings.append(f'wildtype amino acid is K/R')

#             return wildtype_peptides, mutated_peptides, warnings
            
#     except IndexError:
#         # this error is for all other weird errors
#         warnings.append(f'IndexError for ithanet id {id}')
#         return wildtype_peptides, mutated_peptides, warnings
    

def process_info_trypsin(row, hbvar_id_seq_hgvs_hb, hemoglobin_genes_fasta):
    # this version removing leading M in the proteins 
    
    id = row['IthaID']
    check_protein_info = row['Protein Info']
    check_sequence = row['Protein sequence']
    check_common_name = row['Common Name']
    check_hgvs_name = row['HGVS Name']
    check_hb_name = row['Hb Name']

    sequence = []
    mutated_peptides = []
    wildtype_peptides = []
    warnings = []

    ############################################## first make mutation position #################################################################
    # make mutation info, use Common Name first, then Protein Info, (as I found they can diverge e.g. ithanet 3651, CD position is more often correct maybe)
    if pd.notna(check_protein_info) and pd.isna(check_common_name):
        mutation_info = row['Protein Info'] 
        wildtype_mutated_aa = get_wildtype_mutated_aa(string=row['Protein Info']) # aa mutation info in Protein Info is better formatted
        mutation_position = get_position(mutation_info, pattern=r'\b(?:\D*)(\d+)\(')
        
    elif pd.notna(check_common_name) and pd.isna(check_protein_info):
        mutation_info = row['Common Name'] # make a variable to facilitate following steps
        wildtype_mutated_aa = get_wildtype_mutated_aa(string=row['Common Name'])
        mutation_position = get_position(row['Common Name'], pattern=r'CD\s+(\d+)\s+') 

        # when both info is present, mutated aa is more correct in protein info 
    elif pd.notna(check_protein_info) and pd.notna(check_common_name):
        mutation_info = row['Protein Info'] 
        wildtype_mutated_aa = get_wildtype_mutated_aa(string=row['Protein Info']) # aa mutation info in Protein Info is better formatted
        mutation_position_protein_info = get_position(mutation_info, pattern=r'\b(?:\D*)(\d+)\(')
        mutation_position_common_name = get_position(row['Common Name'], pattern=r'[Cc][Dd]\s+(\d+)\s+') 
        # when position info conficts, choose the correct info from common name
        mutation_position = choose_mutation_position(mutation_position_protein_info, mutation_position_common_name)

    # if no mutation info, then skip 
    if pd.isna(check_protein_info) and pd.isna(check_common_name):
        warnings.append('No mutation info')
        return wildtype_peptides, mutated_peptides, warnings  # Skip the rest of the loop for this iteration
        
    
    # if we see mutation_position is None, then we skip this iteration    
    if not isinstance(mutation_position, list):
        warnings.append('No mutation position info')
        return wildtype_peptides, mutated_peptides, warnings
        
    
    
    # if mutated aa is not point mutation, not multiple mutations, not deletion, but has other category as 'X', then skip
    if any([i == 'X' for i in wildtype_mutated_aa[1]]):
        warnings.append('Exceptional mutation case')
        return wildtype_peptides, mutated_peptides, warnings
        
        

    # if mutation position is multiple, and it is deletion, then skip 
    if len(mutation_position) > 1 and any([i == 'del' for i in wildtype_mutated_aa[1]]):
        warnings.append('Complex deletions')
        return wildtype_peptides, mutated_peptides, warnings
        

    # Convert mutation positions to integers 
    mutation_position = [int(pos) - 1 for pos in mutation_position] # minus position by 1, since we will remove the leading M from the sequence 

    # # test mutation condition 
    # 1, point mutation (accept)
    # 2, multiple point mutation (accept) 
    # 3, single deletion (accept)
    # 4, multiple deletion (no)
    # 5, mixture (no)
    
    # if multple position and mutations,and it is not just point mutations, then skip 
    mutated_aa = get_wildtype_mutated_aa(string=mutation_info)[1]
    if len(mutation_position) > 1 and len(mutated_aa) >1:
        if not all([i.isalpha() and len(i) == 1 and i != 'X' for i in mutated_aa]):
            warnings.append('Exceptional mutation case')
            return wildtype_peptides, mutated_peptides, warnings
            
            
    ############################################## then make the sequence ################################################################# 
    # if no sequence info from ithanet, then get it from hbvar based on first hgsv name and then hb name 
    
    if pd.notna(check_protein_info):
        check_deletion = re.findall(r'\b[\w-]+>0$', row['Protein Info'])
    elif pd.notna(check_common_name):
        check_deletion = re.findall(r'\b[\w-]+>0$', row['Common Name'])
    else:
        check_deletion = None
    
    # when no sequence info and is not deletion 
    if pd.isna(check_sequence) and not check_deletion:
        
        # we first get sequence based on hgvs name 
        if pd.notna(check_hgvs_name):
            #sometimes the hgvs name are multiple e.g. HBA1:c.19_21delGAC | HBA2:c.19_21delGAC, so split and loop
            row_name = row['HGVS Name'].split('|')
            for name_idx in range(len(row_name)):
                #name_idx = 0
                row_name_idx = row_name[name_idx].strip()
                matched_seq = hbvar_id_seq_hgvs_hb[hbvar_id_seq_hgvs_hb['hgvs'] == row_name_idx]['sequence'].values
                if matched_seq.size > 0:
                    sequence.extend(matched_seq)
                    
        sequence = list(set(sequence)) #unique sequences in case, multiple hgsv name refer to the same sequence

        # if we cannot find anything from hgvs name, then we use hb name
        if not sequence and pd.notna(check_hb_name):
            row_name = row['Hb Name'].strip()
            matched_seq = hbvar_id_seq_hgvs_hb[hbvar_id_seq_hgvs_hb['hb_name'] == row_name]['sequence'].values
            matched_seq = np.unique(matched_seq)
            if matched_seq.size > 0:
                sequence.extend(matched_seq)
        
    # means it is deletion, we just take the wildtype sequences
    elif check_deletion:
        # get gene from Protein Info
        gene = re.search(r'^[αδβAγGγ]+\d?', row['Protein Info']).group()
        #convert greek letters to english 
        conversion = {'α1': 'HBA1', 'α2': 'HBA2', 'β': 'HBB', 'δ': 'HBD', 'Aγ': 'HBG1', 'Gγ': 'HBG2'}
        try:
            gene_english = conversion[gene]
            sequence = [get_sequence_fasta(hemoglobin_genes_fasta, gene_english)] #list sequence to be consistent
        except:
            warnings.append('Gene name not in the list')
            return wildtype_peptides, mutated_peptides, warnings
            
    
    # if there is protein sequence info already, take it from ithanet         
    else:
        sequence = row['Protein sequence']
        sequence = sequence.replace('_x000D_','').replace('\n','') #I have seen some cases there are extract line or weird things, so clean them 
        sequence = [sequence]

    # make decision based on sequence
    if len(sequence) > 1:
        warnings.append('None unique sequence, likely 2 mutations on the same site')
        return wildtype_peptides, mutated_peptides, warnings
        
    elif len(sequence) == 1:
        sequence = sequence[0]
        # remove the leading M from the protein sequence
        if sequence[0] == 'M':
            sequence = sequence[1:]
    elif not sequence:
        warnings.append('No sequence')
        return wildtype_peptides, mutated_peptides, warnings
        

    ############################################## to get wildtype_peptides, mutated_peptides, warnings ################################################################# 
    try:

        # get mutated aa and compare if info in Protein Info/Common Name and Protein sequences are consistent, if not we use Protein Info/Common Name as standard
        # when it is not deletion 
        if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
            for pos_idx in range(len(wildtype_mutated_aa[1])):
                #pos_idx = 0
                mutated_aa = get_wildtype_mutated_aa(string=mutation_info)[1][pos_idx] # this mutation aa is obtained from protein info
                if mutated_aa == 'K' or mutated_aa == 'R':
                    warnings.append(f'mutated amino acid is K/R')
                if sequence[mutation_position[pos_idx]] != mutated_aa:
                    #print(pos_idx)
                    # when there are multiple mutations, ithanet normally only show one mutation in the sequence, so you always find inconsistency here
                    warnings.append(f'mutation info not consistent between Protein Info and Protein sequence')
                    # however, if only one mutation, but the mutation is not consistent, then it is an issue 
                    if len(mutation_position) == 1:
                        warnings.append(f'Attention!!!, mutation info not consistent between Protein Info and Protein sequence')
                    sequence = sequence[:mutation_position[pos_idx]] + mutated_aa + sequence[mutation_position[pos_idx]+1:]
        
        # Now the sequence is ready to be digested. Ideally the sequences contains all mutations, but does not work for multiple mutation at the same locus
        # when it is not deletion
        # get mutated_peptides
        if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
            # wildtype_mutated_aa[1] is mutated_aa, can be >1, so loop it 
            for mut_idx in range(len(wildtype_mutated_aa[1])):    
                #mut_idx = 0
                miss_cleavage = 0
                tmp_mutated_peptide = find_mutated_peptide(sequence, mutation_position[mut_idx], enzyme='trypsin', missed_cleavages=miss_cleavage, regex=True)[0]
                while tmp_mutated_peptide and len(str(tmp_mutated_peptide)) < 7:
                    miss_cleavage += 1
                    tmp_mutated_peptide = find_mutated_peptide(sequence, mutation_position[mut_idx], missed_cleavages=miss_cleavage)
                    # can be multiple peptides all containing mutated aa
                    if isinstance(tmp_mutated_peptide, list) and len(tmp_mutated_peptide) > 1:
                        tmp_mutated_peptide = select_peptide(peptides=tmp_mutated_peptide) # select peptide having the most mutation number and longest
                mutated_peptides.append(tmp_mutated_peptide) # can be multiple peptides that each has a mutation, or one peptide has multiple mutations
            mutated_peptides = list(dict.fromkeys(mutated_peptides)) # unique lists if 2 peptides both having 2 mutations, so selected multiple times 
            
        # else means to handle the deletion 
        else:
            deletions = []
            for index, cont in enumerate(mutation_position):
                deletions.append((cont+1, wildtype_mutated_aa[0][index]))

            miss_cleavage = 0
            peptide_pair = find_mutated_peptide_deletion(sequence, deletions, miss_cleavage) # returns None, if deletion is K/R
            
            # if we see one deletion is an entire peptide, then we skip 
            if pd.notna(peptide_pair[1]):
                if not peptide_pair[1][0].startswith('-'):
                    if any([pd.notna(i) for i in peptide_pair]):
                        while any([len(i[0]) < 7 for i in peptide_pair]):
                            miss_cleavage = miss_cleavage + 1
                            peptide_pair = find_mutated_peptide_deletion(sequence, deletions, miss_cleavage) 

            wildtype_peptides, mutated_peptides = peptide_pair
            return wildtype_peptides, mutated_peptides, warnings
            
        
            
        # check if mutated position is in same length as mutated_peptides
        if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
            if len(mutated_peptides) != len(mutation_position):
                warnings.append(f'mutated peptides number is different from mutation positions') # when it happens?
                return wildtype_peptides, mutated_peptides, warnings
            
        # get wildtype_peptides
        if not any([i == 'del' for i in wildtype_mutated_aa[1]]):
            for pep_idx in range(len(wildtype_mutated_aa[0])):
                #pep_idx = 0    
                
                peptide = mutated_peptides[pep_idx]
                mut_aa = wildtype_mutated_aa[1][pep_idx].lower()
                wild_aa = wildtype_mutated_aa[0][pep_idx]
                tmp_wildtype_peptide = re.sub(mut_aa, wild_aa, peptide)
                wildtype_peptides.append(tmp_wildtype_peptide)
                # this is to check if the mut_aa is inside the sequences and correctly replaced. if not it can be due to (1) mutation position is 1, but mutations is 2
                check = any([i.islower() for i in tmp_wildtype_peptide]) 
                if check:
                    warnings.append(f'mutation aa is not in the sequence')

                if wild_aa == 'K' or wild_aa == 'R':
                    warnings.append(f'wildtype amino acid is K/R')

            return wildtype_peptides, mutated_peptides, warnings
            
    except IndexError:
        # this error is for all other weird errors
        warnings.append(f'IndexError for ithanet id {id}')
        return wildtype_peptides, mutated_peptides, warnings
        