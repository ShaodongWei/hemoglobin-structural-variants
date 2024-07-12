import pandas as pd
import numpy as np
from tqdm import tqdm
import os

# Add path to scripts 
new_path = os.path.join(os.getcwd(), 'scripts')
sys.path.append(new_path)

import functions  # Load functions
pd.set_option('display.max_columns', None)

##############################################################
### Make hbvar data, which has id, hgvs, hb name, sequence ###

# Read in the hbvar JSON data, downloaded from hbvar
import json
with open('./data_in/hbvar_full.json', 'r') as file:
    js = json.load(file)

# Flatten out the JSON data 
js_df = pd.DataFrame()
for i in range(len(js)):
    sbs = js[i]
    sbs = pd.DataFrame(functions.flatten_dict_hbvar(sbs), columns=['path', 'value'])
    id = sbs[sbs['path'].str.contains('^id$')]['value'].values.astype('str')
    hgvs = sbs[sbs['path'].str.contains('^hgvsName$')]['value'].values 
    hb_name = sbs[sbs['path'].str.contains('^name$')]['value'].values
    tmp = pd.DataFrame({'id': id, 'hgvs': hgvs, 'hb_name': hb_name})
    js_df = pd.concat([js_df, tmp], axis=0)

js_df.to_pickle('./data_out/hbvar_id_hgvs_hb.pkl')

# Read in the protein sequences data, downloaded from hbvar protein sequences for all entries
from Bio import SeqIO

list_fasta = []
for record in SeqIO.parse('./data_in/hbvar_full.fasta', 'fasta'):
    header = record.id
    sequence = str(record.seq)
    list_fasta.append({'header': header, 'sequence': sequence})

df_fasta = pd.DataFrame(list_fasta)

def get_hbvarid(row):
    row_split = row['header'].split(':')[0]
    hbvarid = row_split.split('_')[-1]
    return hbvarid
    
df_fasta['id'] = df_fasta.apply(get_hbvarid, axis=1).astype('str')

# Merge 
hbvar_id_seq_hgvs_hb = pd.merge(js_df, df_fasta, on='id', how='right')  # Only keep those having protein sequences, on right 
hbvar_id_seq_hgvs_hb.to_pickle('./data_out/hbvar_id_seq_hgvs_hb.pkl')

#########################################################################################################
### Make inthanet table, which has IthaID, Common Name, Hb Name, HGVS Name, Genes, Pathogenicity, etc ###

itha = pd.read_csv('./data_in/IthaGenes_structual_Hb.csv')  # Downloaded for all structural variants 
itha

# Above data misses protein sequence and mutation info, now we make it 
ids = itha['IthaID'].values
info = ['Protein sequence', 'Protein Info']

# Now we webscraping info from the webpage of ithanet 
ot = pd.DataFrame()
for i in tqdm(ids):
    tmp = functions.get_ithanet_info(i, info)
    if tmp is not None:
        ot = pd.concat([ot, tmp], axis=0)

itha_prot = pd.merge(itha, ot, left_on='IthaID', right_on='id', how='left')
itha_prot
# I manually checked the protein sequences, and found that ithanet ID 1355, 1373 the protein sequences were wrong, I replace them with the sequences from HbVar 
itha_prot.to_excel('./data_out/itha_protein.xlsx', index=False, header=True)

#####################################################################################################
### Webscraping and text processing to get the unique proteomics peptides for hemoglobin variants ###

hbvar_id_seq_hgvs_hb = pd.read_pickle('./data_out/hbvar_id_seq_hgvs_hb.pkl')
hemoglobin_genes_fasta = './data_in/hemoglobin.genes.fasta'  # Wildtype hemoglobin protein sequences, used for deletion variant 

itha_prot = pd.read_excel('./data_out/itha_protein.xlsx')  # This is from webscraping
itha_prot['trypsin'] = ''
# Re-order columns 
itha_prot.drop(['id'], axis=1, inplace=True)
order = ['IthaID', 'Common Name', 'Hb Name', 'HGVS Name', 'Protein Info', 'trypsin', 'Protein sequence', 'Genes',
       'Functionality', 'Phenotype', 'Locus', 'Position']
itha_prot = itha_prot[order]
itha_prot

itha_prot['trypsin'] = itha_prot.apply(functions.process_info_trypsin, axis=1, args=(hbvar_id_seq_hgvs_hb, hemoglobin_genes_fasta))
#itha_prot.to_excel('ithanet_protein_trypsin2.xlsx', index=False)

# Rename Hb Name if empty, Hb Name is going to be used in the final fasta annotation file. 
itha_prot.loc[pd.isna(itha_prot['Hb Name']), 'Hb Name'] = 'IthaNet-' + itha_prot.loc[pd.isna(itha_prot['Hb Name']), 'IthaID'].astype('str')

###################################################################################################
### Check duplicated peptides for (1) single mutation, (2) double mutation, (3) only duplicates ###

## (1) Single mutation
df = itha_prot[['IthaID', 'trypsin', 'Genes', 'Protein Info', 'Common Name', 'HGVS Name', 'Hb Name']].copy()

df.loc[:, 'wt1'] = df['trypsin'].apply(lambda x: x[0][0] if x[0] else None)
df.loc[:, 'wt2'] = df['trypsin'].apply(lambda x: x[0][1] if (x[0] and len(x[0]) > 1) else None)

df.loc[:, 'mt1'] = df['trypsin'].apply(lambda x: x[1][0] if x[0] else None)
df.loc[:, 'mt2'] = df['trypsin'].apply(lambda x: x[1][1] if x[0] and len(x[0]) > 1 else None)

df = df[pd.notna(df['wt1']) & pd.notna(df['mt1'])]

# Make single mutation df 
df_sg = df[pd.notna(df['wt1']) & pd.isna(df['wt2'])].copy()

# Unique mutation peptide

df_sg['mt1'] = df_sg['mt1'].str.upper()  # Need to upper before hand 
df_sg = df_sg[~df_sg['mt1'].duplicated(keep=False)]
df_sg = df_sg.drop(['wt2', 'mt2'], axis=1)

df['Protein Info'].str.contains('del', case=False).sum()  # Check how many are deletions 

# Compare if any of wt1 to any of mt1 are duplicated
ot = []
for i in df_sg['wt1']:
    for j in df_sg['mt1']:
       if i == j:
           ot.append(j)
            
seq_dup = pd.Series(ot).unique()  # There are 6 unique mt1 sequences that are not usable as they are same as wildtype 
df_sg.shape
df_sg = df_sg[~df_sg['mt1'].isin(seq_dup)]
df_sg.shape  # These 6 are removed 

# Compare if any of mt1 is same as that in the human proteome
import Bio
from Bio import SeqIO

prot = []
for line in SeqIO.parse('./data_in/UP000005640_9606.fasta', 'fasta'):
    header = line.id
    seq = line.seq
    prot.append({'header': header, 'seq': str(seq)})

prot = pd.DataFrame(prot)
prot

# Trypsin digest human proteome 
prot['trypsin'] = prot['seq'].apply(lambda x: functions.cleave_sequence_insensitive(x, missed_cleavages=0))
prot 

# Check each row of df_sg mt1 if it is inside of human proteome 
human_prot_pep = []
for i in prot['trypsin'].tolist():
    for j in i:
        human_prot_pep.append(j)

# Unique peptides
human_prot_pep_set = list(set(human_prot_pep))

def check_duplicate(query, reference):
    return query.strip().upper() in reference 

df_sg['check'] = df_sg['mt1'].apply(lambda x: check_duplicate(query=x, reference=human_prot_pep_set))
df_sg['check'].sum()  # Only 1 duplicate

# Check what is the duplicate
df_sg[df_sg['check']]
the_one = df_sg[df_sg['check']]['mt1'].tolist() 
for i in the_one:
    print(prot[prot['trypsin'].apply(lambda x: i in x)]) 

# Remove duplicates from df_sg
df_sg.shape
df_sg = df_sg[~df_sg['mt1'].isin(the_one)]
df_sg.shape  # Now the final version, 1009 variants 

# Check NaNs
df_sg[pd.isna(df_sg['mt1'])]  # No NaNs 
df_sg[pd.isna(df_sg['Hb Name'])]  # No NaNs 

df_sg  # Ready to use (single mutation peptide)

## (2) Double mutation peptides
df

# Multiple mutation 
df['trypsin'].apply(lambda x: len(x[1]) if x[1] else 0).unique()  # Multiple mutation is all double mutation

df_db = df[pd.notna(df['mt2'])].copy()  # Double mutation
df_db.shape

# To uppercase
df_db.loc[:, 'mt1'] = df_db['mt1'].str.strip().str.upper()
df_db.loc[:, 'mt2'] = df_db['mt2'].str.strip().str.upper()

# Check mutation peptides duplicates 
df_db.loc[:, 'merge'] = df_db['mt1'] + df_db['mt2']  # Concatenate peptides
# Convert to upper case before check duplicates
df_db['merge'].duplicated(keep=False).sum()  # No duplicates for double mutation 

# Check against human proteome 
prot  # This is the human proteome digested (see it above)
def check_duplicate(query, reference):
    return query.strip().upper() in reference 

df_db['check1'] = df_db['mt1'].apply(lambda x: check_duplicate(query=x, reference=human_prot_pep_set))
df_db['check2'] = df_db['mt2'].apply(lambda x: check_duplicate(query=x, reference=human_prot_pep_set))
df_db['check_merge'] = df_db['merge'].apply(lambda x: check_duplicate(query=x, reference=human_prot_pep_set))
df_db['check1'].sum()  # 0 duplicates
df_db['check2'].sum()  # 1 duplicate
df_db['check_merge'].sum()  # 0 duplicates

df_db.shape  # 31 variants

## (3) Only the duplicated peptides
df_dup = df.copy()

# To uppercase
df_dup['mt1'] = df_dup['mt1'].str.strip().str.upper()
df_dup['mt2'] = df_dup['mt2'].str.strip().str.upper()

df_dup = df_dup[df_dup[['mt1', 'mt2']].duplicated(keep=False)]
df_dup.sort_values('mt1', inplace=True)
df_dup  # 259 rows

# Combine mt1 and mt2, then check duplicates, in this way it considers both single/double mutation peptides
def concat(row):
    if pd.notna(row['mt2']):
        ot = row['mt1'] + row['mt2']
    else:
        ot = row['mt1']
    return ot

df_dup['merge'] = df_dup.apply(concat, axis=1)
check_set = set(df_dup['merge'])

# Combine duplicated rows 
def unique_row(series):
    ot = series.dropna().unique()
    return '| '.join(ot)

df_dup2 = pd.DataFrame()
for i in check_set:
    sbs = df_dup[df_dup['merge'] == i].copy()
    sbs['IthaID'] = sbs['IthaID'].astype('str')  # In order to concatenate, convert int to str 
    sbs.drop('trypsin', inplace=True, axis=1)  # Unique does not work for lists, so drop it
    tmp = sbs.apply(unique_row, axis=0).to_frame().T
    if (tmp['Hb Name'] == '').values:
        print(i)
    df_dup2 = pd.concat([df_dup2, tmp], axis=0)

df_dup2  # 125 rows

# Check if duplicates happen for double mutation peptides
np.sum(pd.isna(df_dup2['mt2']))  # No one 

# Compare if any of wt1 to any of mt1
ot = []
for i in df_dup2['wt1']:
    for j in df_dup2['mt1']:
       if i == j:
           ot.append(j)         
ot  # No duplicates

# Check against human proteome 
prot  # This is the human proteome digested (see it above)
def check_duplicate(query, reference):
    return query.strip().upper() in reference 

check = df_dup2['mt1'].apply(lambda x: check_duplicate(query=x, reference=human_prot_pep_set))
np.sum(check)  # 2 duplicates between human proteome 

# Check what are the duplicates
df_dup2[check]
the_ones = df_dup2[check]['mt1'].tolist() 
for i in the_ones:
    print(prot[prot['trypsin'].apply(lambda x: i in x)])  # 3 duplicates in human proteome 

# Remove them from df_dup2
df_dup2.shape
df_dup2 = df_dup2[~df_dup2['mt1'].isin(the_ones)]
df_dup2.shape  # 123 rows

# Rename df_dup2
df_dup = df_dup2

#######################################################################
### Check duplicates between single, double, and duplicated 3 files ###
sg = df_sg.copy()
sg.shape
db = df_db.copy()
db['check'] = db['mt1'] + db['mt2']
db.shape
dup = df_dup.copy()
dup.shape

# Check duplicates between 3 files 
set(sg['mt1']) & set(db['check'])  # No duplicates
set(sg['mt1']) & set(dup['mt1'])  # No duplicates
set(db['check']) & set(dup['mt1'])  # No duplicates, when check mt we use mt1 + mt2

# Check mt against wt 
pooled_wt = sg['wt1'].to_list() + db['wt1'].to_list() + db['wt2'].to_list() + dup['wt1'].to_list()  # When check wildtype, use wt1 and wt2 separately
pooled_wt = list(set(pooled_wt))
pooled_wt = [i for i in pooled_wt if i != '']
len(pooled_wt)

sg_rm = [i for i in sg['mt1'] if i in pooled_wt]
len(sg_rm)  # 8
len([i for i in db['mt1'] if i in pooled_wt])  # 0
db_rm = [i for i in db['mt2'] if i in pooled_wt]
len(db_rm)  # 1
dup_rm = [i for i in dup['mt1'] if i in pooled_wt]
len(dup_rm)  # 0

# Remove duplicates from single/double/duplicates peptides 
df_sg = df_sg[~df_sg['mt1'].isin(sg_rm)]
df_sg.shape  # 1001 rows

df_db = df_db[~df_db['mt2'].isin(db_rm)]
df_db.shape  # 30 rows

###################
### Save tables ###

df_sg[['Hb Name', 'wt1', 'mt1']].to_csv('./data_out/single.mutation.csv', index=False, sep=';')
df_sg.drop(columns=['check']).to_excel('./data_out/single.mutation.xlsx', index=False)

df_db.drop(columns=['check1', 'check2', 'check_merge']).to_excel('./data_out/double.mutation.xlsx', index=False)
df_db[['Hb Name', 'wt1', 'wt2', 'mt1', 'mt2']].to_csv('./data_out/double.mutation.csv', index=False, sep=';')

df_dup.drop(columns=['merge']).to_excel('./data_out/duplicates.xlsx', index=False)
df_dup[['Hb Name', 'wt1', 'mt1']].to_csv('./data_out/duplicates.csv', sep=';', index=False, quoting=False)

####################################################
### Prepare the wildtype peptide for HBB and HBD ###
from Bio import SeqIO
hemoglobin_genes_fasta = './data_in/hemoglobin.genes.fasta'  # Wildtype hemoglobin protein sequences, used for deletion variant 

hemo = []
for record in SeqIO.parse(hemoglobin_genes_fasta, 'fasta'):
    name = record.id
    seq = str(record.seq)
    tmp = {'gene': name, 'seq': seq}
    hemo.append(tmp)

hemo = pd.DataFrame(hemo)
hemo = hemo[hemo['gene'].apply(lambda x: 'HBB' in x or 'HBD' in x)]
hemo['trypsin'] = hemo['seq'].apply(lambda x: functions.cleave_sequence_insensitive(sequence=x, enzyme='trypsin', missed_cleavages=0))
hemo

# First map against our generated mt peptide file 
pooled_mt = df_sg['mt1'].tolist() + (df_db['mt1'] + df_db['mt2']).tolist() + (df_db['mt2'] + df_db['mt1']).tolist() + df_dup['mt1'].tolist()
pooled_mt = list(set(pooled_mt))
len(pooled_mt)  # 1178

# Function to exclude peptide from query that are found in reference
def get_unique_value(query_list, reference_list, invert=False):
    # By default invert is False, means to exclude values that are from reference.
    # When invert is True, means to get unique values from reference. 
    ot = []
    for i in query_list:
        if not invert:
            if i not in reference_list:
                ot.append(i)
        else:
            if i in reference_list:
                ot.append(i)
    return ot

# Exclude peptides that are from pooled mutation peptide 
hemo['unique_mt'] = hemo['trypsin'].apply(lambda x: get_unique_value(query_list=x, reference_list=pooled_mt))

# Next, exclude peptides that are from wildtype. Specifically e.g. for HBB, we exclude peptides that can refer to other genes from the human HBB trypsin cleavage. Same for HBD

# Make a dataframe that contains all wt from all genes 
db_wt2 = df_db[['HGVS Name', 'wt2']].copy()
db_wt2 = db_wt2.rename(columns={'wt2': 'wt1'})
df_wt = pd.concat([df_sg[['HGVS Name', 'wt1']], df_db[['HGVS Name', 'wt1']], db_wt2, df_dup[['HGVS Name', 'wt1']]], axis=0)
pooled_wt = list(set(df_wt['wt1'].values))
len(pooled_wt)

# Find out hbb wt peptide 
wt_hbb = (sg[sg['HGVS Name'].str.contains('HBB')]['wt1'].tolist() + 
db[db['HGVS Name'].str.contains('HBB')]['wt1'].tolist() + 
db[db['HGVS Name'].str.contains('HBB')]['wt2'].tolist() + 
dup[dup['HGVS Name'].str.contains('HBB', na=False)]['wt1'].tolist()
)
wt_hbb = list(set(wt_hbb))
wt_hbb = [i for i in wt_hbb if pd.notna(i)]

# Find out what peptides are shared with other genes
hbb_overlap = []
for i in wt_hbb:
    tmp = df_wt[df_wt['wt1'] == i]
    if np.sum(~tmp['HGVS Name'].str.contains('HBB')) > 0:
        hbb_overlap.append(i)
        
hbb_overlap  # 4 peptides

# Find out hbd wt peptide 
wt_hbd = (sg[sg['HGVS Name'].str.contains('HBD')]['wt1'].tolist() + 
db[db['HGVS Name'].str.contains('HBD')]['wt1'].tolist() + 
db[db['HGVS Name'].str.contains('HBD')]['wt2'].tolist() + 
dup[dup['HGVS Name'].str.contains('HBD', na=False)]['wt1'].tolist()
)
wt_hbd = list(set(wt_hbd))
wt_hbd = [i for i in wt_hbd if pd.notna(i)]

# Find out what peptides are shared with other genes
hbd_overlap = []
for i in wt_hbd:
    tmp = df_wt[df_wt['wt1'] == i]
    if np.sum(~tmp['HGVS Name'].str.contains('HBD')) > 0:
        hbd_overlap.append(i)
        
hbd_overlap  # 9 peptides

# Exclude peptide that are from overlap 
hemo['unique_mt_wt'] = ''
row_index_hbb = hemo.index[hemo['gene'].str.contains('HBB')].to_list()[0]
hemo.at[row_index_hbb, 'unique_mt_wt'] = get_unique_value(
    query_list=hemo.loc[hemo['gene'].str.contains('HBB'), 'unique_mt'].to_list()[0], 
    reference_list=hbb_overlap
)
row_index_hbd = hemo.index[hemo['gene'].str.contains('HBD')].to_list()[0]
hemo.at[row_index_hbd, 'unique_mt_wt'] = get_unique_value(
    query_list=hemo.loc[hemo['gene'].str.contains('HBD'), 'unique_mt'].to_list()[0], 
    reference_list=hbd_overlap
)
hemo['unique_mt'].apply(lambda x: len(x))
hemo['unique_mt_wt'].apply(lambda x: len(x))

# Next excluding peptide coming from entire human proteome 
# Make a peptide set excluding HBB and HBD separately
human_no_hbb = prot[~prot['header'].str.contains(r'\|HBB_HUMAN')]['trypsin'].tolist()
human_no_hbd = prot[~prot['header'].str.contains(r'\|HBD_HUMAN')]['trypsin'].tolist()

def flatten_list(list):
    ot = []
    for i in list:
        for j in i:
            ot.append(j)
    return ot

len(human_no_hbb)
human_no_hbb_set = list(set(flatten_list(human_no_hbb)))
human_no_hbd_set = list(set(flatten_list(human_no_hbd)))
len(human_no_hbb_set)

# Check duplicates against human proteome 
hemo.loc[hemo['gene'].str.contains('HBB'), 'unique_mt_wt_human'] = hemo.loc[hemo['gene'].str.contains('HBB'), 'unique_mt_wt'].apply(lambda x: get_unique_value(query_list=x, reference_list=human_no_hbb_set))
hemo.loc[hemo['gene'].str.contains('HBD'), 'unique_mt_wt_human'] = hemo.loc[hemo['gene'].str.contains('HBD'), 'unique_mt_wt'].apply(lambda x: get_unique_value(query_list=x, reference_list=human_no_hbd_set))

# Sort by length 
hemo['unique_mt_wt_human'] = hemo['unique_mt_wt_human'].apply(lambda x: sorted(x, key=len, reverse=True))
hemo['unique_mt_wt_human'].values  # It is the peptide that represents either HBB or HBD, not duplicates of anything else 

# Keep all unique peptides that can represent wildtype HBB, HBD
set(hemo.loc[2, 'unique_mt_wt_human']) & set(hemo.loc[3, 'unique_mt_wt_human'])  # No duplicates between them 
hemo.loc[2, 'unique_mt_wt_human'] 
hemo.loc[3, 'unique_mt_wt_human'] 

# Check how many times these peptides are in wt
[df_wt[df_wt['wt1'] == i].shape[0] for i in hemo.loc[2, 'unique_mt_wt_human']]  # HBB
[np.sum(~df_wt[df_wt['wt1'] == i]['HGVS Name'].str.contains('HBB')) for i in hemo.loc[2, 'unique_mt_wt_human']]  # No other gene, except HBB
[df_wt[df_wt['wt1'] == i].shape[0] for i in hemo.loc[3, 'unique_mt_wt_human']]  # HBD
[np.sum(~df_wt[df_wt['wt1'] == i]['HGVS Name'].str.contains('HBD')) for i in hemo.loc[3, 'unique_mt_wt_human']]  # No other gene, except HBD

# Keep all these peptides as for wildtype peptide for HBB, HBD
def make_table(row):
    ot = []
    pep = row['unique_mt_wt_human']
    gene = row['gene']
    for idx, value in enumerate(pep):
        header = gene + '_' + 'wt' + str(idx + 1)
        ot_tmp = {'gene': header, 'seq': value}
        ot.append(ot_tmp)
    return ot
    
table = pd.DataFrame(hemo.apply(lambda x: make_table(x), axis=1).sum())
table

table.to_csv('./data_out/hbb.hbd.csv', sep='~', index=False, quoting=False)
hemo.to_excel('./data_out/hbb.hbd.xlsx', index=False)

#################################
### Concatenate 4 fasta files ###

# This part should be done in python 

# Save fasta
awk -F';' 'NR>1 {print ">"$1"_wt""\n"$2; print ">"$1"_mt""\n"$3}' ./data_out/single.mutation.csv > ./data_out/single.mutation.fasta
awk -F';' ' NR > 1 {print ">"$1"_wt1""\n"$2;print ">"$1"_mt1""\n"$4;print ">"$1"_wt2""\n"$3;print ">"$1"_mt2""\n"$5;}'  ./data_out/double.mutation.csv > ./data_out/double.mutation.fasta
awk -F';' 'NR>1 {gsub(/"/,"",$1); print ">"$1"_wt""\n"$2; print">"$1"_mt""\n"$3}' ./data_out/duplicates.csv > ./data_out/duplicates.fasta
awk -F'~' 'NR > 1 {print ">"$1"\n"$2}' ./data_out/hbb.hbd.csv > ./data_out/hbb.hbd.fasta

# Merge fasta 
cat ./data_out/single.mutation.fasta ./data_out/double.mutation.fasta ./data_out/duplicates.fasta ./data_out/hbb.hbd.fasta > ./data_out/structural.variants.fasta  # In total 1164 variants

#####################################
### Unique all wildtype sequences ###

# After this wt is all unique, but mt is not, since mt can be same as one of double mutation mt
from Bio import SeqIO
import pandas as pd 

dat = []
for record in SeqIO.parse('./data_out/structural.variants.fasta', 'fasta'):
    name = record.description
    seq = str(record.seq)
    tmp = {'name': name, 'seq': seq}
    dat.append(tmp)

dat = pd.DataFrame(dat)

# Here duplicates contain both wt and mt, because mt can be same as one of mt1 or mt2, but not for mt1 + mt2
wt_mt_dup = dat[dat['seq'].duplicated(keep=False)]

# Make peptide sequences set 
wt_mt_seq = list(set(wt_mt_dup['seq'].values))
len(wt_mt_seq)

# Check 
wt_mt_uniq = pd.DataFrame()
for i in wt_mt_seq:
    sbs = wt_mt_dup[wt_mt_dup['seq'] == i].copy()
    sbs['name'] = sbs['name'].str.strip().str.replace('Hb ', '')
    sbs['name'] = sbs['name'].str.cat(sep='~')
    sbs = sbs.drop_duplicates()
    wt_mt_uniq = pd.concat([wt_mt_uniq, sbs], axis=0)
    
# Check duplicates for mt 
[i for i in wt_mt_uniq['name'].values if '_mt;' in i]

# Check seq duplicates
wt_mt_uniq['seq'].duplicated(keep=False).sum()  # No more duplicates

# Make a full data 
full_fasta = pd.DataFrame()
for i in dat['seq'].to_list():
    if i not in wt_mt_seq:
        sbs = dat[dat['seq'] == i]
    else:
        sbs = wt_mt_uniq[wt_mt_uniq['seq'] == i]
    full_fasta = pd.concat([full_fasta, sbs], axis=0)

full_fasta = full_fasta.drop_duplicates()
full_fasta  # 1333
dat  # 2380

# Check uniqueness for entire fasta
full_fasta['seq'].duplicated(keep=False).sum()  # No duplicates

# Write to a fasta file 
def dataframe_to_fasta(df, id_col, seq_col, output_file):
    with open(output_file, 'w') as file:
        for idx, row in df.iterrows():
            file.write(f">{row[id_col]}\n")
            file.write(f"{row[seq_col]}\n")

dataframe_to_fasta(full_fasta, 'name', 'seq', output_file='./data_out/structural.variants.fasta')

# Make a fasta where all ids are replaced with integers 
full_fasta['name2'] = [i for i in range(1, full_fasta.shape[0]+1)]
full_fasta  
dataframe_to_fasta(full_fasta, 'name2', 'seq', output_file='./data_out/structural.variants.2.fasta')

# Save data
full_fasta.to_excel('./data_out/full_fasta.xlsx', index=False)
