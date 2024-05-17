import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from thefuzz import process
from thefuzz import fuzz
import numpy as np
import re

filename = "230315_SAG211_8k_annotated_proteome_combined (1).fasta"

# code chunk adapted from David Cain (stackoverflow)
with open(filename) as fasta_file:  # Will close handle cleanly
    identifiers = []    # lists for each output
    descriptions = []
    lengths = []
    sequences = []
    isoelectrics = []
    weights = []
    for title, sequence in SimpleFastaParser(fasta_file):
        identifiers.append(title.split(None, 1)[0])     # First word is ID
        descriptions.append(title.split(None, 1)[1])    # Second word is description
        sequences.append(sequence.replace("*",""))      # Replace asterisks (STOP) for nothing
        x = ProteinAnalysis(sequence.replace("*",""))   # ProteinAnalysis object
        isoelectrics.append(x.isoelectric_point())      # Calculate isoelectric point
        weights.append(x.molecular_weight()/1000)       # Calculate molecular weight, divide by 1000 for kDa
        lengths.append(len(sequence.replace("*","")))   # Caclulate length
    print(str(len(identifiers))+str(" sequences parsed successfully."))     # Terminal output on success

# zip them all together in pandas dataframe
df = pd.DataFrame(zip(identifiers, descriptions, lengths, isoelectrics, weights, sequences), columns=['identifier', 'description', 'length', 'isoelectric','weight (kDa)','sequence'])

chlamy_df = pd.read_csv('chlamy_hits.txt', sep='\t', header=None)
chlamy_df.columns = ['identifier', 'chalmydomonas_hit', 'hit_descriptions']
df2 = pd.merge(df, chlamy_df, on='identifier', how='left')

motif = '[A|P][\D][Q][\D][\D][F][L|M|V][\D][R|K][\D]{0,2}[R|K]'
#motif = ('[A][E][Q][R][E][F][M|L][E][R][K][A]')

motifs = []
number_motifs = []
for entry in df2['sequence']:
    result = re.findall(motif, entry)
    motifs.append(result)
    number_motifs.append(len(result))

df2['number_motifs'] = number_motifs
df2['motifs'] = motifs

targetp_df = pd.read_csv('targetp.tsv', sep = '\t')
targetp_df.rename(columns={'ID': 'identifier'}, inplace=True)   # Rename the identifier column for merge purposes later
df3 = pd.merge(df2, targetp_df, on = 'identifier', how = 'left').fillna('0')       # Merge the dfs, using the identifier column, fill IDs with no value in prot_deg with dashes

tmhmm_df = pd.read_csv('tmhmm.tsv', sep = '\t')
tmhmm_df.rename(columns={'ID': 'identifier'}, inplace=True)   # Rename the identifier column for merge purposes later

lengths = []
for entry in tmhmm_df['Length']:
    length = str(entry.split('len=',1)[1])
    lengths.append(int(length))

numbers = []
for entry in tmhmm_df['Expected Aas']:
    number = str(entry.split('ExpAA=',1)[1])
    numbers.append(int(float(number)))

tmhmm_df['lengths'] = lengths
tmhmm_df['numbers'] = numbers
tmhmm_df['fraction_tm'] = tmhmm_df['numbers'].div(tmhmm_df['lengths'])

transmembrane_likelihood = []
for entry in tmhmm_df['fraction_tm']:
    if entry > 0.1:
        transmembrane_likelihood.append("yes, " + str(round(entry,2)*100)+"%")
    else:
        transmembrane_likelihood.append("no, " + str(round(entry,2)*100)+"%")
tmhmm_df['tm_likelihood'] = transmembrane_likelihood

transmembrane=[]
for entry in tmhmm_df['PredHel']:
    if entry == 'PredHel=0':
        transmembrane.append('no')
    else:
        transmembrane.append("yes, "+str(entry.split('PredHel=',1)[1]))
tmhmm_df['transmembranes'] = transmembrane

df4 = pd.merge(df3, tmhmm_df, on = 'identifier', how = 'left').fillna('0')       # Merge the dfs, using the identifier column, fill IDs with no value in prot_deg with dashes

# this dataset was contaminated with some chlamydomonas proteins, this code chunks removes them
with open('data.csv') as in_file, open('no_chlamy_data.csv', 'w') as out_file:
    for line in in_file:
        line = re.sub('([C][r][e]).+[_][j]', 'j', line)     # Find the Cre gene identifiers, up to the join with Chlorella IDs, remove it
        if not line.startswith('Cre'):                                  # If the line doesn't start with Cre, keep it, otherwise its a Chlamy protein so get rid
            out_file.write(line)

# read in whole cell proteomics data under different conditions
prot_deg_df = pd.read_csv('no_chlamy_data.csv')
prot_deg_df.rename(columns={'Protein.Group': 'prot_deg_id'}, inplace=True)   # Rename the identifier column for merge purposes later

# code chunk adapted from Corralien (stack overflow)
matched_id = lambda x: process.extractOne(x, df['identifier'], scorer=fuzz.partial_ratio)[2]  # here the [2] is the index from the extractOne function, used for mapping, partial_ratio ensures only complete matches are used
prot_deg_df['identifier'] = df.loc[prot_deg_df["prot_deg_id"].map(matched_id).values, 'identifier'].values  # takes the value from the prot_deg_id column, and matches against the identifier column df, then adds that column to df

df5 = pd.merge(df4, prot_deg_df, on = 'identifier', how = 'left').fillna('0')       # Merge the dfs, using the identifier column, fill IDs with no value in prot_deg with dashes

# read in RNAseq DEG data
rna_deg_df = pd.read_csv('data.tsv', sep ='\t')
rna_deg_df.rename(columns={'Name': 'identifier'}, inplace=True)   # Rename the identifier column for merge purposes later
rna_deg_df[['identifier','description']] = rna_deg_df["identifier"].str.split(" ", n=1, expand=True)
rna_deg_df.rename(columns={'Column1': 'assignment'}, inplace=True)   # Rename the identifier column for merge purposes later


df6 = pd.merge(df5, rna_deg_df, on = 'identifier', how = 'left').fillna('0') # Merge the dfs, using the identifier column, fill IDs with no value in prot_deg with dashes

coip_df = pd.read_csv('coip_data.tsv', sep ='\t')
coip_df.rename(columns={'First Protein ID': 'coip_identifier'}, inplace=True)   # Rename the identifier column for merge purposes later

# code chunk adapted from Corralien (stack overflow)
matched_id = lambda x: process.extractOne(x, df['identifier'], scorer=fuzz.partial_ratio)[2]  # here the [2] is the index from the extractOne function, used for mapping, partial_ratio ensures only complete matches are used
coip_df['identifier'] = df.loc[coip_df["coip_identifier"].map(matched_id).values, 'identifier'].values  # takes the value from the prot_deg_id column, and matches against the identifier column df, then adds that column to df

df7 = pd.merge(df6, coip_df, on = 'identifier', how = 'left').fillna('0') # Merge the dfs, using the identifier column, fill IDs with no value in prot_deg with dashes

pyrenoid_df = pd.read_csv('pyrenoid.csv')
pyrenoid_df.rename(columns={'Protein.Group': 'pyrenoid_deg_id'}, inplace=True)   # Rename the identifier column for merge purposes later
matched_id = lambda x: process.extractOne(x, df['identifier'], scorer=fuzz.partial_ratio)[2]  # here the [2] is the index from the extractOne function, used for mapping, partial_ratio ensures only complete matches are used
pyrenoid_df['identifier'] = df.loc[pyrenoid_df["pyrenoid_deg_id"].map(matched_id).values, 'identifier'].values  # takes the value from the prot_deg_id column, and matches against the identifier column df, then adds that column to df

df8 = pd.merge(df7, pyrenoid_df, on = 'identifier', how = 'left').fillna('0') # Merge the dfs, using the identifier column, fill IDs with no value in prot_deg with dashes

E065_df = pd.read_csv('E065.csv')
E065_df.rename(columns={'Protein.Group': 'E065_id'}, inplace=True)
#E065_df = E065_df.convert_objects(convert_numeric=True)
E065_df.replace('#NUM!', 0)

E065_df = pd.read_csv('E065.csv')
E065_df.rename(columns={'Protein.Group': 'E065_id'}, inplace=True)
E065_df = E065_df._convert(numeric=True)
E065_df = E065_df.replace(np.nan, 0)

matched_id = lambda x: process.extractOne(x, df['identifier'], scorer=fuzz.partial_ratio)[2]  # here the [2] is the index from the extractOne function, used for mapping, partial_ratio ensures only complete matches are used
E065_df['identifier'] = df.loc[E065_df["E065_id"].map(matched_id).values, 'identifier'].values  # takes the value from the prot_deg_id column, and matches against the identifier column df, then adds that column to df

df9 = pd.merge(df8, E065_df, on = 'identifier', how = 'left').fillna('0') # Merge the dfs, using the identifier column, fill IDs with no value in prot_deg with dashes

df9.to_csv('data7.csv', index=False)