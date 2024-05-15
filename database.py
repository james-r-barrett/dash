import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis

filename = "230315_SAG211_8k_annotated_proteome_combined (1).fasta"

# code chunk adapted from David Cain (stackoverflow)
with open(filename) as fasta_file:  # Will close handle cleanly
    identifiers = []
    descriptions = []
    lengths = []
    sequences = []
    isoelectrics = []
    weights = []
    for title, sequence in SimpleFastaParser(fasta_file):
        identifiers.append(title.split(None, 1)[0])  # First word is ID
        descriptions.append(title.split(None, 1)[1])  # First word is ID
        sequences.append(sequence.replace("*",""))
        x = ProteinAnalysis(sequence.replace("*",""))
        isoelectrics.append(x.isoelectric_point())
        weights.append(x.molecular_weight()/1000)
        lengths.append(len(sequence.replace("*","")))
    print(str(len(identifiers))+str(" sequences parsed successfully."))

df = pd.DataFrame(zip(identifiers, descriptions, lengths, isoelectrics, weights, sequences), columns=['identifier', 'description', 'length', 'isoelectric','weight (kDa)','sequence'])

df.to_csv('dataframe.csv', index=False)