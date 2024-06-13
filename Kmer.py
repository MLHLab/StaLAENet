###### K-mer Computation
###### Constructed K-mer matrices for each RGI cluster of each pathogen.
###### Need to run the program based on the total number of clusters.

!pip install biopython
import Bio
import pandas as pd
from Bio import SeqIO

def generate_sequences(k, bases=['A', 'C', 'G', 'T']):
    if k == 0:
        return ['']
    sequences = []
    for base in bases:
        for sequence in generate_sequences(k-1, bases):
            sequences.append(base + sequence)
    return sequences

def count_kmers(sequence, k):
    kmer_count = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in kmer_count:
            kmer_count[kmer] += 1
        else:
            kmer_count[kmer] = 1
    return kmer_count

def create_dataframe(k, fasta_file):
    sequences = generate_sequences(k)
    df = pd.DataFrame(columns=sequences)

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        kmer_counts = count_kmers(seq, k)
        df.loc[seq_id] = [kmer_counts.get(kmer, 0) for kmer in sequences]

    return df

k = 4
fasta_file = "/content/drive/MyDrive/NBU/Acinetobacter/cluster_57.fasta"
df = create_dataframe(k, fasta_file)
df.to_csv("/content/drive/MyDrive/NBU/Acinetobacter/cluster_57_4mer.csv")

###### 1. Assigned a column of 'ARO' to each kmer csv of each cluster.
###### 2. Assigned the mode of unique 'ARO' of a cluster to all the members of that cluster.
###### 3. Then, concat all the kmer clusters csv file.
###### 4. Also, concat a cluster of E-Coli with 'ARO' = 'Unclassified'.

import os
import pandas as pd

rgi_result = pd.read_csv('/content/drive/MyDrive/NBU/Acinetobacter/rgi_Acinetobacter.txt', sep = '\t')
csv_path = '/content/drive/MyDrive/NBU/Acinetobacter/4mer/'
ecoli_csv = pd.read_csv('/content/drive/MyDrive/NBU/Acinetobacter/cluster_57_4mer.csv')
output_path = '/content/drive/MyDrive/NBU/Acinetobacter/Acinetobacter_4mer.csv'
m = 257

df1 = rgi_result
folder_path = csv_path
df10 = ecoli_csv

df1.rename(columns = {'Contig':'Sequence_ID'}, inplace = True)
df1['Sequence_ID'] = df1['Sequence_ID'].str.strip()
df1['Sequence_ID'] = df1['Sequence_ID'].str[:-2]

csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]
merged_dfs_list = []

for file_name in csv_files:
    file_path = os.path.join(folder_path, file_name)

    df = pd.read_csv(file_path)
    df.rename(columns={'Unnamed: 0': 'Sequence_ID'}, inplace=True)

    merged_df = pd.merge(df, df1, left_on='Sequence_ID', right_on='Sequence_ID', how='left')
    merged_df = merged_df.iloc[:, :m].join(merged_df['ARO'])

    merged_df['ARO'].replace(merged_df['ARO'].unique(), merged_df['ARO'].mode()[0], inplace=True)
    merged_df['ARO'] = merged_df['ARO'].astype(int)
    unique_ARO_values = merged_df['ARO'].unique()
    merged_dfs_list.append(merged_df)
    print(f"Unique ARO values in {file_name}: {unique_ARO_values}")

concatenated_df = pd.concat(merged_dfs_list, ignore_index=True)


df10.rename(columns = {'Unnamed: 0':'Sequence_ID'}, inplace = True)
cl = []
for i in df10.index:
    cl.append('Unclassified')
df10['ARO'] = cl

final = pd.concat([concatenated_df, df10])
final.to_csv(output_path,index = 0)
