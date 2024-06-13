###### Core Cluster Finding
###### Such clusters which have sequences of all the strains of a specific pathogen.
###### Got cores_with_orthologs.fasta and cores_with_orthologs.tsv files for each cluster
###### Here we have deleted the duplicates and retain the orthologs of each cluster.
###### Each cluster has equal number of sequences which came from all the strains of a specific pathogen.
###### cores_with_orthologs.tsv file includes Cluster_ids, Sequence_ids, Representatives, Identity, Strain_ids etc.
###### cores_with_orthologs.fasta file includes all the sequences which are present in the core clusters of that pathogen.

import os
import sys
import pyfaidx
import pandas as pd
from tqdm import tqdm
from operator import itemgetter

def get_cores(clstr, count):

	# Extracting strain count for given species
	# count = pd.read_csv(count, sep='\t')
	# count = count.shape[0]

	# Change this to shape[0] for species
	count = int(count)

	clstr = pd.read_csv(clstr, sep='\t')
	all_tbls = []

	# Removing cluster_IDs where there are less than `count` members
	clstr = clstr.groupby('Cluster_ID').filter(lambda x: len(x) >= count)

	# Check which of the following clusters contain at least `count` members
	for clust in tqdm(clstr.Cluster_ID.unique()):
		tbl_clust = clstr[clstr.Cluster_ID==clust]

		if len(tbl_clust.Strain_ID.unique()) >= count:
			idx = tbl_clust.groupby('Strain_ID')['Identity'].idxmax()
			result = tbl_clust.loc[idx]
			all_tbls.append(result)

	df = pd.concat(all_tbls)
	# Extract core fasta
	return df

def get_core_fasta(cores, fasta, out):

	core_sequences=list(cores.Sequence_ID.unique())
	core_sequences=[c[1:] for c in core_sequences]
	print(f"Core sequences: {len(core_sequences)}")
	fasta_count=0

	data={}
	records = pyfaidx.Fasta(fasta)
	for record in tqdm(records, ncols=100):
		data[record.name] = record

	# Use itemgetter to extract values based on keys
	value_extractor = itemgetter(*core_sequences)
	subset_values = value_extractor(data)

	# Create a new dictionary using zip
	subset_dict = dict(zip(core_sequences, subset_values))
	with open(out, 'w') as cf:
		for k, v in tqdm(subset_dict.items(), ncols=100):
			fasta_count+=1
			cf.write(">"+v.long_name+"\n")
			cf.write(str(v)+"\n")

	print(f"Fasta count: {fasta_count}")

if __name__=="__main__":

	# TSV of strain_ID tagged clusters
	clstr_df= 'D:VRITIKA_New/DATA/Acinetobacter/clstr_df_strain.tsv'
	# Fasta file extracted from CD-HIT
	clstr_fasta = 'D:VRITIKA_New/DATA/Acinetobacter/Acinetobacter_all_strains.fasta.fai'
	# Metadata of strain count
	strain_count= 228

	cores = get_cores(clstr_df, strain_count)
	cores.to_csv(os.path.join(os.path.dirname(clstr_df), "cores_with_orthologs.tsv"), sep='\t', index=False)
	print(f'Saved to {os.path.join(os.path.dirname(clstr_df), "cores_with_orthologs.tsv")}')

	out = os.path.join(os.path.dirname(clstr_df), "cores_with_orthologs.fasta")
	get_core_fasta(cores, clstr_fasta, out)

