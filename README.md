# Melanoma-Metastasis-Gene-Expression-Vs.-Age-Analysis Code

# (1.) Filtering relevant data from datasets and merging them into merged_data DataFrame  


import pandas as pd

metadata = pd.read_csv('GSE62944_metadata.csv') # clinical metadata
rnadata = pd.read_csv('GSE62944_subsample_log2TPM.csv') # RNA-seq log2(TPM) data

# Normalize rnadata so genes are the index (file has genes in 'Unnamed: 0')
if 'Unnamed: 0' in rnadata.columns:
    rnadata = rnadata.rename(columns={'Unnamed: 0': 'gene'}).set_index('gene')
elif 'gene' in rnadata.columns and rnadata.index.name != 'gene':
    rnadata = rnadata.set_index('gene')

# Filter for just melanoma (SKCM) samples
skcm_metadata = metadata[metadata['cancer_type'] == 'SKCM'].copy()

# Prefer age_at_diagnosis if age_at_initial_pathologic_diagnosis is missing
if ('age_at_initial_pathologic_diagnosis' not in skcm_metadata.columns
        or skcm_metadata['age_at_initial_pathologic_diagnosis'].isna().all()):
    if 'age_at_diagnosis' in skcm_metadata.columns:
        skcm_metadata['age_at_initial_pathologic_diagnosis'] = skcm_metadata['age_at_diagnosis']

# Standardize sample id and sample type column names used downstream
if 'sample' in skcm_metadata.columns:
    skcm_metadata = skcm_metadata.rename(columns={'sample': 'sample_id'})
if 'submitted_tumor_site' in skcm_metadata.columns:
    skcm_metadata = skcm_metadata.rename(columns={'submitted_tumor_site': 'sample_type'})

# Filter for the clinical features we need (now using standardized names)
skcm_metadata = skcm_metadata[['sample_id', 'age_at_initial_pathologic_diagnosis', 'sample_type']]

# Filter for relevant genes and check if they are in the rnadata index
genes = [
    'BRAF','NRAS','AXL','MITF','MMP2','MMP9','MMP14','FN1','VIM','TGFB1',
    'CASP8','KISS1','PTEN','TERT','CDKN2A','NEDD9'
]
present_genes = [gene for gene in genes if gene in rnadata.index]
print(f"Genes present in RNA-seq data: {present_genes}")
missing_genes = [gene for gene in genes if gene not in rnadata.index]
print(f"Genes missing in RNA-seq data: {missing_genes}")

# Adjust columns and rows of rnadata to merge the metadata with the RNA-seq data
common_samples = sorted(set(skcm_metadata['sample_id']) & set(rnadata.columns))
# Use only present_genes (avoids KeyError if some genes missing)
rnadata_subset = rnadata.loc[present_genes, common_samples].T.reset_index().rename(columns={'index': 'sample_id'})
merged_data = pd.merge(skcm_metadata, rnadata_subset, on='sample_id', how='inner')

print() # For spacing
print(merged_data.head())
