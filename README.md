# Glyphosate Resistant


## Prepare Environment
```bash
module load conda/latest
conda create -n deepectransformer -f environment.yml
conda activate deepectransformer
pip install -e .

```
### Install DeepECTransformer
```bash
git clone https://github.com/jinghuazhao/DeepECTransformer.git
cd DeepECTransformer
pip install -e .
``` 

## Extract Top Species

### General top species (Without Metadata)

process_top_n_species(
    species_abundance=abundance_df,
    top_n=100,
    top_n_output_file="results/top_species"
)

### Top species split by Status


```python
python script/run_preprocess.py extract-top-species --help
```


Example:
python script/run_preprocess.py extract-top-species   
--input-file data/raw/Internal2_metaphlan_merged_abundance_table.txt   
--species-output-file data/processed/species_full_abundance_table.csv   
--metadata-file data/metadata/Internal_meta_diet.txt   
--top-n-output-file data/processed/top_specie

## Download from NCBI

```python
python script/run_preprocess.py download-top-taxa --help
```

### Convert FNA to FAA using Prodigal

```python
python script/run_prodigal.py --help
```

### Run deep EC transformer

```bash
cd DeepProZymed
sbatch sbatch_ec.sh --input ../data/downloads/general_taxa --output ../data/result/general
```

## Downstream analysis and visualization
check notebook postprocess.ipynb