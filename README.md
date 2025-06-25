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
!python script/run_preprocess.py extract-top-species --help
```

    [1m                                                                                [0m
    [1m [0m[1;33mUsage: [0m[1mrun_preprocess.py extract-top-species [OPTIONS][0m[1m                        [0m[1m [0m
    [1m                                                                                [0m
     Extract species abundance and compute the top 100 abundant species.            
                                                                                    
     [2mIf a metadata file is provided, the top 100 will be computed per status group.[0m 
     [2mOtherwise, the top 100 species will be computed from all samples combined.[0m     
     [2mParameters: - input_file: Tab-separated file containing raw species abundance [0m 
     [2mdata. - species_output_file: Where to save the cleaned species abundance data.[0m 
     [2m- metadata_file: Optional file with 'SampleID' and 'Status' to split samples. [0m 
     [2m- top_n_output_file: Output file prefix for top 100 species result(s).[0m         
                                                                                    
    [2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────[0m[2m─╮[0m
    [2m│[0m [31m*[0m  [1;36m-[0m[1;36m-input[0m[1;36m-file[0m                 [1;33mTEXT[0m  Path to the input file containing      [2m│[0m
    [2m│[0m                                       species abundance data.                [2m│[0m
    [2m│[0m                                       [2m[default: None]                       [0m [2m│[0m
    [2m│[0m                                       [2;31m[required]                            [0m [2m│[0m
    [2m│[0m    [1;36m-[0m[1;36m-species[0m[1;36m-output-file[0m        [1;33mTEXT[0m  Path to save filtered species          [2m│[0m
    [2m│[0m                                       abundance data.                        [2m│[0m
    [2m│[0m                                       [2m[default: None]                       [0m [2m│[0m
    [2m│[0m    [1;36m-[0m[1;36m-metadata[0m[1;36m-file[0m              [1;33mTEXT[0m  Path to metadata file with 'SampleID'  [2m│[0m
    [2m│[0m                                       and 'Status'.                          [2m│[0m
    [2m│[0m                                       [2m[default: None]                       [0m [2m│[0m
    [2m│[0m    [1;36m-[0m[1;36m-top[0m[1;36m-n-output-file[0m          [1;33mTEXT[0m  Prefix for output file(s) with top     [2m│[0m
    [2m│[0m                                       abundant species.                      [2m│[0m
    [2m│[0m                                       [2m[default: top_species]                [0m [2m│[0m
    [2m│[0m    [1;36m-[0m[1;36m-help[0m                       [1;33m    [0m  Show this message and exit.            [2m│[0m
    [2m╰──────────────────────────────────────────────────────────────────────────────╯[0m
    


Example:
python script/run_preprocess.py extract-top-species   
--input-file data/raw/Internal2_metaphlan_merged_abundance_table.txt   
--species-output-file data/processed/species_full_abundance_table.csv   
--metadata-file data/metadata/Internal_meta_diet.txt   
--top-n-output-file data/processed/top_specie

## Download from NCBI


```python
!python script/run_preprocess.py download-top-taxa --help
```

    [1m                                                                                [0m
    [1m [0m[1;33mUsage: [0m[1mrun_preprocess.py download-top-taxa [OPTIONS][0m[1m                          [0m[1m [0m
    [1m                                                                                [0m
     Download genome assemblies for the top N species.                              
                                                                                    
                                                                                    
    [2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────[0m[2m─╮[0m
    [2m│[0m [1;36m-[0m[1;36m-top[0m[1;36m-n-species-file[0m        [1;33mTEXT[0m  Path to the file containing top N species. [2m│[0m
    [2m│[0m                                   [2m[default: None]                           [0m [2m│[0m
    [2m│[0m [1;36m-[0m[1;36m-abundance[0m[1;36m-file[0m            [1;33mTEXT[0m  Path to the file containing species        [2m│[0m
    [2m│[0m                                   abundance data if top N species file is    [2m│[0m
    [2m│[0m                                   empty.                                     [2m│[0m
    [2m│[0m                                   [2m[default: None]                           [0m [2m│[0m
    [2m│[0m [1;36m-[0m[1;36m-output[0m[1;36m-dir[0m                [1;33mTEXT[0m  Directory to save downloaded files.        [2m│[0m
    [2m│[0m                                   [2m[default: ./data/downloads]        [0m        [2m│[0m
    [2m│[0m [1;36m-[0m[1;36m-help[0m                      [1;33m    [0m  Show this message and exit.                [2m│[0m
    [2m╰──────────────────────────────────────────────────────────────────────────────╯[0m
    


### Convert FNA to FAA using Prodigal


```python
!python script/run_prodigal.py --help
```

    [1m                                                                                [0m
    [1m [0m[1;33mUsage: [0m[1mrun_prodigal.py [OPTIONS] TAXA_DIR[0m[1m                                     [0m[1m [0m
    [1m                                                                                [0m
    [2m╭─[0m[2m Arguments [0m[2m─────────────────────────────────────────────────────────────────[0m[2m─╮[0m
    [2m│[0m [31m*[0m    taxa_dir      [1;33mTEXT[0m  Path to the directory containing species subfolders [2m│[0m
    [2m│[0m                          [2m[default: None]                                    [0m [2m│[0m
    [2m│[0m                          [2;31m[required]                                         [0m [2m│[0m
    [2m╰──────────────────────────────────────────────────────────────────────────────╯[0m
    [2m╭─[0m[2m Options [0m[2m───────────────────────────────────────────────────────────────────[0m[2m─╮[0m
    [2m│[0m [1;36m-[0m[1;36m-install[0m[1;36m-completion[0m          Install completion for the current shell.      [2m│[0m
    [2m│[0m [1;36m-[0m[1;36m-show[0m[1;36m-completion[0m             Show completion for the current shell, to copy [2m│[0m
    [2m│[0m                               it or customize the installation.              [2m│[0m
    [2m│[0m [1;36m-[0m[1;36m-help[0m                        Show this message and exit.                    [2m│[0m
    [2m╰──────────────────────────────────────────────────────────────────────────────╯[0m
    


### Run deep EC transformer

```bash
cd DeepProZymed
sbatch sbatch_ec.sh --input ../data/downloads/general_taxa --output ../data/result/general
```

## Downstream analysis and visualization
check notebook postprocess.ipynb