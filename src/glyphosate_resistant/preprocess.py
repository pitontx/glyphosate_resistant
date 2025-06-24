import pandas as pd
import numpy as np
import typer
import re
from typing import Optional, Dict

typer_app = typer.Typer()

@typer_app.command()
def expand_taxonomic_column(df, clade_col= "clade_name", sep="|"):
    """
    Detects the column containing taxonomy strings, renames it to 'clade',
    and expands it into taxonomic rank columns.
    
    Parameters:
    - df: pandas DataFrame
    - sep: separator used in taxonomic strings (default: '|')
    
    Returns:
    - Expanded DataFrame with columns: clade, kingdom, phylum, class, order, family, genus, species, strain
    """
    rank_map = {
        "k__": "kingdom",
        "p__": "phylum",
        "c__": "class",
        "o__": "order",
        "f__": "family",
        "g__": "genus",
        "s__": "species",
        "t__": "strain"
    }


    # Step 3: Split and map taxonomic levels
    tax_split = df[clade_col].str.split(sep, expand=True)
    tax_df = pd.DataFrame()

    for col in tax_split.columns:
        for prefix, rank in rank_map.items():
            mask = tax_split[col].str.startswith(prefix, na=False)
            if mask.any():
                tax_df.loc[mask, rank] = tax_split.loc[mask, col].str.replace(prefix, '', regex=False)

    print(tax_df)
    # Step 4: Combine
    df_expanded = pd.concat([df, tax_df], axis=1)
    return df_expanded
    
def extract_species_abundance(input_file: str, output_file: str) -> pd.DataFrame:
    """
    Extracts species-level abundance data from a taxonomy file,
    expands full taxonomic levels, and writes to an output file.

    Parameters:
        input_file (str): Path to input taxonomic profile (e.g. from MetaPhlAn).
        output_file (str): Path to write the filtered + expanded data.

    Returns:
        pd.DataFrame: Expanded DataFrame including abundance and taxonomic columns.
    """
    try:
        # Read the input file into DataFrame (ignore comment lines)
        df = pd.read_csv(input_file, sep='\t', comment='#', header=0)

        # Expand taxonomic column
        df_expanded = expand_taxonomic_column(df.copy(), clade_col="clade_name")

        # Keep only species-level entries (must have s__ defined but no t__)
        is_species = df_expanded["species"].notna() & df_expanded["strain"].isna()
        df_filtered = df_expanded[is_species].copy()

        # Rename clade_name to "Species" and extract just species name (without prefix)
        df_filtered.rename(columns={"clade_name": "Species"}, inplace=True)
        df_filtered["Species"] = df_filtered["Species"].str.extract(r's__(\S+)', expand=False)

        # Save to file
        df_filtered.to_csv(output_file, sep='\t', index=False)
        print(f"Species-level abundance with taxonomy written to '{output_file}'.")

        return df_filtered

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return pd.DataFrame()
    except Exception as e:
        print(f"An error occurred: {e}")
        return pd.DataFrame()


from typing import Optional
import pandas as pd

def process_top_n_species(
    species_abundance: pd.DataFrame,
    top_n: int,
    top_n_output_file: str,
    metadata_file: Optional[str] = None
) -> Dict[str, pd.DataFrame]:
    """
    Process species abundance data to extract the top N abundant species.
    If metadata is provided, extract top N species per Status group.

    Parameters:
        species_abundance (pd.DataFrame): Species abundance table. First column must be 'Species', rest are sample columns.
        top_n (int): Number of top species to extract.
        top_n_output_file (str): Output file prefix. If metadata is used, appends _<Status>.tsv.
        metadata_file (str, optional): Metadata file with 'sample_name' and 'Status' columns.

    Returns:
        dict: Dictionary mapping status group (or "all") to the top N species DataFrame.
    """
    result_dict = {}

    if metadata_file:
        metadata = pd.read_csv(metadata_file, sep="\t")

        if "sample_name" not in metadata.columns or "Status" not in metadata.columns:
            raise ValueError("Metadata must contain 'sample_name' and 'Status' columns.")

        # Create a mapping: simplify abundance column names by extracting prefix before first "_"
        simplified_columns = []
        sample_mapping = {}

        for col in species_abundance.columns:
            if col == "Species":
                simplified_columns.append("Species")
            else:
                sample_prefix = col.split("_")[0]
                simplified_columns.append(sample_prefix)
                sample_mapping[sample_prefix] = col  # map simplified -> full name

        species_abundance.columns = simplified_columns

        for status, group in metadata.groupby("Status"):
            sample_ids = group["sample_name"].tolist()
            valid_sample_ids = [sid for sid in sample_ids if sid in species_abundance.columns]

            if not valid_sample_ids:
                print(f"Warning: No valid samples found for status '{status}'. Skipping.")
                continue

            df_subset = species_abundance[["Species"] + valid_sample_ids].copy()
            # Compute average instead of total
            df_subset["Average_Abundance"] = df_subset.iloc[:, 1:].mean(axis=1)
        
            top_species = df_subset[["Species", "Average_Abundance"]] \
                .sort_values(by="Average_Abundance", ascending=False).head(top_n)

            out_file = f"{top_n_output_file}_{status}.tsv"
            top_species.to_csv(out_file, sep="\t", index=False)

            print(f"Top {top_n} species for status '{status}' saved to {out_file}")
            result_dict[status] = top_species

    else:
        df = species_abundance.copy()
        df["Average_Abundance"] = df.iloc[:, 1:].mean(axis=1)

        top_species = df[["Species", "Average_Abundance"]] \
            .sort_values(by="Average_Abundance", ascending=False).head(top_n)

        out_file = f"{top_n_output_file}_all.tsv"
        top_species.to_csv(out_file, sep="\t", index=False)

        print(f"Top {top_n} species across all samples saved to {out_file}")
        result_dict["all"] = top_species

    return result_dict