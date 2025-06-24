from glyphosate_resistant.preprocess import extract_species_abundance, process_top_n_species
from glyphosate_resistant.uniprot_search import search_uniprot_combined
from glyphosate_resistant.download_ncbi import download_genome
import pandas as pd
import typer
from io import StringIO
import concurrent.futures
from loguru import logger
app = typer.Typer()
from typing import Optional


enzymes = {
    "aroF_G_H": "2.5.1.54",
    "aroB": "4.2.3.4",
    "aroD": "4.2.1.10",
    "aroE": "1.1.1.25",
    "aroK_L": "2.7.1.71",
    "aroA": "2.5.1.19",
    "aroC": "4.2.3.5"
}

@app.command()
def extract_top_species(
    input_file: str = typer.Option(..., help="Path to the input file containing species abundance data."),
    species_output_file: Optional[str] = typer.Option(None, help="Path to save filtered species abundance data."),
    metadata_file: Optional[str] = typer.Option(None, help="Path to metadata file with 'SampleID' and 'Status'."),
    top_n_output_file: Optional[str] = typer.Option("top_species", help="Prefix for output file(s) with top abundant species.")
):
    """
    Extract species abundance and compute the top 100 abundant species.

    If a metadata file is provided, the top 100 will be computed per status group.
    Otherwise, the top 100 species will be computed from all samples combined.

    Parameters:
    - input_file: Tab-separated file containing raw species abundance data.
    - species_output_file: Where to save the cleaned species abundance data.
    - metadata_file: Optional file with 'SampleID' and 'Status' to split samples.
    - top_n_output_file: Output file prefix for top 100 species result(s).
    """
    try:
        logger.info(f"Reading species abundance data from {input_file}")
        species_abundance = extract_species_abundance(input_file, species_output_file)

        if species_abundance.empty:
            logger.warning("No data to process after filtering.")
            raise typer.Exit(code=1)

        logger.info("Computing top 100 abundant species...")
        result = process_top_n_species(
            species_abundance=species_abundance,
            top_n=100,
            top_n_output_file=top_n_output_file,
            metadata_file=metadata_file
        )

        logger.info(f"Top 100 species results saved with prefix '{top_n_output_file}'")
        return result

    except Exception as e:
        logger.error(f"Failed to extract top species: {e}")
        raise typer.Exit(code=1)


@app.command()
def combine_outputs(
    top_n_species_file: str = typer.Option(...,help="Path to the file containing top N species."),
    output_file: str = typer.Option(None, help="Path to save the concatenated UniProt search results."),
    abundance_file: str = typer.Option(..., help="Path to the file containing species abundance data.")
):
    """
    Loop through species and enzymes, search UniProt, and concatenate results.

    Args:
        top_n_species_file (str): Path to the file containing top N species.
        output_file (str): Path to save the concatenated UniProt search results.
    """
    logger.info(f"Loading top N species from {top_n_species_file}")
    top_n_species = pd.read_csv(top_n_species_file, sep="\t")

    all_results = pd.DataFrame()

    for _, row in top_n_species.iterrows():
        species = row["Species"]
        for enzyme, ec_number in enzymes.items():
            logger.info(f"Searching UniProt for species: {species}, enzyme: {enzyme}, EC: {ec_number}")
            try:
                result_df = search_uniprot_combined(organism=species, ec_number=ec_number)
                if not result_df.empty:
                    result_df["Species"] = species
                    result_df["Enzyme"] = enzyme
                    all_results = pd.concat([all_results, result_df], ignore_index=True)
            except Exception as e:
                logger.error(f"Error searching UniProt for {species}, {enzyme}, {ec_number}: {e}")

    if output_file:
        all_results.to_csv(output_file, sep="\t", index=False)
        logger.info(f"Concatenated UniProt search results have been saved to '{output_file}'.")
    return all_results

def load_or_generate_results(
    result_file: str,
    top_n_species_file: str,
    abundance_file: str,
    output_file: str
) -> pd.DataFrame:
    """
    Load an existing result file or generate it using combine_outputs.

    Args:
        result_file (str): Path to the input file containing UniProt search results.
        top_n_species_file (str): Path to the file containing top N species.
        abundance_file (str): Path to the file containing species abundance data.
        output_file (str): Path to save the generated results.

    Returns:
        pd.DataFrame: The resulting DataFrame.
    """
    if result_file:
        logger.info(f"Loading existing result file: {result_file}")
        return pd.read_csv(result_file, sep="\t")
    else:
        logger.info("No result file provided. Generating result file...")
        return combine_outputs(
            top_n_species_file=top_n_species_file,
            output_file=output_file,
            abundance_file=abundance_file
        )

@app.command()
def add_ec_from_uniref(
    result_file: str = typer.Option(None, help="Path to the input file containing UniProt search results with UniRef90 IDs if it exists."),
    top_n_species_file: str = typer.Option(..., help="Path to the file containing top N species, required if there's no UniProt search results."),
    abundance_file: str = typer.Option(None, help="Path to the file containing species abundance data, required if there's no UniProt search results."),
    output_file: str = typer.Option(..., help="Path to save the results with UniRef90 EC numbers."),
    max_workers: int = typer.Option(8, help="Number of threads to use for UniProt queries.")
):
    """
    Add EC numbers from UniRef90 to the result file, using multithreading for faster API queries.
    """
    initial_df = load_or_generate_results(
        result_file=result_file,
        top_n_species_file=top_n_species_file,
        abundance_file=abundance_file,
        output_file=output_file
    )

    if "UniRef90" in initial_df.columns:
        species_list = initial_df["Species"].unique().tolist()
        uniRef90_list = initial_df["UniRef90"].unique().tolist()
        required_ecs = set(enzymes.values())

        filtered_species = [
            species for species in species_list
            if not required_ecs.issubset(
                set(initial_df.loc[initial_df["Species"] == species, "EC number"].dropna().unique())
            )
        ]

        all_results = []

        def worker(args):
            species, uniref90 = args
            try:
                logger.info(f"Searching UniProt for species: {species}, uniref90: {uniref90}")
                df = search_uniprot_combined(organism=species, uniref90=uniref90)
                if not df.empty:
                    df["Species"] = species
                return df
            except Exception as e:
                logger.error(f"Error searching UniProt for {species}, {uniref90}: {e}")
                return pd.DataFrame()

        tasks = [(species, uniref90) for species in filtered_species for uniref90 in uniRef90_list]

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            for result_df in executor.map(worker, tasks):
                if result_df is not None and not result_df.empty:
                    all_results.append(result_df)

        if all_results:
            final_df = pd.concat(all_results, ignore_index=True)
        else:
            final_df = pd.DataFrame()

        final_df.to_csv(output_file, sep="\t", index=False)
        logger.info(f"Concatenated UniProt search results have been saved to '{output_file}'.")

@app.command()
def download_top_taxa(
    top_n_species_file: str = typer.Option(None, help="Path to the file containing top N species."),
    abundance_file: str = typer.Option(None, help="Path to the file containing species abundance data if top N species file is empty."),
    output_dir: str = typer.Option("./data/downloads", help="Directory to save downloaded files.")
):
    """
    Download genome assemblies for the top N species.
    """
    if not top_n_species_file:
        logger.info(f"Extracting top species from abundance file: {abundance_file}")
        top_n_species = extract_top_species(input_file = abundance_file)
    else:
        logger.info(f"Loading top N species from {top_n_species_file}")
        top_n_species = pd.read_csv(top_n_species_file, sep="\t")
    
    for _, row in top_n_species.iterrows():
        species_name = row["Species"]
        logger.info(f"Downloading genome assembly for {species_name}...")
        try:
            download_genome(
                species_query=species_name,
                output_dir=output_dir
            )
            logger.info(f"Downloaded genome assembly for {species_name} to {output_dir}")
        except Exception as e:
            logger.error(f"Failed to download genome assembly for {species_name}: {e}")

if __name__ == "__main__":
    app()
