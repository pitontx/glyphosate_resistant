import os
import gzip
import shutil
import subprocess
import typer
from loguru import logger

app = typer.Typer()

def decompress_gz(input_path: str, output_path: str):
    with gzip.open(input_path, 'rt') as f_in, open(output_path, 'w') as f_out:
        shutil.copyfileobj(f_in, f_out)

def compress_gz(input_path: str, output_path: str):
    with open(input_path, 'rb') as f_in, gzip.open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def run_prodigal(input_fna: str, output_faa: str):
    prodigal_cmd = [
        "prodigal",
        "-i", input_fna,
        "-a", output_faa,
        "-p", "meta"
    ]
    subprocess.run(prodigal_cmd, check=True)
    
def process_species_folder(species_dir: str):
    faa_exists = False
    fna_path = None
    fna_gz_path = None

    for file in os.listdir(species_dir):
        if file.endswith(".faa") or file.endswith(".faa.gz"):
            faa_exists = True
        elif file.endswith(".fna.gz"):
            fna_gz_path = os.path.join(species_dir, file)
        elif file.endswith(".fna"):
            fna_path = os.path.join(species_dir, file)

    if faa_exists:
        logger.info(f".faa or .faa.gz already exists in {species_dir}, skipping.")
        return

    # If both .fna and .fna.gz are present, prioritize .fna
    if not fna_path and fna_gz_path:
        logger.info(f"Decompressing {fna_gz_path} ...")
        base_name = os.path.splitext(os.path.splitext(os.path.basename(fna_gz_path))[0])[0]
        fna_path = os.path.join(species_dir, base_name + ".fna")
        decompress_gz(fna_gz_path, fna_path)

    if not fna_path:
        logger.warning(f"No .fna or .fna.gz file found in {species_dir}, skipping.")
        return

    logger.info(f"Running Prodigal on {fna_path} ...")

    base_name = os.path.splitext(os.path.basename(fna_path))[0]
    faa_path = os.path.join(species_dir, base_name + ".faa")
    faa_gz_path = faa_path + ".gz"

    run_prodigal(fna_path, faa_path)
    compress_gz(faa_path, faa_gz_path)

    logger.info(f"Created: {faa_gz_path}")



@app.command()
def main(taxa_dir: str = typer.Argument(..., help="Path to the directory containing species subfolders")):
    for species_name in os.listdir(taxa_dir):
        species_dir = os.path.join(taxa_dir, species_name)
        if os.path.isdir(species_dir):
            process_species_folder(species_dir)

if __name__ == "__main__":
    app()
