import os
import re
import urllib.request
from Bio import Entrez, SeqIO
import typer
from typer import Typer
from typing import Optional
from loguru import logger
import ssl, certifi
import time  # <-- Add this line

app = Typer()

# Set your email for NCBI Entrez          
Entrez.email = "hangyin1993@gmail.com"
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

def resolve_species_name(query):
    """
    Try to find the closest matching species from NCBI Taxonomy
    that also has genome assemblies available.
    """
    query = query.replace("_", " ")
    logger.info(f"Resolving species name for query: {query}")

    # Search taxonomy for the query
    handle = Entrez.esearch(db="taxonomy", term=query, retmode="xml", retmax=20)
    records = Entrez.read(handle)
    time.sleep(0.5)
    handle.close()

    if not records["IdList"]:
        logger.error(f"No taxonomy match found for query: {query}")
        return None, None

    for taxid in records["IdList"]:
        # Check if this taxid has at least one assembly
        handle = Entrez.esearch(
            db="assembly",
            term=f"txid{taxid}[Organism:exp] AND latest[filter] AND all[filter]",
            retmode="xml",
            retmax=1
        )
        assembly_check = Entrez.read(handle)
        time.sleep(0.5)
        handle.close()

        if assembly_check["IdList"]:
            # Fetch the resolved scientific name
            summary_handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            summary = Entrez.read(summary_handle)
            time.sleep(0.5)
            summary_handle.close()

            species_name = summary[0]["ScientificName"]
            logger.info(f"Resolved species with available assembly: {species_name} (TaxID: {taxid})")
            return species_name, taxid

    logger.error(f"No matching species with assemblies found for query: {query}")
    return None, None



def download_genome(species_query: str, output_dir: str):
    output_species_dir = os.path.join(output_dir, species_query.replace(" ", "_"))
    if os.path.exists(output_species_dir):
        print(f"Species '{species_query}' already downloaded in {output_species_dir}, skipping.")
        return

    species_name, taxid = resolve_species_name(species_query)
    if not species_name:
        return

    print(f"Resolved species name: {species_name} (TaxID: {taxid})")

    # Step 1: Search for assemblies by taxid
    term = f"txid{taxid}[Organism:exp] AND latest[filter] AND all[filter]"
    handle = Entrez.esearch(db="assembly", term=term, retmax=1)
    record = Entrez.read(handle)
    time.sleep(0.5)
    handle.close()

    if not record["IdList"]:
        print(f"No assemblies found for taxid {taxid}")
        return

    assembly_uid = record["IdList"][0]
    print(f"Assembly UID: {assembly_uid}")

    # Step 2: Fetch assembly summary to get FTP link
    summary_handle = Entrez.esummary(db="assembly", id=assembly_uid, report="full", retmode="xml")
    summary_record = Entrez.read(summary_handle)
    time.sleep(0.5)
    summary_handle.close()

    docsum = summary_record["DocumentSummarySet"]["DocumentSummary"][0]
    assembly_accession = docsum["AssemblyAccession"]
    assembly_name = docsum["AssemblyName"]
    ftp_path = docsum.get("FtpPath_RefSeq") or docsum.get("FtpPath_GenBank")

    if not ftp_path:
        print("No FTP path found for the selected assembly.")
        return

    print(f"Using assembly: {assembly_accession} ({assembly_name})")
    print(f"FTP path: {ftp_path}")

    # Step 3: Scan FTP directory for desired files
    base_url = ftp_path + "/"
    try:
        with urllib.request.urlopen(base_url) as response:
            html = response.read().decode()
    except Exception as e:
        print(f"Failed to access FTP directory: {e}")
        return

    pattern = re.compile(r'(GC[AF]_[^/]+_(genomic|protein|cds_from_genomic)[^/]*\.gz)')
    files_to_download = []
    for line in html.splitlines():
        match = pattern.search(line)
        if match:
            files_to_download.append(match.group(0))

    if not files_to_download:
        print("No genomic or protein files found in FTP directory.")
        return

    print(f"Found {len(files_to_download)} file(s) to download.")

    os.makedirs(output_species_dir, exist_ok=True)

    for filename in files_to_download:
        url = base_url + filename
        local_path = os.path.join(output_species_dir, filename)
        print(f"Downloading {filename} from {url} ...")
        try:
            urllib.request.urlretrieve(url, local_path)
            print(f"Saved to {local_path}")
        except Exception as e:
            print(f"Failed to download {filename}: {e}")

    print("Download completed.")


@app.command()
def cli_download_genome(species_name: str, output_dir: str = "./data/downloads"):
    """
    Download the genome assembly for a given species and save it to a dynamically named file.

    Args:
        species_name (str): The species name or query to search for.
        output_dir (str): Directory to save the downloaded files.
    """
    download_genome(species_name, output_dir)

if __name__ == "__main__":
    app()