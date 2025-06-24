import requests
import pandas as pd
from io import StringIO
import typer
import re
import urllib.request
import os
from typing import Optional

app = typer.Typer()

def get_uniref90_for_accession(accession: str) -> str:
    """
    Fetch the UniRef90 cluster ID for a given UniProt accession.
    """
    url = "https://rest.uniprot.org/uniref/search"
    params = {
        "query": f'uniprotkb:{accession} AND identity:0.9',
        "format": "tsv",
        "fields": "id"
    }

    response = requests.get(url, params=params)
    if response.ok and response.text.strip():
        df = pd.read_csv(StringIO(response.text), sep="\t")
        if not df.empty and "Cluster ID" in df.columns:
            return df.iloc[0]["Cluster ID"]
        else:
            print(f"[Warning] 'Cluster ID' not found in response for accession {accession}")
            return ""

def search_uniprot_api(query: str, fields: str, uniref90: str = None) -> pd.DataFrame:
    """
    Helper function to query the UniProt API and process the response.

    Args:
        query (str): The query string for the UniProt API.
        fields (str): The fields to include in the response.
        uniref90 (str, optional): UniRef90 cluster ID to add as a column.

    Returns:
        pd.DataFrame: The resulting DataFrame from the API query.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "tsv",
        "fields": fields
    }

    response = requests.get(url, params=params)

    if response.ok and response.text.strip():
        # Read the response into a DataFrame
        df = pd.read_csv(StringIO(response.text), sep="\t")

        # Add the UniRef90 cluster ID as a new column if provided
        if uniref90:
            df["UniRef90"] = uniref90

        return df
    else:
        print(f"Error for query: {query}: {response.status_code}, {response.text}")
        return pd.DataFrame()

@app.command()
def search_uniprot_combined(
    organism: str,
    ec_number: Optional[str] = None,
    uniref90: Optional[str] = None
):
    """
    Search UniProt for entries matching the given organism, EC number, or UniRef90 cluster ID.

    Args:
        organism (str): The organism name to search for.
        ec_number (str, optional): The EC number to search for.
        uniref90 (str, optional): The UniRef90 cluster ID to search for.
    """
    # Construct the query based on the provided arguments
    query = f'organism_name:"{organism}"'
    if ec_number:
        query += f' AND ec:{ec_number}'
        print(f"Searching for EC number: {ec_number}")
    elif uniref90:
        query += f' AND uniref_cluster_90:{uniref90}'
        print(f"Searching for UniRef90 cluster ID: {uniref90}")

    # Define the fields to include in the response
    fields = "accession,id,protein_name,organism_name,ec"

    # Call the helper function
    df = search_uniprot_api(query=query, fields=fields, uniref90=uniref90)

    # Print and return the resulting DataFrame
    print(df)
    return df

def main(
    organism: str = typer.Option(...),
    ec_number: Optional[str] = typer.Option(None),
    uniref90: Optional[str] = typer.Option(None)
):
    search_uniprot_combined(organism, ec_number, uniref90)

if __name__ == "__main__":
    typer.run(main)

