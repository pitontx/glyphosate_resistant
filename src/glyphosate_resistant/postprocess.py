import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import matplotlib.pyplot as plt
from glyphosate_resistant.preprocess import extract_species_abundance,expand_taxonomic_column

def summarize_ec_completeness(combined_df: pd.DataFrame, enzymes: dict, status: str) -> pd.DataFrame:
    """
    Calculate and summarize EC number detection completeness for each species.

    Parameters:
        combined_df (pd.DataFrame): DataFrame containing at least 'Species' and 'EC number' columns.
        enzymes (dict): Dictionary where values are EC numbers of interest.

    Returns:
        pd.DataFrame: Summary DataFrame with species, detected ECs, Status, and completeness percentage.
    """
    target_ecs = set(enzymes.values())

    # Filter for relevant EC numbers
    df_filtered = combined_df[combined_df['EC number'].isin(target_ecs)]

    # Group by species and collect EC sets
    species_ec_sets = df_filtered.groupby('Species')['EC number'].apply(set)

    # Calculate completeness per species
    completeness = species_ec_sets.apply(lambda ecs: len(ecs) / len(target_ecs) * 100)
    Status = status

    # Assemble result DataFrame
    summary = pd.DataFrame({
        'Detected ECs': species_ec_sets,
        'Completeness (%)': completeness,
        'Status': Status
    }).reset_index()

    # Sort by completeness
    summary = summary.sort_values(by='Completeness (%)', ascending=False)

    return summary

def load_deepec_results(path_pattern: str) -> pd.DataFrame:
    """
    Load and combine DeepECv2 result files from multiple species directories.

    Parameters:
        path_pattern (str): Glob pattern to match DeepECv2_result.txt files.

    Returns:
        pd.DataFrame: Combined DataFrame with EC number and species annotations.
    """
    file_paths = glob.glob(path_pattern)
    all_results = []

    for file_path in file_paths:
        species_name = os.path.basename(os.path.dirname(file_path))
        df = pd.read_csv(file_path, sep="\t")

        # Extract EC number
        df['EC number'] = df['prediction'].str.split(":").str[1]

        # Add species column
        df['Species'] = species_name

        all_results.append(df)

    if not all_results:
        raise ValueError(f"No DeepEC result files found with pattern: {path_pattern}")

    combined_result = pd.concat(all_results, ignore_index=True)
    return combined_result


def summarize_taxa_completion(df, taxa_col="genus", completeness_col="Completeness (%)", threshold=100.0):
    """
    Summarize how many species per genus have complete vs incomplete pathway.

    Returns a DataFrame like:
    | genus        | Complete | Incomplete |
    |--------------|----------|------------|
    | Bacteroides  |     3    |     1      |
    | Roseburia    |     2    |     0      |
    """
    # Clean
    df = df.dropna(subset=[taxa_col])
    df["Status"] = df[completeness_col].apply(lambda x: "Complete" if x >= threshold else "Incomplete")

    # Count species per genus per status
    summary = (
        df.groupby([taxa_col, "Status"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    # Ensure both columns exist
    for col in ["Complete", "Incomplete"]:
        if col not in summary.columns:
            summary[col] = 0

    return summary


def plot_taxa_completion_horizontal(summary_df, taxa_col="genus"):
    """
    Plot horizontal grouped bar chart: Complete vs Incomplete counts per genus.
    """
    # Sort by total number of species per genus
    summary_df["Total"] = summary_df["Complete"] + summary_df["Incomplete"]
    summary_df = summary_df.sort_values("Total", ascending=True)

    # Plot
    ax = summary_df.set_index(taxa_col)[["Complete", "Incomplete"]].plot(
        kind="barh", figsize=(8, 8), color=["#4caf50", "#f44336"]
    )

    plt.xlabel("Number of Species")
    plt.ylabel("Genus")
    plt.title("Pathway Completion per taxa")
    plt.legend(title="Pathway Status", loc="lower right")
    plt.tight_layout()
    plt.show()

def stacked_bar_multiple_group(df, level="phylum", abundance_col="Average_Abundance", group_col="Status", top_n=10):
    """
    Create a stacked bar plot comparing taxonomic composition across multiple groups.
    
    Parameters:
        df: DataFrame with taxonomic levels, abundance, and group info
        level: Taxonomic level to summarize ('phylum', 'genus', etc.)
        abundance_col: Column with abundance values
        group_col: Column that defines group (e.g. 'Status')
        top_n: Number of top taxa (by total abundance across all groups) to display individually
    """
    # Sum by taxon across all groups to find top N taxa
    top_taxa = (
        df.groupby(level)[abundance_col].sum()
        .sort_values(ascending=False)
        .head(top_n)
        .index
        .tolist()
    )

    # Label non-top taxa as "Other"
    df[level] = df[level].apply(lambda x: x if x in top_taxa else "Other")

    # Group again: now collapsed minor taxa
    grouped = (
        df.groupby([group_col, level])[abundance_col].sum()
        .reset_index()
    )

    # Pivot to wide format: rows = groups, columns = taxa
    pivot_df = (
        grouped.pivot(index=group_col, columns=level, values=abundance_col)
        .fillna(0)
    )

    # Sort columns by total abundance for legend consistency
    pivot_df = pivot_df[pivot_df.sum().sort_values(ascending=False).index]

    # Plot
    pivot_df.plot(kind="bar", stacked=True, figsize=(10, 5), colormap="tab20")
    plt.ylabel("Average Abundance")
    plt.title(f"{level.title()} Composition of top species by Group")
    plt.xticks(rotation=0)
    plt.legend(title=level.title(), bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()

