"""Ideogram plotting executables."""
from enum import Enum
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

from ideogram.plotting_utils import plot_ideogram


class GENOME(Enum):
    HG38 = "hg38"
    CHM13 = "chm13"
    HS1 = "hs1"


COLOUR_LOOKUP = {
    "gneg": (1.0, 1.0, 1.0),
    "gpos25": (0.6, 0.6, 0.6),
    "gpos50": (0.4, 0.4, 0.4),
    "gpos75": (0.2, 0.2, 0.2),
    "gpos100": (0.0, 0.0, 0.0),
    # 'acen': (.8, .4, .4),
    # Set acen to be white as we use a
    #   polygon to add it in later
    "acen": (1.0, 1.0, 1.0),
    "gvar": (0.8, 0.8, 0.8),
    "stalk": (0.9, 0.9, 0.9),
}
STATIC_PATH = Path(__file__).parent / "static"


def get_cytobands(genome: GENOME) -> Path:
    """
    Return the cytobands file for the given genome.

    :param genome: The genome variant to get the cytobands file for.
    :return: The path to the cytobands file associated with the provided genome variant.
    :raises ValueError: If the provided genome variant is not recognized.

    >>> get_cytobands(GENOME.HG38)
    PosixPath('/.../ideogram/src/ideogram/static/cytobands_HG38.bed')
    >>> get_cytobands(GENOME.CHM13)
    PosixPath('/.../ideogram/src/ideogram/static/cytobands_chm13.bed')
    >>> get_cytobands(GENOME.HS1)
    PosixPath('/.../ideogram/src/ideogram/static/cytobands_chm13.bed')
    >>> get_cytobands("invalid_genome")
    Traceback (most recent call last):
    ...
    ValueError: Unknown genome: invalid_genome
    """
    match genome:
        case GENOME.HG38:
            return STATIC_PATH / "cytobands_HG38.bed"
        case GENOME.CHM13:
            return STATIC_PATH / "cytobands_chm13.bed"
        case GENOME.HS1:
            return STATIC_PATH / "cytobands_chm13.bed"
        case _:
            msg = f"Unknown genome: {genome}"
            raise ValueError(msg)


def get_cytoband_df(genome: GENOME) -> pd.DataFrame:
    """
    Convert the cytogram file for the given genome into a dataframe.

    :param genome: The genome to plot the ideogram for.
    """
    cytobands = pd.read_csv(
        get_cytobands(genome), sep="\t", names=["chrom", "chromStart", "chromEnd", "name", "gieStain"]
    )
    cytobands["arm"] = cytobands["name"].str[0]
    cytobands["colour"] = cytobands["gieStain"].map(COLOUR_LOOKUP)
    cytobands["width"] = cytobands["chromEnd"] - cytobands["chromStart"]
    return cytobands


if __name__ == "__main__":
    fig, axes = plt.subplots(
        ncols=1,
        nrows=22,
        figsize=(11, 11),
        facecolor="white",
    )
    cytoband_df = get_cytoband_df(GENOME.CHM13)
    for ax, contig_name in zip(axes, range(1, 23)):
        chromosome = f"chr{contig_name}"
        plot_ideogram(ax, cytoband_df, chromosome)
    fig.savefig("ideogram.png", dpi=300)
