from enum import Enum
from pathlib import Path

# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.collections import BrokenBarHCollection, PolyCollection
# from matplotlib.patches import PathPatch, Rectangle
# from matplotlib.path import Path
# from matplotlib.ticker import FuncFormatter, NullFormatter


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
    'cytoBand_HG38.tsv'
    >>> get_cytobands(GENOME.CHM13)
    'cytoBand_CHM13.tsv'
    >>> get_cytobands(GENOME.HS1)
    'cytoBand_HS1.tsv'
    >>> get_cytobands("invalid_genome")
    Traceback (most recent call last):
    ...
    ValueError: Unknown genome: invalid_genome
    """
    match genome:
        case GENOME.HG38:
            return STATIC_PATH / "cytoband_HG38.bed"
        case GENOME.CHM13:
            return STATIC_PATH / "cytoband_chm13.bed"
        case GENOME.HS1:
            return STATIC_PATH / "cytoband_chm13.bed"
        case _:
            msg = f"Unknown genome: {genome}"
            raise ValueError(msg)
