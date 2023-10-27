"""Ideogram plotting executables."""
from enum import Enum
from pathlib import Path
from typing import Optional

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.collections import BrokenBarHCollection
from matplotlib.patches import PathPatch, Rectangle

from ideogram.plotting_utils import set_xmargin


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


def plot_ideogram(
    ax: Axes,
    target: str,
    genome: GENOME = GENOME.HG38,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    lower_anchor: int = 0,
    height: int = 1,
    curve: float = 0.05,
    y_margin: float = 0.05,
    right_margin: float = 0.005,
    left_margin: float = 0.25,
    target_region_extent: float = 0.3,
):
    """
    Plot a chromosome ideogram with cytobands and optionally highlight a specific region.

    :param ax: Matplotlib axis object where the ideogram will be plotted.
    :param cytobands_df: DataFrame containing cytoband data with columns "chrom", "chromStart",
      "chromEnd", "gieStain", and "colour".
    :param target: Target chromosome to filter and plot.
    :param start: Starting base pair position for the region of interest (optional).
    :param stop: Ending base pair position for the region of interest (optional).
    :param lower_anchor: Lower anchor point for the ideogram, for outline.
    :param height: Height of the ideogram.
    :param curve: Curve factor for the ideogram edges.
    :param y_margin: Margin for the y-axis.
    :param right_margin: Margin for the right side of the x-axis.
    :param left_margin: Margin for the left side of the x-axis.
    :param target_region_extent: Extent of the target region highlight.

    :return: Updated axis object with the plotted ideogram.

    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> data = {
    ...     "chrom": ["chr1", "chr1", "chr1"],
    ...     "chromStart": [0, 100, 200],
    ...     "chromEnd": [100, 200, 300],
    ...     "gieStain": ["p", "cen", "q"],
    ...     "colour": ["#FFFFFF", "#FF0000", "#FFFFFF"],
    ...     "width": [100, 100, 100]
    ... }
    >>> df = pd.DataFrame(data)
    >>> fig, ax = plt.subplots()
    >>> ax = plot_ideogram(ax, df, "chr1", start=50, stop=250)
    >>> ax.get_xlim()  # To test if the ideogram was plotted (not a direct measure but gives an idea)
    (0.0, 1.0)
    """
    # TODO: various kwds params for passing through to other methods
    df = get_cytoband_df(genome)
    df = df[df["chrom"].eq(target)]
    yrange = (lower_anchor, height)  # lower anchor, height
    ymid = (max(yrange) - min(yrange)) / 2
    xrange = df[["chromStart", "width"]].values
    chr_len = df["chromEnd"].max()

    # Stain lines
    bars = BrokenBarHCollection(xrange, yrange, facecolors=df["colour"])
    ax.add_collection(bars)

    # Define and draw the centromere using the rows marked as 'cen' in the 'gieStain' column
    cen_df = df[df["gieStain"].str.contains("cen")]
    cen_start = cen_df["chromStart"].min()
    cen_end = cen_df["chromEnd"].max()

    cen_poly = [
        (cen_start, lower_anchor),
        (cen_start, height),
        (cen_end, lower_anchor),
        (cen_end, height),
    ]
    cen_patch = PathPatch(Path(cen_poly), facecolor=(0.8, 0.4, 0.4), lw=0)
    ax.add_patch(cen_patch)
    # Define and draw the chromosome outline, taking into account the shape around the centromere
    chr_move, chr_poly = zip(
        *[
            (Path.MOVETO, (lower_anchor, height)),
            # Top left, bottom right: ‾\_
            (Path.LINETO, (cen_start, height)),
            (Path.LINETO, (cen_end, lower_anchor)),
            (Path.LINETO, (chr_len, lower_anchor)),
            # Right telomere: )
            (Path.LINETO, (chr_len, lower_anchor)),
            (Path.CURVE3, (chr_len + chr_len * curve, ymid)),
            (Path.LINETO, (chr_len, height)),
            # Top right, bottom left: _/‾
            (Path.LINETO, (cen_end, height)),
            (Path.LINETO, (cen_start, lower_anchor)),
            (Path.LINETO, (lower_anchor, lower_anchor)),
            # Left telomere: (
            (Path.CURVE3, (lower_anchor - chr_len * curve, ymid)),
            (Path.LINETO, (lower_anchor, height)),
            (Path.CLOSEPOLY, (lower_anchor, height)),
        ]
    )
    chr_patch = PathPatch(Path(chr_poly, chr_move), fill=None, joinstyle="round")
    ax.add_patch(chr_patch)

    # If start and stop positions are provided, draw a rectangle to highlight this region
    if start is not None and stop is not None:
        r = Rectangle(
            (start, lower_anchor - target_region_extent),
            width=stop - start,
            height=height + 2 * target_region_extent,
            fill=False,
            edgecolor="r",
            linewidth=1,
            joinstyle="round",
        )
        ax.add_patch(r)

    # Adjust x-axis margins
    set_xmargin(ax, left=left_margin, right=right_margin)
    ax.set_ymargin(y_margin)

    # Remove axis spines and ticks for a cleaner look
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    # Add chromosome name to the plot
    x0, x1 = ax.get_xlim()
    name = f"Chromosome {target.lstrip('chr')}"
    ax.text(x0, ymid, name, fontsize="x-large", va="center")

    return ax


if __name__ == "__main__":
    fig, axes = plt.subplots(
        ncols=1,
        nrows=22,
        figsize=(11, 11),
        facecolor="white",
    )
    genome = GENOME.CHM13
    for ax, contig_name in zip(axes, range(1, 23)):
        chromosome = f"chr{contig_name}"
        plot_ideogram(ax, target=chromosome, genome=genome)
    fig.savefig("ideogram.png", dpi=300)
