"""
Plotting function helpers for the ideogram plots.
"""
from typing import Optional

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.collections import BrokenBarHCollection
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path
from matplotlib.ticker import FuncFormatter


def set_xmargin(ax: Axes, left: float = 0.1, right: float = 0.0) -> None:
    """
    Adjust the x-axis margin of a given Axes object.

    :param ax: The Axes object to modify.
    :param left: Factor by which to expand the left x-axis limit.
    :param right: Factor by which to expand the right x-axis limit.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.set_xlim(0, 10)
    (0.0, 10.0)
    >>> set_xmargin(ax)
    >>> ax.get_xlim()
    (-1.0, 10.0)
    """

    # https://stackoverflow.com/a/49382894
    ax.set_xmargin(0)
    ax.autoscale_view()
    lim = ax.get_xlim()
    delta = lim[1] - lim[0]
    left = lim[0] - delta * left
    right = lim[1] + delta * right
    ax.set_xlim(left, right)


def set_ymargin(ax: Axes, top: float = 0.1, bottom: float = 0.0) -> None:
    """
    Adjust the y-axis margin of a given Axes object.

    :param ax: The Axes object to modify.
    :param top: Factor by which to expand the top y-axis limit.
    :param bottom: Factor by which to expand the bottom y-axis limit.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.set_ylim(0, 10)
    (0.0, 10.0)
    >>> set_ymargin(ax)
    >>> ax.get_ylim()
    (0.0, 11.0)
    """

    # https://stackoverflow.com/a/49382894
    ax.set_ymargin(0)
    ax.autoscale_view()
    lim = ax.get_ylim()
    delta = lim[1] - lim[0]
    bottom = lim[0] - delta * bottom
    top = lim[1] + delta * top
    ax.set_ylim(bottom, top)


def human_readable_yield(num: int, factor: int = 1000, suffix: str = "B") -> str:
    """
    Return a human readable string of a large number using SI unit prefixes.

    :param num: A number to convert to decimal form.
    :type num: int
    :param factor: The SI factor, use 1000 for SI units and 1024 for binary multiples.
    :type factor: int, optional
    :param suffix: The suffix to place after the SI prefix, for example use B for SI units and iB for binary multiples.
    :type suffix: str, optional
    :return: Returns the input number formatted to two decimal places with the SI unit and suffix.
    :rtype: str
    >>> human_readable_yield(7)
    '7.00 B'
    >>> human_readable_yield(1200)
    '1.20 kB'
    >>> human_readable_yield(1200, 1024, "iB")
    '1.17 kiB'
    >>> human_readable_yield(1500000)
    '1.50 MB'
    """
    for unit in ["", "k", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < factor:
            return f"{num:3.2f} {unit}{suffix}"
        num /= factor
    return "{n:3.2f} {u}{s}".format(n=num, u="Y", s=suffix)


@FuncFormatter
def format_genomics_lengths(x, _pos):
    return f"{human_readable_yield(x)}"


def plot_scale(ax: Axes, start: int, stop: int, scale_loc: float = 0.1, length_loc: float = 0.75):
    """
    Plot a scale bar on a given axis.

    :param ax: Matplotlib axis object to draw the scale on.
    :param start: Starting point of the scale bar.
    :param stop: Ending point of the scale bar.
    :param scale_loc: Location for the scale axis.
    :param length_loc: Location for the length of the scale bar.
    :return: Updated axis object with the scale bar.


    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax = plot_scale(ax, 0, 1000)
    >>> ax.get_xlim()  # To test if the xlim is set correctly
    (0.0, 1000.0)
    """
    ax.set_xlim(start, stop)
    ax.yaxis.set_visible(False)
    ax.set_ylim(0, 1)
    for side in ("right", "top", "left"):
        ax.spines[side].set_visible(False)
    ax.set_xlim(start, stop)
    ax.spines.bottom.set_position(("axes", scale_loc))

    scale = stop - start
    mid_point = start + scale / 2

    for point in (start, stop):
        ax.annotate(
            f"{human_readable_yield(scale)}",
            xytext=(mid_point, length_loc),
            xycoords="data",
            xy=(point, length_loc),
            textcoords="data",
            size="x-large",
            va="center",
            ha="center",
            arrowprops={"arrowstyle": "-|>", "facecolor": "black"},
        )

    ax.xaxis.set_major_formatter(format_genomics_lengths)
    return ax


def plot_bed_regions(ax: Axes, df: pd.DataFrame, target: str, start: int, stop: int) -> Axes:
    """
    Plot BED regions from a DataFrame on a given axis.
    Plot genomic regions in a given DataFrame onto a given Axes object.
    It takes in an Axes object, a pandas DataFrame containing genomic data, a target chromosome,
    and start and stop positions for the region of interest.
    It then filters the DataFrame to only include rows corresponding to the
    target chromosome and the specified region,
    and plots the resulting regions as rectangles on the given Axes object.
    Finally, it removes the spines and axes from the plot and returns the modified Axes object.

    :param ax: Matplotlib axis object where the regions will be plotted.
    :param df: DataFrame containing BED-like formatted data with columns "chromosome", "start", and "end".
    :param target: Target chromosome to filter and plot.
    :param start: Starting base pair position for the region of interest.
    :param stop: Ending base pair position for the region of interest.
    :return: Updated axis object with the plotted BED regions.


    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> data = {
    ...     "chromosome": ["chr1", "chr1", "chr2"],
    ...     "start": [100, 1500, 200],
    ...     "end": [500, 2000, 400]
    ... }
    >>> df = pd.DataFrame(data)
    >>> fig, ax = plt.subplots()
    >>> ax = plot_bed_regions(ax, df, "chr1", 0, 2500)
    >>> ax.get_xlim()  # To test if the regions were plotted correctly (not a direct measure but gives an idea)
    (0.0, 1.0)
    """
    df = df[df["chromosome"].eq(target) & (df["start"].between(start, stop) | df["end"].between(start, stop))]
    ylim = ax.get_ylim()
    for _, start_, stop_ in df.itertuples(False):
        r = Rectangle(
            (start_, min(ylim)),
            width=stop_ - start_,
            height=max(ylim),
            fill=True,
        )
        ax.add_patch(r)

    for side in ("right", "top", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    return ax


def plot_ideogram(
    ax: Axes,
    cytobands_df: pd.DataFrame,
    target: str,
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
    df = cytobands_df[cytobands_df["chrom"].eq(target)]
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
