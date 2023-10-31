"""
Plotting function helpers for the ideogram plots.
"""

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
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
    >>> human_readable_yield(1e26)
    '100.00 YB'
    """
    for unit in ["", "k", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < factor:
            return f"{num:3.2f} {unit}{suffix}"
        num /= factor
    return "{n:3.2f} {u}{s}".format(n=num, u="Y", s=suffix)


@FuncFormatter
def format_genomics_lengths(x: int | float, _pos) -> str:
    """
    This function formats genomic lengths into a human-readable format.

    :param x: The genomic length to be formatted.
    :param _pos: A placeholder parameter required by the FuncFormatter decorator. It's not used in this function.
    :return: Returns the input genomic length formatted into a human-readable string.

    .. doctest::

        >>> format_genomics_lengths(7, 0)
        '7.00 B'
        >>> format_genomics_lengths(1200, 0)
        '1.20 kB'
        >>> format_genomics_lengths(1500000, 0)
        '1.50 MB'
        >>> format_genomics_lengths(1e25, 0)
        '10.00 YB'
    """
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
