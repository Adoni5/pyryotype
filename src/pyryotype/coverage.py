import pandas as pd
from matplotlib.axes import Axes
from matplotlib.ticker import NullFormatter

from pyryotype.plotting_utils import set_ymargin


def _invert_regions(start, stop, regions):
    """
    Return inverted regions within a specified interval.

    Given an overall region defined by `start` and `stop`, and a list of sub-regions
    within this overall region, this function returns a list of inverted regions
    that are not covered by the sub-regions.

    Parameters:
    - start (int): Start of the overall region.
    - stop (int): End of the overall region.
    - regions (List[Tuple[int, int]]): List of sub-regions as (start, stop) tuples.

    Returns:
    - List[Tuple[int, int]]: List of inverted regions as (start, stop) tuples.

    Examples:

    >>> _invert_regions(0, 100, [(40, 50), (20, 30), (70, 80)])
    [(0, 20), (30, 40), (50, 70), (80, 100)]

    >>> _invert_regions(0, 50, [(10, 20), (30, 40)])
    [(0, 10), (20, 30), (40, 50)]

    >>> _invert_regions(0, 100, [])
    [(0, 100)]

    >>> _invert_regions(0, 100, [(0, 50), (50, 100)])
    []
    """
    # Sort regions by start
    regions = sorted(regions, key=lambda x: x[0])

    # Add the overall start and stop to the regions for easy calculation
    regions = [(start, start), *regions, (stop, stop)]

    inverted = []
    for i in range(len(regions) - 1):
        inverted_start = regions[i][1]
        inverted_stop = regions[i + 1][0]
        if inverted_start != inverted_stop:
            inverted.append((inverted_start, inverted_stop))
    return inverted


def plot_coverage(
    ax: Axes,
    df: pd.DataFrame,
    target: str,
    regions: list[tuple[int, int]] | None = None,
    start: int | None = None,
    stop: int | None = None,
    ylabel: str | None = "Coverage",
    top_margin=0.1,
    bottom_margin=0,
    **kwargs,
):
    """
    Plot coverage data on a given axis.

    This function plots coverage data from a provided DataFrame between specified start and stop positions
    on a specific target chromosome.

    :param ax: The matplotlib axis object to plot on.
    :param df: The DataFrame containing coverage data. Expected columns: "chromosome", "start", "end", "value".
    :param target: The chromosome target to filter the DataFrame by.
    :param regions: The regions to highlight on the plot.
    :param start: The starting position to filter the DataFrame by.
    :param stop: The ending position to filter the DataFrame by.
    :param ylabel: The label for the y-axis.
    :param top_margin: The top margin for y-axis, default is 0.1.
    :param bottom_margin: The bottom margin for y-axis, default is 0.
    :param kwargs: Additional keyword arguments to pass to ax.fill_between().

    :return: The matplotlib axis object with the plotted data.

    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> df = pd.DataFrame({
    ...     'chromosome': ['chr1', 'chr1', 'chr2'],
    ...     'start': [1, 2, 1],
    ...     'end': [2, 3, 2],
    ...     'value': [5, 10, 15]
    ... })
    >>> fig, ax = plt.subplots()
    >>> plot_coverage(ax, df, 'chr1' )
    <...>
    >>> plot_coverage(ax, df, 'chr1', start=2)
    <...>
    >>> plot_coverage(ax, df, 'chr1', stop=2)
    <...>
    >>> plot_coverage(ax, df, 'chr1', regions=[(1, 2)])
    <...>
    >>> plot_coverage(ax, df, 'chr1', regions=[(1, 2), (2, 3)])
    <...>
    >>> plot_coverage(ax, df, 'chr1', regions=[(1, 2), (2, 3)], start=2, stop=3)
    <...>
    >>> plot_coverage(ax, df, 'chr1', ylabel='Read Depth')
    <...>
    """
    df = df[df["chromosome"].eq(target)]
    if start is not None:
        df = df[df["start"].ge(start)]
    if stop is not None:
        df = df[df["end"].le(stop)]
    start = start if start is not None else df["start"].min()
    stop = stop if stop is not None else df["end"].max()

    if regions:
        inverted_regions = _invert_regions(start, stop, regions)
        for start, stop in inverted_regions:
            df_region = df[(df["start"].between(start, stop) | df["end"].between(start, stop))]
            ax.fill_between(df_region["start"], df_region["value"], step="post", interpolate=True, color="C1")
        for start, stop in regions:
            df_region = df[(df["start"].between(start, stop) | df["end"].between(start, stop))]
            ax.fill_between(df_region["start"], df_region["value"], step="post", interpolate=True, color="C0")
    else:
        ax.fill_between(
            df["start"],
            df["value"],
            step="post",
            interpolate=True,
            **kwargs,
        )
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_tick_params(length=0)
    set_ymargin(ax, top=top_margin, bottom=bottom_margin)
    ax.grid(True, which="both", axis="both")
    ax.set_ylabel(ylabel, size="x-large")
    ax.yaxis.set_tick_params(labelsize="x-large")
    return ax
