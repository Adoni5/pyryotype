from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import EngFormatter
from pyryotype.coverage import plot_coverage


def test_plot_coverage():
    test_cov_file = Path(__file__).parent / "static" / "chr1_cov.regions.bed.gz"

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(11, 2),
    )
    df = pd.read_csv(test_cov_file, sep="\t", names=["chromosome", "start", "end", "value"])
    ax = plot_coverage(ax, df, "chr1", regions=[(0, 100000000)], ylabel="Coverage", color="black")
    ax.set_xlabel("Genomic Position (bp)")

    # ax.yaxis.set_visible(False)
    ax.set_yticks([])
    ax.set_yticklabels([])
    for side in ("right", "top", "left"):
        ax.spines[side].set_visible(False)
    formatter = EngFormatter(unit="b", places=1)
    # Set formatter for the x-axis
    ax.xaxis.set_major_formatter(formatter)

    fig.savefig("example_outputs/test_coverage.png", dpi=300, bbox_inches="tight")
