import unittest
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.ticker import EngFormatter
from pyryotype.paf_plotting import PAFProtocol, PlotMode, plot_paf_alignments
from readpaf import parse_paf


class TestPAFProtocol(unittest.TestCase):
    class PAFProtocolImpl(PAFProtocol):
        def __str__(self) -> str:
            return "Test String"

        def _fmt_tags(self) -> str:
            return "Test Tags"

        def blast_identity(self) -> float:
            return 0.8

    def test_protocol_methods(self):
        paf = self.PAFProtocolImpl()
        self.assertEqual(str(paf), "Test String")
        self.assertEqual(paf._fmt_tags(), "Test Tags")
        self.assertEqual(paf.blast_identity(), 0.8)


def test_drawing_paf_alignments_strict():
    test_paf = Path(__file__).parent / "static" / "test.paf"
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(11, 1),
    )

    ax = plot_paf_alignments(
        ax,
        parse_paf(test_paf.open()),
        target="chr1",
        mapq_filter=0,
        strict=PlotMode.STRICT,
        contig_colours=PlotMode.UNIQUE_COLOURS,
    )
    ax.set_xlabel("")

    # ax.yaxis.set_visible(False)
    ax.set_yticks([])
    ax.set_yticklabels([])
    for side in ("right", "top", "left"):
        ax.spines[side].set_visible(False)
    formatter = EngFormatter(unit="b", places=1)
    # Set formatter for the x-axis
    ax.xaxis.set_major_formatter(formatter)
    fig.savefig("example_outputs/test_paf_plotting.png", dpi=300, bbox_inches="tight")


def test_drawing_paf_alignments_chill():
    test_paf = Path(__file__).parent / "static" / "test.paf"
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(11, 1),
    )

    ax = plot_paf_alignments(
        ax,
        parse_paf(test_paf.open()),
        target="chr1",
        mapq_filter=0,
        strict=PlotMode.CHILL,
        contig_colours=PlotMode.UNIQUE_COLOURS,
    )
    ax.set_xlabel("")

    # ax.yaxis.set_visible(False)
    ax.set_yticks([])
    ax.set_yticklabels([])
    for side in ("right", "top", "left"):
        ax.spines[side].set_visible(False)
    fig.savefig("example_outputs/test_paf_plotting_chill.png", dpi=300, bbox_inches="tight")
