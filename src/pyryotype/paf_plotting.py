import secrets
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from enum import Enum
from typing import Protocol

import matplotlib.pyplot as plt  # noqa: F401
from matplotlib import patches
from matplotlib.axes import Axes

FIELDS = [
    "query_name",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target_name",
    "target_length",
    "target_start",
    "target_end",
    "residue_matches",
    "alignment_block_length",
    "mapping_quality",
    "tags",
]
NA_VALUES = ["*"]
INTERVAL_MIN_SIZE = 2


class PAFProtocol(Protocol):
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    residue_matches: int
    alignment_block_length: int
    mapping_quality: int
    tags: dict[str, str]

    def __str__(self) -> str:
        ...

    def _fmt_tags(self) -> str:
        ...

    def blast_identity(self) -> float:
        ...


@dataclass
class PAF(PAFProtocol):
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    residue_matches: int
    alignment_block_length: int
    mapping_quality: int
    tags: dict[str, str]

    def __str__(self):
        """Formats a record as a PAF line for writing to a file"""
        return "{}\t{}".format(
            "\t".join(
                map(
                    str,
                    [
                        self.query_name,
                        self.query_length,
                        self.query_start,
                        self.query_end,
                        self.strand,
                        self.target_name,
                        self.target_length,
                        self.target_start,
                        self.target_end,
                        self.residue_matches,
                        self.alignment_block_length,
                        self.mapping_quality,
                    ],
                )
            ),
            self._fmt_tags(),
        )

    def _fmt_tags(self):
        """Format tag dict as SAM style"""
        return "\t".join("{}:{}:{}".format(k, "Z", v) for k, v in self.tags.items())

    def blast_identity(self):
        """BLAST identity, see:
        https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
        """
        return self.residue_matches / self.alignment_block_length

    @classmethod
    def from_protocol(cls, obj: PAFProtocol) -> "PAF":
        if isinstance(obj, tuple):  # assuming namedtuple is essentially a tuple
            return cls(**obj._asdict())
        elif isinstance(obj, cls):
            return cls(**asdict(obj))
        else:
            msg = "Unsupported type"
            raise TypeError(msg)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """
    Merge the given intervals, returning collapsed intervals.
    Assumes the intervals are sorted.

    Parameters:
    - intervals (list[tuple[int, int]]): A list of sorted intervals represented as tuples.

    Returns:
    - list[tuple[int, int]]: A list of merged intervals.

    Examples:
    >>> merge_intervals([])
    []

    >>> merge_intervals([(1,3)])
    [(1, 3)]

    >>> merge_intervals([(1, 3), (2, 4)])
    [(1, 4)]

    >>> merge_intervals([(1, 2), (3, 4)])
    [(1, 2), (3, 4)]

    >>> merge_intervals([(1, 4), (2, 3)])
    [(1, 4)]

    >>> merge_intervals([(1, 5), (6, 10), (7, 8)])
    [(1, 5), (6, 10)]

    >>> merge_intervals([(1, 5), (5, 10)])
    [(1, 10)]

    Note:
    - It's crucial to pass the intervals sorted for the function to work correctly.
    """
    if len(intervals) < INTERVAL_MIN_SIZE:
        return intervals

    collapsed_intervals = []
    (curr_start, curr_end) = intervals[0]

    for start, end in intervals[1:]:
        if start > curr_end:  # We have a new non-overlapping start
            collapsed_intervals.append((curr_start, curr_end))
            curr_start, curr_end = start, end
        else:  # Start is within the current range
            curr_end = max(curr_end, end)

    collapsed_intervals.append((curr_start, curr_end))
    return collapsed_intervals


def _collapse_multiple_mappings(alignments: Iterator[PAFProtocol]) -> Iterator[PAFProtocol]:
    """
    Collapse multiple mappings of sequences into a single mapping by selecting the mapping
    that covers the largest proportion of the sequence.

    This function takes an iterator of alignments and yields the largest proportion alignments
    for each query sequence.


    :param alignments: An iterator of alignments, typically sorted by query_name.

    :yields PAFProtocol: The alignment with the largest mapping proportion for each query sequence.

    Example:
    >>> alignment1 = PAF(query_name='seq1', query_length=120, query_start=10, query_end=40, strand='+',
    ...                  target_name='targetA', target_length=200, target_start=50, target_end=80,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignment2 = PAF(query_name='seq1', query_length=120, query_start=20, query_end=50, strand='+',
    ...                  target_name='targetB', target_length=250, target_start=60, target_end=90,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignment3 = PAF(query_name='seq1', query_length=120, query_start=90, query_end=120, strand='+',
    ...                  target_name='targetA', target_length=200, target_start=130, target_end=160,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignment4 = PAF(query_name='seq1', query_length=120, query_start=90, query_end=120, strand='-',
    ...                  target_name='targetA', target_length=200, target_start=130, target_end=160,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignment5 = PAF(query_name='seq2', query_length=120, query_start=90, query_end=120, strand='+',
    ...                  target_name='targetA', target_length=200, target_start=130, target_end=160,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignments = [alignment1, alignment2, alignment3]
    >>> result = list(_collapse_multiple_mappings(iter(alignments)))
    >>> len(result)
    1
    >>> result[0].target_name
    'targetA'
    >>> result[0].target_start
    50
    >>> result[0].target_end
    160

    # Now add 2 contigs (seq1 and seq2) and different strands on seq1
    >>> alignments = [alignment1, alignment2, alignment3, alignment4, alignment5]
    >>> result = list(_collapse_multiple_mappings(iter(alignments)))
    >>> len(result)
    2

    Note:
    The above example demonstrates that for the 'seq1' query, the alignment to 'targetA' is chosen
    because it covers a larger proportion of the query sequence.
    """
    curr_id = None

    #  Store the proportion of the contig that maps to a given contig
    contig_mapping_proportion = defaultdict(list)

    for alignment in alignments:
        if curr_id is None:
            curr_id = alignment.query_name
        if curr_id != alignment.query_name:
            # We have a new query, yield the previous one
            # First we check which contig/strand had the largest amount of alignment to it
            # Then we yield the largest contig/strand
            max_key = max(
                contig_mapping_proportion,
                key=lambda k: sum(item.query_end - item.query_start for item in contig_mapping_proportion[k]),
            )
            largest_proportion_alignments = sorted(contig_mapping_proportion[max_key], key=lambda x: x.target_start)

            paf = PAF.from_protocol(largest_proportion_alignments[0])
            paf.target_end = largest_proportion_alignments[-1].target_end
            yield paf
            curr_id = alignment.query_name
            contig_mapping_proportion.clear()

        contig_mapping_proportion[(alignment.target_name, alignment.strand)].append(alignment)
    # Yield the last contig alignment
    max_key = max(
        contig_mapping_proportion,
        key=lambda k: sum(item.query_end - item.query_start for item in contig_mapping_proportion[k]),
    )
    largest_proportion_alignments = sorted(contig_mapping_proportion[max_key], key=lambda x: x.target_start)

    paf = PAF.from_protocol(largest_proportion_alignments[0])
    paf.target_end = largest_proportion_alignments[-1].target_end
    yield paf


def generate_random_color():
    """
    Generates a random RGB color.

    :return: A tuple of three integers representing RGB values.

    Example:
    >>> color = generate_random_color()
    >>> len(color)
    3
    >>> all(0 <= val <= 255 for val in color)
    True
    """
    return (secrets.randbelow(255), secrets.randbelow(255), secrets.randbelow(255))


class PlotMode(Enum):
    """Plotting modes for PAF alignments"""

    STRICT = 0
    CHILL = 1
    UNIQUE_COLOURS = 2
    STRAND_COLOURS = 3


def plot_paf_alignments(
    ax: Axes,
    alignments: Iterator[PAFProtocol],
    target: str,
    mapq_filter: int = 0,
    strict: PlotMode = PlotMode.CHILL,
    contig_colours: PlotMode = PlotMode.STRAND_COLOURS,
) -> Axes:
    """
    Plots Pairwise Alignment Format (PAF) alignments as rectangles on a matplotlib axis.

    :param ax: Matplotlib axis
    :param alignments: Iterator provided by parsepaf of PAF named tuples
    :param target: The target contig name to plot
    :param filter: Map q filter, filter out any alignments under this threshold. Defaults to 0 (All alignments).
    :param strict: If STRICT, Collapse contigs with multiple primary mappings into one
      contiguous block. Defaults to CHILL, where multiple mappings are plotted as separate blocks.
    :contig_colours: If UNIQUE_COLOURS, use a unique colour for each contig. Defaults to STRAND_COLOURS,
      where contigs coloured by strand alignment

    :return: None

    :note:
        Mapq filter is applied before collapsing multiple mappings, if strict is True.

    Example:
    >>> fig, ax = plt.subplots()
    >>> alignments = [PAF(query_name='seq1', query_length=120, query_start=10, query_end=40, strand='+',
    ...                  target_name='chr1', target_length=200, target_start=50, target_end=80,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})]
    >>> ax = plot_paf_alignments(ax, alignments, target="chr1")
    >>> ax.get_xlim()
    (0.0, 200.0)
    >>> ax.get_ylim()
    (0.0, 1.0)

    # Now add 2 contigs (seq1 and seq2) and set strict to false ( collapse contigs with multiple mappings )
    >>> import matplotlib.pyplot as plt
    >>> from collections import namedtuple

    >>> PAF = namedtuple('PAF', ['query_name', 'query_length', 'query_start', 'query_end', 'strand', 'target_name',
    ...                          'target_length', 'target_start', 'target_end', 'residue_matches',
    ...                          'alignment_block_length', 'mapping_quality', 'tags'])

    >>> fig, ax = plt.subplots()
    >>> alignment1 = PAF(query_name='seq1', query_length=120, query_start=10, query_end=40, strand='+',
    ...                  target_name='targetA', target_length=200, target_start=50, target_end=80,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignment2 = PAF(query_name='seq1', query_length=120, query_start=20, query_end=50, strand='+',
    ...                  target_name='targetB', target_length=250, target_start=60, target_end=90,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignment3 = PAF(query_name='seq1', query_length=120, query_start=90, query_end=120, strand='+',
    ...                  target_name='targetA', target_length=200, target_start=130, target_end=160,
    ...                  residue_matches=30, alignment_block_length=30, mapping_quality=60, tags={})
    >>> alignments = [alignment1, alignment2, alignment3]
    >>> ax = plot_paf_alignments(ax, alignments, target="targetA", strict=PlotMode.STRICT)
    >>> ax.get_xlim()
    (0.0, 200.0)
    >>> ax.get_ylim()
    (0.0, 1.0)
    """

    iterable = filter(
        lambda x: x.target_name == target and x.mapping_quality >= mapq_filter,
        alignments,
    )
    if strict:
        iterable = _collapse_multiple_mappings(iterable)
    for alignment in iterable:
        if contig_colours == PlotMode.UNIQUE_COLOURS:
            colour = generate_random_color()
        else:
            colour = "#332288" if alignment.strand == "+" else "#882255"

        target_len = alignment.target_length
        target_rect = patches.Rectangle(
            (alignment.target_start, 0),
            alignment.target_end - alignment.target_start,
            0.9,
            edgecolor=None,
            facecolor=colour,
            fill=True,
            visible=True,
        )
        ax.add_patch(target_rect)

    ax.set_xlim((0, target_len))
    ax.set_ylim((0, 1))
    ax.set_xlabel("Position")
    # ax.legend()
    return ax
