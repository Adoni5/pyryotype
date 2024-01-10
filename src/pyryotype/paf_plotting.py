"""
Module to manage and process Pairwise Alignment Format (PAF) data.

This module provides utilities to:
- Parse PAF formatted data and represent it using data classes.
- Perform various operations on PAF data such as collapsing multiple mappings, computing BLAST identities,
  and merging overlapping intervals.
- Visualize PAF alignments along a given target contig using Matplotlib.

Constants:
- FIELDS: A list of field names in the PAF format.
- NA_VALUES: A list of not available values in the PAF format.
- INTERVAL_MIN_SIZE: The minimum size of an interval to be considered for merging.

Classes:
- PAFProtocol: A protocol defining the attributes and methods for a PAF record.
- PAF: A dataclass that implements the PAFProtocol and represents a single PAF record.

Functions:
- merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]: Merge overlapping intervals.
- _collapse_multiple_mappings(alignments: Iterator[PAFProtocol]) -> Iterator[PAFProtocol]:
 Collapse multiple mappings into one.
- generate_random_color() -> tuple[int, int, int]: Generate a random RGB color.
- plot_paf_alignments(ax: Axes, alignments: Iterator[PAFProtocol],
target: str, mapq_filter: int = 0, strict=False, unique_contig_colours=False) -> Axes:
  Plot PAF alignments on a given matplotlib axis.

Examples and usage details are provided in the docstrings of the respective functions and classes.
"""

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
CHEVRON_FORWARD = ">"
CHEVRON_REVERSE = "<"
CHEVRON_FONTSIZE = 20
MAX_TRACKS = 10
LABEL_FONTSIZE = 11


class PAFProtocol(Protocol):
    """
    This class represents the PAF (Pairwise mApping Format) protocol. It is used to store alignment information
    between sequences in bioinformatics.

    :ivar query_name: The name of the query sequence.
    :ivar query_length: The length of the query sequence.
    :ivar query_start: The start position of the query in the alignment.
    :ivar query_end: The end position of the query in the alignment.
    :ivar strand: The strand of the target sequence.
    :ivar target_name: The name of the target sequence.
    :ivar target_length: The length of the target sequence.
    :ivar target_start: The start position of the target in the alignment.
    :ivar target_end: The end position of the target in the alignment.
    :ivar residue_matches: The number of residue matches in the alignment.
    :ivar alignment_block_length: The length of the alignment block.
    :ivar mapping_quality: The mapping quality score.
    :ivar tags: A dictionary of tag-value pairs.

    :method __str__: Returns a string representation of the PAFProtocol object.
    :method _fmt_tags: Returns a formatted string of the tags.
    :method blast_identity: Returns the BLAST identity of the alignment.

    """

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
    """
    A dataclass that represents a PAF (Pairwise mApping Format) record.

    :param query_name: The name of the query sequence.
    :param query_length: The total length of the query sequence.
    :param query_start: Start position of alignment in query.
    :param query_end: End position of alignment in query.
    :param strand: Relative strand: "+" or "-".
    :param target_name: The name of the target sequence.
    :param target_length: The total length of the target sequence.
    :param target_start: Start position of alignment in target.
    :param target_end: End position of alignment in target.
    :param residue_matches: Number of residue matches.
    :param alignment_block_length: Block length of the alignment.
    :param mapping_quality: Mapping quality score.
    :param tags: A dictionary of tags in SAM format.

    .. doctest::

       >>> paf = PAF(
       ...     query_name="seq1", query_length=1000, query_start=0, query_end=100,
       ...     strand="+", target_name="seq2", target_length=2000, target_start=0,
       ...     target_end=100, residue_matches=90, alignment_block_length=100,
       ...     mapping_quality=60, tags={"NM": "10", "MD": "100"}
       ... )
       >>> str(paf).replace("\\t", " ")
       'seq1 1000 0 100 + seq2 2000 0 100 90 100 60 NM:Z:10 MD:Z:100'
       >>> paf.blast_identity()
       0.9
       >>> paf._fmt_tags().replace("\\t", " ")
       'NM:Z:10 MD:Z:100'
       >>> PAF.from_protocol(["hello"])
       Traceback (most recent call last):
       TypeError: Unsupported type
    """

    # ... [rest of the class as you've written]

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


def filter_extraneous_mappings(alignments: Iterator[PAFProtocol]) -> Iterator[PAFProtocol]:
    """
    Iterates through PAF alignments and yields each alignment that is not fully contained within another.

    :param alignments: An iterator over PAF alignment records.

    :yields: Each alignment that is not fully contained within another alignment.

    Examples:
    >>> from collections import namedtuple
    >>> PAF = namedtuple('PAF', ['query_name', 'query_length', 'query_start', 'query_end', 'strand',
    ... 'target_name', 'target_length', 'target_start', 'target_end', 'residue_matches',
    ... 'alignment_block_length', 'mapping_quality', 'tags'])

    # Simple case: Two non-overlapping alignments
    >>> alignments = iter([PAF('seq1', 120, 10, 100, '+', 'chr1', 200, 50, 80, 30, 70, 60, {}),
    ... PAF('seq2', 120, 10, 100, '+', 'chr1', 200, 90, 120, 30, 70, 60, {})])
    >>> list(filter_extraneous_mappings(alignments))
    [PAF(query_name='seq1', ...), PAF(query_name='seq2', ...)]

    # Case where the second alignment is contained within the first
    >>> alignments = iter([PAF('seq1', 120, 10, 100, '+', 'chr1', 200, 50, 120, 30, 70, 60, {}),
    ... PAF('seq2', 120, 10, 50, '+', 'chr1', 200, 70, 90, 20, 40, 60, {})])
    >>> list(filter_extraneous_mappings(alignments))
    [PAF(query_name='seq1', ...)]

    # Case with multiple alignments, some overlapping, some contained
    >>> alignments = iter([PAF('seq1', 120, 10, 50, '+', 'chr1', 200, 50, 80, 30, 70, 60, {}),
    ... PAF('seq2', 120, 20, 70, '+', 'chr1', 200, 60, 90, 30, 70, 60, {}),
    ... PAF('seq3', 120, 30, 40, '+', 'chr1', 200, 70, 100, 20, 70, 60, {})])
    >>> list(filter_extraneous_mappings(alignments))
    [PAF(query_name='seq1', ...), PAF(query_name='seq2', ...), PAF(query_name='seq3', ...)]

    # More complex case with multiple overlaps and contained alignments
    >>> alignments = iter([PAF('seq1', 120, 10, 50, '+', 'chr1', 200, 50, 80, 30, 70, 60, {}),
    ... PAF('seq2', 120, 60, 100, '+', 'chr1', 200, 90, 120, 30, 70, 60, {}),
    ... PAF('seq3', 120, 20, 40, '+', 'chr1', 200, 70, 90, 20, 70, 60, {}),
    ... PAF('seq4', 120, 110, 150, '+', 'chr1', 200, 130, 160, 30, 70, 60, {})])
    >>> list(filter_extraneous_mappings(alignments))
    [PAF(query_name='seq1', ...), PAF(query_name='seq2', ...), PAF(query_name='seq4', ...)]
    """

    accepted_mappings = []  # List to store mappings that are not fully contained within others

    for alignment in alignments:
        start, end = alignment.target_start, alignment.target_end

        # Check if this alignment is fully contained within any accepted mapping
        if any(start >= alignment.target_start and end <= alignment.target_end for alignment in accepted_mappings):
            continue  # Skip this alignment

        # Add this alignment to the list of accepted mappings
        accepted_mappings.append(alignment)
    yield from accepted_mappings


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
    Collapse multiple supplementary mappings on the same contig/strand combination into one contiguous block,
    and yield a PAF for each block.

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
    2
    >>> result[0].target_name
    'targetA'
    >>> result[0].target_start
    50
    >>> result[0].target_end
    160
    >>> result[0].strand
    '+'

    # Now add 2 contigs (seq1 and seq2) and different strands on seq1
    >>> alignments = [alignment1, alignment2, alignment3, alignment4, alignment5]
    >>> result = list(_collapse_multiple_mappings(iter(alignments)))
    >>> len(result)
    3

    Note:
    We get three alignments back from the final test, as there are three query name, contig combinations
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
            # Yield supplementary alignments collapsed as well
            for query_name, contig in contig_mapping_proportion:
                largest_proportion_alignments = sorted(
                    contig_mapping_proportion[(query_name, contig)], key=lambda x: x.target_start
                )
                paf = PAF.from_protocol(largest_proportion_alignments[0])
                paf.target_end = largest_proportion_alignments[-1].target_end
                yield paf
            curr_id = alignment.query_name
            contig_mapping_proportion.clear()

        contig_mapping_proportion[(alignment.query_name, alignment.target_name)].append(alignment)
    # Yield the last contig alignment
    for query_name, contig in contig_mapping_proportion:
        largest_proportion_alignments = sorted(
            contig_mapping_proportion[(query_name, contig)], key=lambda x: x.target_start
        )
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
    return (secrets.randbelow(255) / 255, secrets.randbelow(255) / 255, secrets.randbelow(255) / 255)


class PlotMode(Enum):
    """Plotting modes for PAF alignments

    :ivar STRICT:
        Collapse contigs with multiple primary mappings into one contiguous block.
    :ivar CHILL:
        Plot multiple mappings as separate blocks.
    :ivar UNIQUE_COLOURS:
        Use a unique colour for each aligned block.
    :ivar STRAND_COLOURS:
        Colour blocks by strand alignment.
    :ivar CHEVRON:
        Plot chevron arrows representing strandedness for alignments.
    :ivar FILTER_DOWN:
        Filter out alignments that are fully contained within another alignment.
    :ivar EXPAND_OVERLAPS:
        Expand overlapping alignments into their own unique track
    """

    STRICT = 0
    CHILL = 1
    UNIQUE_COLOURS = 2
    STRAND_COLOURS = 3
    CHEVRON = 4
    FILTER_DOWN = 5
    EXPAND_OVERLAPS = 6


def _choose_track_with_min_overlap(start: int, end: int, tracks: list[list[tuple[int, int]]]) -> int:
    """
    Chooses a track with minimum overlap for a new rectangle defined by start and end.

    If a track with no overlap is found, the rectangle is placed in this track. Otherwise,
    the track with the least overlap is chosen.

    Parameters:
    start (int): Start position of the new rectangle.
    end (int): End position of the new rectangle.
    tracks (list[list[tuple[int, int]]]): A list of tracks, each track being a list of rectangles
                                          represented as (start, end) tuples.

    Returns:
    int: The index of the chosen track.

    Examples:
    >>> tracks = [[(0, 10), (20, 30)], [(10, 15)], [(30, 40)]]
    >>> _choose_track_with_min_overlap(12, 22, tracks)
    2
    >>> tracks  # The new rectangle (12, 22) is placed in track index 2 as it has the zero overlap
    [[(0, 10), (20, 30)], [(10, 15)], [(30, 40), (12, 22)]]

    >>> tracks = [[(0, 10)], [(10, 20)], [(20, 30)]]
    >>> _choose_track_with_min_overlap(5, 25, tracks)
    0
    >>> tracks  # The new rectangle (5, 25) is placed in track 0 as it the first track with the least overlap
    [[(0, 10), (5, 25)], [(10, 20)], [(20, 30)]]
    """
    # Initialize tracks (each track is a list of (start, end) tuples)
    track_found = False
    # First, try to find a track without any overlap
    for track_index, track in enumerate(tracks):  # noqa: B007
        if not any(start < existing_end and end > existing_start for existing_start, existing_end in track):
            # Place the rectangle in this track without overlap
            track.append((start, end))
            track_found = True
            break
    # If no free track is found, select track with minimum overlap
    if not track_found:
        min_overlap = float("inf")
        selected_track = 0

        for track_index, track in enumerate(tracks):
            overlap = sum(
                max(0, min(end, existing_end) - max(start, existing_start)) for existing_start, existing_end in track
            )

            if overlap < min_overlap:
                min_overlap = overlap
                selected_track = track_index

        # Add to the track with minimum overlap
        tracks[selected_track].append((start, end))
        track_index = selected_track
    return track_index


def plot_paf_alignments(
    ax: Axes,
    alignments: Iterator[PAFProtocol],
    target: str,
    mapq_filter: int = 0,
    strict: PlotMode = PlotMode.CHILL,
    contig_colours: PlotMode = PlotMode.STRAND_COLOURS,
    chevron: PlotMode | None = None,
    filter_down: PlotMode | None = None,
    expand_overlaps: PlotMode | None = None,
    label_n_largest: int | None = None,
    **kwargs,
) -> Axes:
    """
    Plots Pairwise Alignment Format (PAF) alignments as rectangles on a matplotlib axis.

    :param ax: Matplotlib axis
    :param alignments: Iterator provided by parsepaf of PAF named tuples
    :param target: The target contig name to plot
    :param filter: Map q filter, filter out any alignments under this threshold. Defaults to 0 (All alignments).
    :param strict: If STRICT, Collapse contigs with multiple primary mappings into one
      contiguous block. Defaults to CHILL, where multiple mappings are plotted as separate blocks.
    :param contig_colours: If UNIQUE_COLOURS, use a unique colour for each contig. Defaults to STRAND_COLOURS,
      where contigs coloured by strand alignment
    :param chevron: If True, plot chevron arrows representing strandedness
      for alignments. Defaults to None, where chevrons are not plotted.
    :param filter_down: If True, filter out alignments that are fully contained within another alignment.
    :param label_n_largest: Label n alignments by size of aligned block with their query name
    :param **kwargs: Additional keyword arguments in order to override the default chevron
        symbol. Options include 'chevron_symbol', 'label_fontsize' and 'chevron_fontsize'.
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
    >>> import matplotlib
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
    >>> ax = plot_paf_alignments(ax, alignments, target="targetA", strict=PlotMode.UNIQUE_COLOURS)

    >>> fig, ax = plt.subplots(figsize=(4,1))
    >>> ax.set_xlim((0,500))
    (0.0, 500.0)
    >>> alignments = [
    ...     PAF(query_name='seq1', query_length=120, query_start=10, query_end=100, strand='+',
    ... target_name='chr1', target_length=200, target_start=50, target_end=150, residue_matches=30,
    ... alignment_block_length=70, mapping_quality=60, tags={}),
    ...     PAF(query_name='seq2', query_length=120, query_start=10, query_end=11, strand='-',
    ... target_name='chr1', target_length=1000, target_start=90, target_end=91, residue_matches=10,
    ... alignment_block_length=5, mapping_quality=60, tags={})
    ... ]
    >>> plot_paf_alignments(ax, iter(alignments), "chr1", chevron=PlotMode.CHEVRON)
    <Axes: xlabel='Position'>

    # Check chevron is drawn for seq1 (enough space)
    >>> text_objs = [obj for obj in ax.get_children() if isinstance(obj, matplotlib.text.Text)]
    >>> any(">" in obj.get_text() for obj in text_objs)
    True

    # Check chevron is not drawn for seq2 (not enough space)
    >>> any("<" in obj.get_text() for obj in text_objs)
    False
    """

    iterable = filter(
        lambda x: x.target_name == target and x.mapping_quality >= mapq_filter,
        alignments,
    )
    if strict == PlotMode.STRICT:
        iterable = _collapse_multiple_mappings(iterable)
    if kwargs.get("sorted", False):
        iterable = sorted(iterable, key=lambda x: x.target_start)
    if filter_down == PlotMode.FILTER_DOWN:
        iterable = filter_extraneous_mappings(iterable)

    # Upper limit on the track
    max_y = 0.9
    max_tracks = kwargs.get("max_tracks", MAX_TRACKS)
    # How much to step by
    dy = max_y / MAX_TRACKS if expand_overlaps == PlotMode.EXPAND_OVERLAPS else max_y
    # Track which start stop are in use for each track
    tracks = [[] for _ in range(max_tracks)]
    target_len = 1
    for alignment in iterable:
        if contig_colours == PlotMode.UNIQUE_COLOURS:
            colour = generate_random_color()
        else:
            colour = "#332288" if alignment.strand == "+" else "#882255"
        start, end = alignment.target_start, alignment.target_end
        target_len = alignment.target_length
        track_index = 0
        if expand_overlaps == PlotMode.EXPAND_OVERLAPS:
            track_index = _choose_track_with_min_overlap(start, end, tracks)
        target_rect = patches.Rectangle(
            (start, dy * track_index),
            end - start,
            dy,
            edgecolor=None,
            facecolor=colour,
            fill=True,
            visible=True,
            zorder=1,
        )
        ax.add_patch(target_rect)
        # Add chevron if there's enough space
        if chevron == PlotMode.CHEVRON:
            strand = alignment.strand
            chevron_symbol = ">" if strand == "+" else "<"
            fig = ax.get_figure()
            # Calculate the rectangle width in inches lol
            font_size_pt = kwargs.get("chevron_fontsize", CHEVRON_FONTSIZE)  # 10 points

            # Convert font size to inches (1 point = 1/72 inches)
            font_size_inch = font_size_pt / 72
            # print(fig.get_size_inches())
            # print(f"font_size_inch: {font_size_inch}")
            # print(fig.dpi_scale_trans.inverted().transform(ax.transData.transform((alignment.target_end, 0))))
            rect_width_inches = (
                fig.dpi_scale_trans.inverted().transform(ax.transData.transform((alignment.target_end, 0)))
                - fig.dpi_scale_trans.inverted().transform(ax.transData.transform((alignment.target_start, 0)))
            )[0]
            # print(f"rect_width_inches: {rect_width_inches}")
            # Check if there's enough space for the chevron
            if rect_width_inches > font_size_inch + 0.05:
                # print("drawing chevron")
                ax.text(
                    (alignment.target_start + alignment.target_end) / 2,  # X position (center of the rectangle)
                    0.45,  # Y position (roughly the middle of the rectangle in height)
                    chevron_symbol,
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=font_size_pt,
                    zorder=2,
                )
        if label_n_largest:
            fig = ax.get_figure()
            # Calculate the rectangle width in inches lol
            font_size_pt = kwargs.get("label_fontsize", LABEL_FONTSIZE)  # 10 points

            # Convert font size to inches (1 point = 1/72 inches)
            font_size_inch = font_size_pt / 72
            # print(fig.get_size_inches())
            # print(f"font_size_inch: {font_size_inch}")
            # print(fig.dpi_scale_trans.inverted().transform(ax.transData.transform((alignment.target_end, 0))))
            rect_width_inches = (
                fig.dpi_scale_trans.inverted().transform(ax.transData.transform((alignment.target_end, 0)))
                - fig.dpi_scale_trans.inverted().transform(ax.transData.transform((alignment.target_start, 0)))
            )[0]
            # print(f"rect_width_inches: {rect_width_inches}")
            # Check if there's enough space for the chevron
            if rect_width_inches > font_size_inch * len(alignment.query_name) + 0.05:
                ax.text(
                    (alignment.target_start + alignment.target_end) / 2,  # X position (center of the rectangle)
                    0.9,  # Y position (roughly the middle of the rectangle in height)
                    alignment.query_name,
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=font_size_pt,
                    zorder=3,
                )
    ax.set_xlim((0, target_len))
    ax.set_ylim((0, 1))
    ax.set_xlabel("Position")
    # ax.legend()
    return ax
