import csv
from collections import defaultdict
from collections.abc import Iterator
from itertools import chain
from pathlib import Path

from pyryotype.paf_plotting import PAF, PAFProtocol


def _choose_largest_alignments_block(alignments: Iterator[PAFProtocol]) -> Iterator[PAFProtocol]:
    """

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
    >>> result = list(_choose_largest_alignments_block(iter(alignments)))
    >>> len(result)
    1
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
    >>> result = list(_choose_largest_alignments_block(iter(alignments)))
    >>> len(result)
    2

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
            most_contig = defaultdict(int)
            for (_query_name, contig), als in contig_mapping_proportion.items():
                for al in als:
                    most_contig[contig] += al.target_end - al.target_start
            biggest_al_block_target = max(most_contig, key=most_contig.get)
            # print(f"query name {_query_name}, biggest contig {biggest_al_block_target}")
            # print(contig_mapping_proportion)
            # We have a new query, yield the previous one
            # First we check which contig/strand had the largest amount of alignment to it
            # Then we yield the largest contig/strand
            # Yield supplementary alignments collapsed as well
            for query_name, contig in filter(
                lambda kv: kv[1] == biggest_al_block_target, contig_mapping_proportion.keys()
            ):
                # print(f"hello {alignment}")
                largest_proportion_alignments = sorted(
                    contig_mapping_proportion[(query_name, contig)], key=lambda x: x.target_start
                )
                paf = PAF.from_protocol(largest_proportion_alignments[0])
                largest_proportion_alignments = sorted(
                    contig_mapping_proportion[(query_name, contig)], key=lambda x: x.target_end
                )
                paf.target_end = largest_proportion_alignments[-1].target_end
                yield paf
            curr_id = alignment.query_name
            contig_mapping_proportion.clear()

        contig_mapping_proportion[(alignment.query_name, alignment.target_name)].append(alignment)
    # Yield the last contig alignment
    most_contig = defaultdict(int)
    for (_query_name, contig), als in contig_mapping_proportion.items():
        for al in als:
            most_contig[contig] += al.target_end - al.target_start
    biggest_al_block_target = max(most_contig, key=most_contig.get)
    for query_name, contig in filter(lambda kv: kv[1] == biggest_al_block_target, contig_mapping_proportion.keys()):
        largest_proportion_alignments = sorted(
            contig_mapping_proportion[(query_name, contig)], key=lambda x: x.target_start
        )
        paf = PAF.from_protocol(largest_proportion_alignments[0])
        largest_proportion_alignments = sorted(
            contig_mapping_proportion[(query_name, contig)], key=lambda x: x.target_end
        )
        paf.target_end = largest_proportion_alignments[-1].target_end
        yield paf


def get_metadata_csv(alignments: Iterator[PAFProtocol], filename: str | Path, mapq_threshold: int = 0) -> None:
    """
    Write out a CSV for use with bandage, which parses PAF files and
    puts out metadata for labelling various contigs in bandage, such as
    which target on the reference the contig aligns to

    :param alignments: An iterator which yields objects that implement the PAFProtocol methods
    :param filename: The filename to write the metadata to. Raises FileExistsError if the file already exists.
    :param mapq_threshold: Filter out alignments under this threshold. Defaults 0.
    """
    headers = [
        "query_name",
        "query_length",
        "mapping_quality",
        "target_name",
        "target_start",
        "target_end",
        "target_length",
    ]
    with open(filename, "x", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([*headers, "blast_identity", "fraction_target_length"])
        alignments = sorted(
            filter(
                lambda x: x.mapping_quality >= mapq_threshold,
                alignments,
            ),
            key=lambda x: x.query_name,
        )
        # print(alignments)
        for alignment in _choose_largest_alignments_block(alignments):
            writer.writerow(
                chain(
                    (getattr(alignment, h, "*") for h in headers),
                    (alignment.blast_identity(), alignment.query_length / alignment.target_length),
                )
            )
