import csv
from collections.abc import Iterator
from itertools import chain
from pathlib import Path

from pyryotype.paf_plotting import PAFProtocol, _collapse_multiple_mappings


def get_metadata_csv(alignments: Iterator[PAFProtocol], filename: str | Path, mapq_threshold: int = 0) -> None:
    """
    Write out a CSV for use with bandage, which parses PAF files and
    puts out metadata for labelling various contigs in bandage, such as
    which target on the reference the contig aligns to

    :param alignments: An iterator which yields objects that implement the PAFProtocol methods
    :param filename: The filename to write the metadata to. Raises FileExistsError if the file already exists.
    :param mapq_threshold: Filter out alignments under this threshold. Defaults 0.
    """
    headers = ["query_name", "target_name", "target_start", "target_end", "mapping_quality"]
    with open(filename, "x", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([*headers, "blast_identity"])
        for alignment in _collapse_multiple_mappings(alignments):
            if not alignment.mapping_quality < mapq_threshold:
                writer.writerow(chain((getattr(alignment, h, "*") for h in headers), (alignment.blast_identity(),)))
