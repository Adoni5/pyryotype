import csv
from collections.abc import Iterator
from itertools import chain
from pathlib import Path

from pyryotype.paf_plotting import PAFProtocol


def get_metadata_csv(alignments: Iterator[PAFProtocol], filename: str | Path) -> None:
    """
    Write out a CSV for use with bandage, which parses PAF files and
    puts out metadata for labelling various contigs in bandage, such as
    which target on the reference the contig aligns to

    :param alignments: An iterator which yields objects that implement the PAFProtocol methods
    :param filename: The filename to write the metadata to. Will
    """
    headers = ["query_name", "target_name", "target_start", "target_end", "mapping_quality"]
    with open(filename, "x", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([*headers, "blast_identity"])
        for alignment in alignments:
            writer.writerow(chain((getattr(alignment, h, "*") for h in headers), (alignment.blast_identity(),)))
