from pathlib import Path

from pyryotype import get_metadata_csv
from readpaf import parse_paf


def test_paf_metadata_csv():
    test_paf = Path(__file__).parent / "static" / "test_extraneous.paf"
    get_metadata_csv(parse_paf(test_paf.open()), Path("test_metadata.csv"))
