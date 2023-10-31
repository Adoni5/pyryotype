import unittest

from pyryotype.paf_plotting import PAFProtocol


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
