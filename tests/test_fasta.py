import unittest
from src.data_handler import FastaProcessor
import io

class TestFastaProcessor(unittest.TestCase):
    def test_parse_fasta(self):
        fasta_content = """>seq1
MSTNGER
>seq2
MKLLSI"""
        uploaded_file = io.BytesIO(fasta_content.encode())
        processor = FastaProcessor(uploaded_file)
        sequences = processor.sequences
        self.assertEqual(len(sequences), 2)
        self.assertEqual(sequences[0][0], "seq1")
        self.assertEqual(str(sequences[0][1]), "MSTNGER")

if __name__ == '__main__':
    unittest.main()
