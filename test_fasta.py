# -*- coding: utf-8 -*-
"""
Testes básicos para validação de leitura de FASTA.
"""

import unittest
import pandas as pd
import io
import sys
from pathlib import Path

# Adiciona src ao path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_processor import FastaProcessor


class TestFastaProcessor(unittest.TestCase):
    """Testes para FastaProcessor."""
    
    def setUp(self):
        """Configuração inicial."""
        self.processor = FastaProcessor(min_length=8, max_length=14)
    
    def test_read_fasta_from_bytes(self):
        """Testa leitura de FASTA a partir de bytes."""
        # Cria FASTA de exemplo
        fasta_content = """>seq1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
>seq2
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL
"""
        file_bytes = fasta_content.encode('utf-8')
        
        # Testa leitura
        df = self.processor.read_fasta_from_bytes(file_bytes)
        
        # Verifica resultado
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn('peptide', df.columns)
        self.assertGreater(len(df), 0)
    
    def test_validate_peptides(self):
        """Testa validação de peptídeos."""
        # Cria DataFrame de teste
        test_data = pd.DataFrame({
            'peptide': [
                'MKTAYIAK',  # Válido (8 resíduos)
                'MKTAYIAKQRQISFVK',  # Válido (16 resíduos, mas será filtrado por tamanho)
                'MKTAYIAKQR',  # Válido (10 resíduos)
                'MKTAYIA',  # Inválido (7 resíduos)
                'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL',  # Inválido (muito longo)
                'MKTAYIAKX',  # Inválido (contém X)
                'MKTAYIAK',  # Duplicata
            ]
        })
        
        # Valida
        df_valid = self.processor.validate_peptides(test_data)
        
        # Verifica resultado
        self.assertIsInstance(df_valid, pd.DataFrame)
        self.assertIn('peptide', df_valid.columns)
        # Deve ter pelo menos 1 válido (MKTAYIAKQR)
        self.assertGreaterEqual(len(df_valid), 1)
        # Não deve ter duplicatas
        self.assertEqual(len(df_valid), len(df_valid.drop_duplicates()))
    
    def test_calculate_physchem_properties(self):
        """Testa cálculo de propriedades físico-químicas."""
        peptide = "MKTAYIAKQR"
        properties = self.processor.calculate_physchem_properties(peptide)
        
        # Verifica propriedades
        self.assertIn('length', properties)
        self.assertIn('mw', properties)
        self.assertIn('pI', properties)
        self.assertIn('gravy', properties)
        self.assertEqual(properties['length'], len(peptide))


if __name__ == '__main__':
    unittest.main()
