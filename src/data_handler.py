import io
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import logging

logger = logging.getLogger(__name__)

class FastaProcessor:
    def __init__(self, uploaded_file):
        self.uploaded_file = uploaded_file
        self.sequences = self.parse_fasta()

    def parse_fasta(self):
        """Lê o arquivo FASTA e retorna uma lista de tuplas (header, sequência)."""
        sequences = []
        try:
            # Lê o conteúdo do arquivo carregado
            content = self.uploaded_file.getvalue().decode("utf-8")
            fasta_io = io.StringIO(content)
            for record in SeqIO.parse(fasta_io, "fasta"):
                sequences.append((record.id, record.seq))
            logger.info(f"Total de sequências lidas: {len(sequences)}")
        except Exception as e:
            logger.error(f"Erro ao processar o arquivo FASTA: {e}")
            raise
        return sequences

    def generate_peptides(self, sequence, lengths):
        """Gera peptídeos de um determinado comprimento a partir de uma sequência."""
        peptides = []
        for length in lengths:
            for i in range(len(sequence) - length + 1):
                peptide = sequence[i:i+length]
                peptides.append(peptide)
        return peptides

    def calculate_physicochemical(self, sequence):
        """Calcula propriedades físico-químicas de uma sequência proteica."""
        try:
            analysis = ProteinAnalysis(str(sequence))
            iep = analysis.isoelectric_point()
            aromaticity = analysis.aromaticity()
            instability = analysis.instability_index()
            gravy = analysis.gravy()  # Hidrofobicidade média (escala de Kyte-Doolittle)
            mw = analysis.molecular_weight()
            return {
                'isoelectric_point': iep,
                'aromaticity': aromaticity,
                'instability_index': instability,
                'gravy': gravy,
                'molecular_weight': mw
            }
        except Exception as e:
            logger.error(f"Erro ao calcular propriedades físico-químicas: {e}")
            return None

    def get_fasta_data(self):
        """Retorna os dados do FASTA em formato estruturado."""
        data = []
        for header, seq in self.sequences:
            physchem = self.calculate_physicochemical(seq)
            data.append({
                'id': header,
                'sequence': str(seq),
                'length': len(seq),
                'isoelectric_point': physchem['isoelectric_point'] if physchem else None,
                'aromaticity': physchem['aromaticity'] if physchem else None,
                'instability_index': physchem['instability_index'] if physchem else None,
                'gravy': physchem['gravy'] if physchem else None,
                'molecular_weight': physchem['molecular_weight'] if physchem else None
            })
        return data
