# -*- coding: utf-8 -*-
"""
Módulo de processamento de dados FASTA.
Classe profissional para leitura e validação de sequências.
"""

import re
import io
import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Union
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Configuração de logging
logger = logging.getLogger(__name__)


class FastaProcessor:
    """Classe para processamento de arquivos FASTA e validação de peptídeos."""
    
    # Escala de hidrofobicidade Kyte-Doolittle (1982)
    KD_SCALE = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5,
        'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
        'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
    
    def __init__(self, min_length: int = 8, max_length: int = 14):
        """
        Inicializa processador de FASTA.
        
        Args:
            min_length: Tamanho mínimo do peptídeo (padrão: 8 para MHC-I)
            max_length: Tamanho máximo do peptídeo (padrão: 14 para MHC-I)
        """
        self.min_length = min_length
        self.max_length = max_length
        self.aa_regex = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")
        logger.info(f"FastaProcessor inicializado (tamanho: {min_length}-{max_length})")
    
    def read_fasta_from_bytes(self, file_bytes: bytes) -> pd.DataFrame:
        """
        Lê arquivo FASTA a partir de bytes.
        
        Args:
            file_bytes: Conteúdo do arquivo em bytes
            
        Returns:
            DataFrame com coluna 'peptide'
            
        Raises:
            ValueError: Se nenhuma sequência válida for encontrada
        """
        try:
            file_string = file_bytes.decode('utf-8')
            file_io = io.StringIO(file_string)
            
            peptides = []
            for record in SeqIO.parse(file_io, "fasta"):
                seq = str(record.seq).upper()
                seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq)
                if seq:
                    peptides.append(seq)
            
            if not peptides:
                raise ValueError("Nenhuma sequência válida encontrada no arquivo FASTA")
            
            logger.info(f"Lidos {len(peptides)} sequências do arquivo FASTA")
            return pd.DataFrame({'peptide': peptides})
            
        except Exception as e:
            logger.error(f"Erro ao ler FASTA: {e}")
            raise ValueError(f"Erro ao ler arquivo FASTA: {e}")
    
    def read_tsv_csv_from_bytes(self, file_bytes: bytes, filename: str, 
                                sep: Optional[str] = None) -> pd.DataFrame:
        """
        Lê arquivo TSV/CSV a partir de bytes.
        
        Args:
            file_bytes: Conteúdo do arquivo em bytes
            filename: Nome do arquivo (para detectar formato)
            sep: Separador (None para detecção automática)
            
        Returns:
            DataFrame com coluna 'peptide'
        """
        try:
            if filename.lower().endswith('.tsv') or (sep is None and '\t' in file_bytes.decode('utf-8', errors='ignore')[:1000]):
                df = pd.read_csv(io.BytesIO(file_bytes), sep='\t', header=None)
            else:
                df = pd.read_csv(io.BytesIO(file_bytes), sep=sep if sep else ',', header=None)
            
            peptides_flat = df.values.flatten()
            df_clean = pd.DataFrame(peptides_flat, columns=['peptide'])
            df_clean = df_clean.dropna()
            df_clean['peptide'] = df_clean['peptide'].astype(str).str.strip()
            df_clean = df_clean[df_clean['peptide'] != '']
            
            logger.info(f"Lidos {len(df_clean)} peptídeos do arquivo {filename}")
            return df_clean
            
        except Exception as e:
            logger.error(f"Erro ao ler TSV/CSV: {e}")
            raise ValueError(f"Erro ao ler arquivo TSV/CSV: {e}")
    
    def read_excel_from_bytes(self, file_bytes: bytes, 
                              sheet_name: Optional[Union[str, int]] = 0) -> pd.DataFrame:
        """
        Lê arquivo Excel a partir de bytes.
        
        Args:
            file_bytes: Conteúdo do arquivo em bytes
            sheet_name: Nome ou índice da planilha
            
        Returns:
            DataFrame com coluna 'peptide'
        """
        try:
            df = pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet_name)
            
            if 'peptide' in df.columns:
                df_clean = df[['peptide']].copy()
            else:
                for col in df.columns:
                    if df[col].astype(str).str.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$').any():
                        df_clean = pd.DataFrame({'peptide': df[col]})
                        break
                else:
                    df_clean = pd.DataFrame({'peptide': df.iloc[:, 0]})
            
            df_clean = df_clean.dropna()
            df_clean['peptide'] = df_clean['peptide'].astype(str).str.strip()
            df_clean = df_clean[df_clean['peptide'] != '']
            
            logger.info(f"Lidos {len(df_clean)} peptídeos do arquivo Excel")
            return df_clean
            
        except Exception as e:
            logger.error(f"Erro ao ler Excel: {e}")
            raise ValueError(f"Erro ao ler arquivo Excel: {e}")
    
    def load_peptides_from_bytes(self, file_bytes: bytes, filename: str) -> pd.DataFrame:
        """
        Carrega peptídeos a partir de bytes (para uso com Streamlit file_uploader).
        
        Args:
            file_bytes: Conteúdo do arquivo em bytes
            filename: Nome do arquivo (para detectar formato)
            
        Returns:
            DataFrame com coluna 'peptide'
        """
        filename_lower = filename.lower()
        
        if filename_lower.endswith(('.fasta', '.fa', '.faa')):
            return self.read_fasta_from_bytes(file_bytes)
        elif filename_lower.endswith(('.tsv', '.txt')):
            return self.read_tsv_csv_from_bytes(file_bytes, filename, sep='\t')
        elif filename_lower.endswith('.csv'):
            return self.read_tsv_csv_from_bytes(file_bytes, filename, sep=',')
        elif filename_lower.endswith(('.xlsx', '.xls')):
            return self.read_excel_from_bytes(file_bytes)
        else:
            logger.warning(f"Formato desconhecido {filename}, tentando como TSV")
            return self.read_tsv_csv_from_bytes(file_bytes, filename, sep='\t')
    
    def validate_peptides(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Valida e filtra peptídeos baseado em critérios biológicos.
        
        Args:
            df: DataFrame com coluna 'peptide'
            
        Returns:
            DataFrame filtrado com apenas peptídeos válidos
        """
        if 'peptide' not in df.columns:
            raise ValueError("DataFrame deve conter coluna 'peptide'")
        
        initial_count = len(df)
        df = df.copy()
        
        # Padroniza para maiúsculas
        df['peptide'] = df['peptide'].astype(str).str.upper().str.strip()
        
        # Valida sequências
        df["valid"] = df['peptide'].apply(lambda x: bool(self.aa_regex.match(x)))
        
        # Remove inválidos e duplicatas
        df_valid = df[df["valid"]].drop(columns=["valid"]).drop_duplicates(subset=['peptide']).copy()
        
        # Filtro de tamanho
        df_valid = df_valid[df_valid['peptide'].str.len().between(self.min_length, self.max_length)]
        
        final_count = len(df_valid)
        logger.info(f"Validação: {initial_count} → {final_count} peptídeos válidos")
        
        return df_valid.reset_index(drop=True)
    
    @staticmethod
    def calculate_kyte_doolittle_hydrophobicity(peptide: str) -> float:
        """
        Calcula hidrofobicidade média usando escala de Kyte-Doolittle.
        
        Args:
            peptide: Sequência do peptídeo
            
        Returns:
            Hidrofobicidade média (valores positivos = hidrofóbico)
        """
        if not peptide:
            return 0.0
        
        try:
            hydrophobicity_values = [FastaProcessor.KD_SCALE.get(aa, 0.0) for aa in peptide.upper()]
            return round(sum(hydrophobicity_values) / len(hydrophobicity_values), 3)
        except Exception as e:
            logger.warning(f"Erro ao calcular hidrofobicidade para {peptide}: {e}")
            return 0.0
    
    def calculate_physchem_properties(self, peptide: str) -> Dict[str, float]:
        """
        Calcula propriedades físico-químicas de um peptídeo.
        
        Args:
            peptide: Sequência do peptídeo
            
        Returns:
            Dicionário com propriedades calculadas
        """
        try:
            pa = ProteinAnalysis(peptide)
            
            properties = {
                "length": len(peptide),
                "mw": round(pa.molecular_weight(), 2),
                "pI": round(pa.isoelectric_point(), 2),
                "gravy": round(pa.gravy(), 3),
                "kd_hydrophobicity": self.calculate_kyte_doolittle_hydrophobicity(peptide),
            }
            
            # Índice de instabilidade
            try:
                properties["instability_index"] = round(pa.instability_index(), 2)
            except:
                properties["instability_index"] = 0.0
            
            # Índice alifático
            try:
                aa_percent = pa.get_amino_acids_percent()
                properties["aliphatic_index"] = round(
                    (aa_percent.get('A', 0) + 2.9 * aa_percent.get('V', 0) + 
                     3.9 * (aa_percent.get('I', 0) + aa_percent.get('L', 0))) * 100, 2
                )
            except:
                properties["aliphatic_index"] = 0.0
            
            return properties
            
        except Exception as e:
            logger.warning(f"Erro ao calcular propriedades para {peptide}: {e}")
            return {
                "length": len(peptide),
                "mw": 0.0,
                "pI": 0.0,
                "gravy": 0.0,
                "kd_hydrophobicity": 0.0,
                "instability_index": 0.0,
                "aliphatic_index": 0.0
            }
    
    def add_physchem_properties(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Adiciona propriedades físico-químicas ao DataFrame.
        
        Args:
            df: DataFrame com coluna 'peptide'
            
        Returns:
            DataFrame com colunas adicionais de propriedades
        """
        if 'peptide' not in df.columns:
            raise ValueError("DataFrame deve conter coluna 'peptide'")
        
        logger.info(f"Calculando propriedades físico-químicas para {len(df)} peptídeos...")
        
        props = df['peptide'].apply(self.calculate_physchem_properties).apply(pd.Series)
        
        # Remove colunas duplicadas antes de concatenar
        cols_to_remove = [c for c in props.columns if c in df.columns]
        if cols_to_remove:
            df = df.drop(columns=cols_to_remove)
        
        df_result = pd.concat([df.reset_index(drop=True), props.reset_index(drop=True)], axis=1)
        
        logger.info("Propriedades físico-químicas calculadas com sucesso")
        return df_result
