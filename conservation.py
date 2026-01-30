# -*- coding: utf-8 -*-
"""
Módulo de análise de conservação de sequências.
Calcula conservação usando múltiplos métodos.
"""

import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from collections import Counter

# Configuração de logging
logger = logging.getLogger(__name__)

# Tenta importar Bio.Align
try:
    from Bio import Align
    from Bio.Align import substitution_matrices
    BIO_ALIGN_AVAILABLE = True
except ImportError:
    BIO_ALIGN_AVAILABLE = False
    logger.warning("Bio.Align não disponível. Análise de conservação limitada.")


class ConservationAnalyzer:
    """Classe para análise de conservação de sequências."""
    
    def __init__(self):
        """Inicializa analisador de conservação."""
        logger.info("ConservationAnalyzer inicializado")
        if BIO_ALIGN_AVAILABLE:
            self.aligner = Align.PairwiseAligner()
            self.aligner.mode = 'local'
            try:
                self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            except:
                logger.warning("Não foi possível carregar matriz BLOSUM62")
    
    @staticmethod
    def calculate_shannon_entropy(sequences: List[str]) -> float:
        """
        Calcula entropia de Shannon para uma posição alinhada.
        
        Args:
            sequences: Lista de sequências alinhadas
            
        Returns:
            Entropia de Shannon
        """
        if not sequences:
            return 0.0
        
        try:
            counter = Counter(sequences)
            total = len(sequences)
            
            entropy = 0.0
            for count in counter.values():
                p = count / total
                if p > 0:
                    entropy -= p * np.log2(p)
            
            return entropy
        except Exception as e:
            logger.warning(f"Erro ao calcular entropia: {e}")
            return 0.0
    
    def calculate_conservation_score(self, peptide: str, reference_sequences: List[str]) -> Dict[str, float]:
        """
        Calcula score de conservação de um peptídeo em relação a sequências de referência.
        
        Args:
            peptide: Sequência do peptídeo
            reference_sequences: Lista de sequências de referência
            
        Returns:
            Dicionário com scores de conservação
        """
        if not reference_sequences:
            return {
                'conservation_score': 0.0,
                'max_identity': 0.0,
                'mean_identity': 0.0,
                'matches': 0
            }
        
        if not BIO_ALIGN_AVAILABLE:
            logger.warning("Bio.Align não disponível. Retornando scores padrão.")
            return {
                'conservation_score': 0.0,
                'max_identity': 0.0,
                'mean_identity': 0.0,
                'matches': 0
            }
        
        identities = []
        matches = 0
        
        try:
            for ref_seq in reference_sequences:
                try:
                    alignments = self.aligner.align(peptide, ref_seq)
                    if alignments:
                        alignment = alignments[0]
                        identity = alignment.score / (len(peptide) * 2)
                        identities.append(identity)
                        
                        if identity > 0.8:
                            matches += 1
                except Exception as e:
                    logger.debug(f"Erro ao alinhar {peptide} com {ref_seq}: {e}")
                    continue
        except Exception as e:
            logger.warning(f"Erro na análise de conservação: {e}")
            return {
                'conservation_score': 0.0,
                'max_identity': 0.0,
                'mean_identity': 0.0,
                'matches': 0
            }
        
        if not identities:
            return {
                'conservation_score': 0.0,
                'max_identity': 0.0,
                'mean_identity': 0.0,
                'matches': 0
            }
        
        max_identity = max(identities)
        mean_identity = np.mean(identities)
        conservation_score = (max_identity * 0.6 + mean_identity * 0.4)
        
        return {
            'conservation_score': round(conservation_score, 4),
            'max_identity': round(max_identity, 4),
            'mean_identity': round(mean_identity, 4),
            'matches': matches
        }
    
    def calculate_positional_conservation(self, peptides: List[str]) -> pd.DataFrame:
        """
        Calcula conservação posicional entre múltiplos peptídeos.
        
        Args:
            peptides: Lista de sequências de peptídeos
            
        Returns:
            DataFrame com conservação por posição
        """
        if not peptides:
            return pd.DataFrame()
        
        try:
            max_len = max(len(p) for p in peptides)
            aligned = [p.ljust(max_len, '-') for p in peptides]
            
            positions = []
            for pos in range(max_len):
                column = [seq[pos] for seq in aligned if seq[pos] != '-']
                if column:
                    entropy = self.calculate_shannon_entropy(column)
                    conservation = 1 - (entropy / 4.32) if entropy > 0 else 1.0
                    positions.append({
                        'position': pos + 1,
                        'entropy': round(entropy, 4),
                        'conservation': round(conservation, 4),
                        'most_common_aa': Counter(column).most_common(1)[0][0] if column else '-'
                    })
            
            return pd.DataFrame(positions)
        except Exception as e:
            logger.error(f"Erro ao calcular conservação posicional: {e}")
            return pd.DataFrame()
    
    def add_conservation_to_dataframe(self, df: pd.DataFrame, 
                                     reference_sequences: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Adiciona scores de conservação ao DataFrame.
        
        Args:
            df: DataFrame com coluna 'peptide'
            reference_sequences: Sequências de referência (opcional)
            
        Returns:
            DataFrame com colunas de conservação adicionadas
        """
        if 'peptide' not in df.columns:
            raise ValueError("DataFrame deve conter coluna 'peptide'")
        
        df_result = df.copy()
        peptides = df_result['peptide'].tolist()
        
        logger.info(f"Calculando conservação para {len(peptides)} peptídeos...")
        
        # Conservação posicional entre peptídeos
        try:
            pos_conservation = self.calculate_positional_conservation(peptides)
            if not pos_conservation.empty:
                mean_conservation = pos_conservation['conservation'].mean()
                df_result['positional_conservation'] = mean_conservation
            else:
                df_result['positional_conservation'] = 0.0
        except Exception as e:
            logger.warning(f"Erro ao calcular conservação posicional: {e}")
            df_result['positional_conservation'] = 0.0
        
        # Conservação em relação a sequências de referência (se fornecidas)
        if reference_sequences:
            try:
                conservation_scores = []
                for peptide in peptides:
                    scores = self.calculate_conservation_score(peptide, reference_sequences)
                    conservation_scores.append(scores['conservation_score'])
                df_result['reference_conservation'] = conservation_scores
            except Exception as e:
                logger.warning(f"Erro ao calcular conservação de referência: {e}")
                df_result['reference_conservation'] = 0.0
        else:
            df_result['reference_conservation'] = 0.0
        
        logger.info("Conservação calculada com sucesso")
        return df_result
