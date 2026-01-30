# -*- coding: utf-8 -*-
"""
Módulo de análise de conservação de sequências.
Calcula conservação usando múltiplos métodos.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from collections import Counter
from Bio import Align
from Bio.Align import substitution_matrices


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
    
    # Conta frequência de cada aminoácido
    counter = Counter(sequences)
    total = len(sequences)
    
    # Calcula entropia
    entropy = 0.0
    for count in counter.values():
        p = count / total
        if p > 0:
            entropy -= p * np.log2(p)
    
    return entropy


def calculate_conservation_score(peptide: str, reference_sequences: List[str]) -> Dict[str, float]:
    """
    Calcula score de conservação de um peptídeo em relação a sequências de referência.
    
    Args:
        peptide: Sequência do peptídeo
        reference_sequences: Lista de sequências de referência (proteoma, etc.)
        
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
    
    # Alinha peptídeo com cada sequência de referência
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    
    identities = []
    matches = 0
    
    for ref_seq in reference_sequences:
        try:
            alignments = aligner.align(peptide, ref_seq)
            if alignments:
                alignment = alignments[0]
                # Calcula identidade
                identity = alignment.score / (len(peptide) * 2)  # Normalizado
                identities.append(identity)
                
                if identity > 0.8:  # 80% de identidade
                    matches += 1
        except:
            continue
    
    if not identities:
        return {
            'conservation_score': 0.0,
            'max_identity': 0.0,
            'mean_identity': 0.0,
            'matches': 0
        }
    
    # Scores
    max_identity = max(identities)
    mean_identity = np.mean(identities)
    
    # Score de conservação (0-1): combina max e mean
    conservation_score = (max_identity * 0.6 + mean_identity * 0.4)
    
    return {
        'conservation_score': round(conservation_score, 4),
        'max_identity': round(max_identity, 4),
        'mean_identity': round(mean_identity, 4),
        'matches': matches
    }


def calculate_positional_conservation(peptides: List[str]) -> pd.DataFrame:
    """
    Calcula conservação posicional entre múltiplos peptídeos.
    
    Args:
        peptides: Lista de sequências de peptídeos
        
    Returns:
        DataFrame com conservação por posição
    """
    if not peptides:
        return pd.DataFrame()
    
    # Encontra tamanho máximo
    max_len = max(len(p) for p in peptides)
    
    # Alinha sequências (padding com gaps se necessário)
    aligned = []
    for p in peptides:
        aligned.append(p.ljust(max_len, '-'))
    
    # Calcula conservação por posição
    positions = []
    for pos in range(max_len):
        column = [seq[pos] for seq in aligned if seq[pos] != '-']
        if column:
            entropy = calculate_shannon_entropy(column)
            # Converte entropia para score de conservação (0-1)
            # Entropia máxima para 20 AA = log2(20) ≈ 4.32
            conservation = 1 - (entropy / 4.32) if entropy > 0 else 1.0
            positions.append({
                'position': pos + 1,
                'entropy': round(entropy, 4),
                'conservation': round(conservation, 4),
                'most_common_aa': Counter(column).most_common(1)[0][0] if column else '-'
            })
    
    return pd.DataFrame(positions)


def add_conservation_to_dataframe(df: pd.DataFrame, reference_sequences: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Adiciona scores de conservação ao DataFrame.
    
    Args:
        df: DataFrame com coluna 'peptide'
        reference_sequences: Sequências de referência (opcional)
        
    Returns:
        DataFrame com colunas de conservação adicionadas
    """
    if 'peptide' not in df.columns:
        return df
    
    df_result = df.copy()
    
    # Conservação posicional entre peptídeos
    peptides = df_result['peptide'].tolist()
    pos_conservation = calculate_positional_conservation(peptides)
    
    # Score médio de conservação posicional
    if not pos_conservation.empty:
        mean_conservation = pos_conservation['conservation'].mean()
        df_result['positional_conservation'] = mean_conservation
    else:
        df_result['positional_conservation'] = 0.0
    
    # Conservação em relação a sequências de referência (se fornecidas)
    if reference_sequences:
        conservation_scores = []
        for peptide in peptides:
            scores = calculate_conservation_score(peptide, reference_sequences)
            conservation_scores.append(scores['conservation_score'])
        df_result['reference_conservation'] = conservation_scores
    else:
        df_result['reference_conservation'] = 0.0
    
    return df_result
