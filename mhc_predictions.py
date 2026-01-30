# -*- coding: utf-8 -*-
"""
Módulo de predições MHC (MHC-I e MHC-II).
Integra MHCflurry, NetMHCpan e outras ferramentas.
"""

import pandas as pd
import numpy as np
from typing import List, Optional
import warnings
warnings.filterwarnings('ignore')

# Tenta importar MHCflurry
try:
    from mhcflurry import Class1PresentationPredictor, Class2PresentationPredictor
    MHCFLURRY_AVAILABLE = True
except ImportError:
    MHCFLURRY_AVAILABLE = False
    print("⚠️ MHCflurry não disponível. Instale com: pip install mhcflurry")


def run_mhcflurry_class1_predictions(df: pd.DataFrame, alleles: List[str]) -> pd.DataFrame:
    """
    Executa predições MHC-I usando MHCflurry.
    
    Args:
        df: DataFrame com coluna 'peptide'
        alleles: Lista de alelos HLA-I (ex: ['HLA-A*02:01'])
        
    Returns:
        DataFrame enriquecido com predições MHC-I
    """
    if not MHCFLURRY_AVAILABLE:
        print("⚠️ MHCflurry não disponível. Pulando predições MHC-I.")
        return df
    
    try:
        predictor = Class1PresentationPredictor.load()
        
        preds = []
        for allele in alleles:
            results = predictor.predict(
                peptides=df['peptide'].values,
                alleles=[allele],
                include_affinity_percentile=True,
                verbose=0
            )
            results['allele'] = allele
            preds.append(results)
        
        # Combina e pega melhor alelo
        full_preds = pd.concat(preds)
        best_preds = full_preds.loc[full_preds.groupby('peptide')['presentation_score'].idxmax()]
        
        # Merge
        df_result = pd.merge(
            df,
            best_preds[['peptide', 'affinity', 'presentation_score', 'allele']],
            on='peptide',
            how='left'
        )
        
        # Renomeia
        if 'affinity' in df_result.columns:
            df_result['affinity_nm'] = df_result['affinity']
        if 'presentation_score' in df_result.columns:
            df_result['mhc_score'] = df_result['presentation_score']
        
        # Percentile rank
        if 'presentation_percentile' in best_preds.columns:
            df_result['percentile_rank'] = df_result['peptide'].map(
                best_preds.set_index('peptide')['presentation_percentile']
            )
        else:
            df_result['percentile_rank'] = df_result['mhc_score'].rank(pct=True, ascending=False)
        
        return df_result
        
    except Exception as e:
        print(f"⚠️ Erro nas predições MHC-I: {e}")
        return df


def run_mhcflurry_class2_predictions(df: pd.DataFrame, alleles: List[str]) -> pd.DataFrame:
    """
    Executa predições MHC-II usando MHCflurry.
    
    Args:
        df: DataFrame com coluna 'peptide'
        alleles: Lista de alelos HLA-II (ex: ['HLA-DRB1*01:01'])
        
    Returns:
        DataFrame enriquecido com predições MHC-II
    """
    if not MHCFLURRY_AVAILABLE:
        print("⚠️ MHCflurry não disponível. Pulando predições MHC-II.")
        return df
    
    try:
        predictor = Class2PresentationPredictor.load()
        
        preds = []
        for allele in alleles:
            results = predictor.predict(
                peptides=df['peptide'].values,
                alleles=[allele],
                include_affinity_percentile=True,
                verbose=0
            )
            results['allele'] = allele
            preds.append(results)
        
        # Combina e pega melhor alelo
        full_preds = pd.concat(preds)
        best_preds = full_preds.loc[full_preds.groupby('peptide')['presentation_score'].idxmax()]
        
        # Merge
        df_result = pd.merge(
            df,
            best_preds[['peptide', 'affinity', 'presentation_score', 'allele']],
            on='peptide',
            how='left',
            suffixes=('', '_mhcii')
        )
        
        # Renomeia para MHC-II
        if 'affinity' in df_result.columns:
            df_result['affinity_mhcii_nm'] = df_result['affinity']
        if 'presentation_score' in df_result.columns:
            df_result['mhcii_score'] = df_result['presentation_score']
        if 'allele' in df_result.columns:
            df_result['allele_mhcii'] = df_result['allele']
        
        return df_result
        
    except Exception as e:
        print(f"⚠️ Erro nas predições MHC-II: {e}")
        return df


def calculate_final_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calcula scores finais normalizados para ranking.
    
    Args:
        df: DataFrame com predições
        
    Returns:
        DataFrame com scores finais calculados
    """
    df = df.copy()
    
    from sklearn.preprocessing import MinMaxScaler
    scaler = MinMaxScaler()
    
    # Normaliza score MHC-I
    if 'mhc_score' in df.columns:
        try:
            df['norm_mhc'] = scaler.fit_transform(df[['mhc_score']].fillna(0))
        except:
            df['norm_mhc'] = 0.5
    else:
        df['norm_mhc'] = 0.0
    
    # Normaliza score MHC-II (se disponível)
    if 'mhcii_score' in df.columns:
        try:
            df['norm_mhcii'] = scaler.fit_transform(df[['mhcii_score']].fillna(0))
        except:
            df['norm_mhcii'] = 0.5
    else:
        df['norm_mhcii'] = 0.0
    
    # Score final combinado (prioriza MHC-I)
    if 'norm_mhc' in df.columns and 'norm_mhcii' in df.columns:
        df['final_rank_score'] = 0.7 * df['norm_mhc'] + 0.3 * df['norm_mhcii']
    elif 'norm_mhc' in df.columns:
        df['final_rank_score'] = df['norm_mhc']
    elif 'norm_mhcii' in df.columns:
        df['final_rank_score'] = df['norm_mhcii']
    else:
        df['final_rank_score'] = 0.0
    
    # Ordena por score final
    df = df.sort_values(by='final_rank_score', ascending=False).reset_index(drop=True)
    
    return df
