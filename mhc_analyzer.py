# -*- coding: utf-8 -*-
"""
Módulo de análise MHC usando MHCflurry.
Gerencia download automático de modelos e predições MHC-I/II.
"""

import logging
import pandas as pd
import numpy as np
from typing import List, Optional, Dict
from sklearn.preprocessing import MinMaxScaler

# Configuração de logging
logger = logging.getLogger(__name__)

# Tenta importar MHCflurry
try:
    from mhcflurry import Class1PresentationPredictor, Class2PresentationPredictor
    from mhcflurry.downloads import get_path
    MHCFLURRY_AVAILABLE = True
except ImportError:
    MHCFLURRY_AVAILABLE = False
    logger.warning("MHCflurry não disponível. Instale com: pip install mhcflurry")


class MHCAnalyzer:
    """Classe para análise de predições MHC usando MHCflurry."""
    
    def __init__(self):
        """Inicializa analisador MHC."""
        self.class1_predictor: Optional[Class1PresentationPredictor] = None
        self.class2_predictor: Optional[Class2PresentationPredictor] = None
        self._models_loaded = {'class1': False, 'class2': False}
        logger.info("MHCAnalyzer inicializado")
    
    def _ensure_class1_models(self) -> bool:
        """
        Garante que modelos MHC-I estão disponíveis, baixando se necessário.
        
        Returns:
            True se modelos estão disponíveis, False caso contrário
        """
        if not MHCFLURRY_AVAILABLE:
            logger.error("MHCflurry não está instalado")
            return False
        
        if self._models_loaded['class1']:
            return True
        
        try:
            logger.info("Carregando modelos MHC-I...")
            self.class1_predictor = Class1PresentationPredictor.load()
            self._models_loaded['class1'] = True
            logger.info("Modelos MHC-I carregados com sucesso")
            return True
        except Exception as e:
            logger.error(f"Erro ao carregar modelos MHC-I: {e}")
            logger.info("Tentando baixar modelos MHC-I...")
            try:
                # Tenta baixar modelos automaticamente
                import subprocess
                result = subprocess.run(
                    ['mhcflurry-downloads', 'fetch', 'models_class1_presentation'],
                    capture_output=True,
                    text=True,
                    timeout=300
                )
                if result.returncode == 0:
                    logger.info("Modelos MHC-I baixados com sucesso")
                    self.class1_predictor = Class1PresentationPredictor.load()
                    self._models_loaded['class1'] = True
                    return True
                else:
                    logger.error(f"Falha ao baixar modelos: {result.stderr}")
                    return False
            except Exception as download_error:
                logger.error(f"Erro ao baixar modelos MHC-I: {download_error}")
                return False
    
    def _ensure_class2_models(self) -> bool:
        """
        Garante que modelos MHC-II estão disponíveis, baixando se necessário.
        
        Returns:
            True se modelos estão disponíveis, False caso contrário
        """
        if not MHCFLURRY_AVAILABLE:
            return False
        
        if self._models_loaded['class2']:
            return True
        
        try:
            logger.info("Carregando modelos MHC-II...")
            self.class2_predictor = Class2PresentationPredictor.load()
            self._models_loaded['class2'] = True
            logger.info("Modelos MHC-II carregados com sucesso")
            return True
        except Exception as e:
            logger.warning(f"Modelos MHC-II não disponíveis: {e}")
            return False
    
    def predict_class1(self, df: pd.DataFrame, alleles: List[str]) -> pd.DataFrame:
        """
        Executa predições MHC-I usando MHCflurry.
        
        Args:
            df: DataFrame com coluna 'peptide'
            alleles: Lista de alelos HLA-I (ex: ['HLA-A*02:01'])
            
        Returns:
            DataFrame enriquecido com predições MHC-I
        """
        if not self._ensure_class1_models():
            logger.warning("Modelos MHC-I não disponíveis. Retornando DataFrame original.")
            return df
        
        if 'peptide' not in df.columns:
            raise ValueError("DataFrame deve conter coluna 'peptide'")
        
        try:
            logger.info(f"Executando predições MHC-I para {len(df)} peptídeos com alelos {alleles}")
            
            preds = []
            for allele in alleles:
                results = self.class1_predictor.predict(
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
            
            # Renomeia colunas
            if 'affinity' in df_result.columns:
                df_result['affinity_nm'] = df_result['affinity']
            if 'presentation_score' in df_result.columns:
                df_result['mhc_score'] = df_result['presentation_score']
            if 'allele' in df_result.columns:
                df_result['allele_class1'] = df_result['allele']
            
            # Percentile rank
            if 'presentation_percentile' in best_preds.columns:
                df_result['percentile_rank'] = df_result['peptide'].map(
                    best_preds.set_index('peptide')['presentation_percentile']
                )
            else:
                df_result['percentile_rank'] = df_result['mhc_score'].rank(pct=True, ascending=False)
            
            logger.info("Predições MHC-I concluídas com sucesso")
            return df_result
            
        except Exception as e:
            logger.error(f"Erro nas predições MHC-I: {e}")
            return df
    
    def predict_class2(self, df: pd.DataFrame, alleles: List[str]) -> pd.DataFrame:
        """
        Executa predições MHC-II usando MHCflurry.
        
        Args:
            df: DataFrame com coluna 'peptide'
            alleles: Lista de alelos HLA-II (ex: ['HLA-DRB1*01:01'])
            
        Returns:
            DataFrame enriquecido com predições MHC-II
        """
        if not self._ensure_class2_models():
            logger.warning("Modelos MHC-II não disponíveis. Retornando DataFrame original.")
            return df
        
        if 'peptide' not in df.columns:
            raise ValueError("DataFrame deve conter coluna 'peptide'")
        
        try:
            logger.info(f"Executando predições MHC-II para {len(df)} peptídeos com alelos {alleles}")
            
            preds = []
            for allele in alleles:
                results = self.class2_predictor.predict(
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
            
            # Merge com sufixos para evitar conflitos
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
                df_result['allele_class2'] = df_result['allele']
            
            logger.info("Predições MHC-II concluídas com sucesso")
            return df_result
            
        except Exception as e:
            logger.error(f"Erro nas predições MHC-II: {e}")
            return df
    
    def calculate_final_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calcula scores finais normalizados para ranking.
        
        Args:
            df: DataFrame com predições
            
        Returns:
            DataFrame com scores finais calculados
        """
        df = df.copy()
        scaler = MinMaxScaler()
        
        # Normaliza score MHC-I
        if 'mhc_score' in df.columns:
            try:
                df['norm_mhc'] = scaler.fit_transform(df[['mhc_score']].fillna(0))
            except Exception as e:
                logger.warning(f"Erro ao normalizar MHC-I: {e}")
                df['norm_mhc'] = 0.5
        else:
            df['norm_mhc'] = 0.0
        
        # Normaliza score MHC-II (se disponível)
        if 'mhcii_score' in df.columns:
            try:
                df['norm_mhcii'] = scaler.fit_transform(df[['mhcii_score']].fillna(0))
            except Exception as e:
                logger.warning(f"Erro ao normalizar MHC-II: {e}")
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
        
        logger.info("Scores finais calculados")
        return df
