# -*- coding: utf-8 -*-
"""
Módulo de gerenciamento de APIs externas (IEDB, UniProt).
Implementa processamento paralelo com tratamento robusto de erros.
"""

import logging
import time
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Optional, Tuple
import pandas as pd
import numpy as np

# Configuração de logging
logger = logging.getLogger(__name__)


class APIManager:
    """Gerenciador de APIs externas com processamento paralelo."""
    
    def __init__(self, max_workers: int = 5, request_delay: float = 0.1, timeout: int = 30):
        """
        Inicializa gerenciador de APIs.
        
        Args:
            max_workers: Número máximo de threads para processamento paralelo
            request_delay: Delay entre requisições (segundos) para evitar rate limiting
            timeout: Timeout em segundos para requisições HTTP
        """
        self.max_workers = max_workers
        self.request_delay = request_delay
        self.timeout = timeout
        
        # Sessões HTTP reutilizáveis
        self.iedb_session = requests.Session()
        self.iedb_session.headers.update({
            'User-Agent': 'Bioinformatics-Pipeline/2.0'
        })
        
        self.uniprot_session = requests.Session()
        self.uniprot_session.headers.update({
            'User-Agent': 'Bioinformatics-Pipeline/2.0'
        })
        
        logger.info(f"APIManager inicializado (workers: {max_workers}, delay: {request_delay}s)")
    
    def _predict_immunogenicity_iedb(self, peptide: str, allele: str) -> Optional[float]:
        """
        Prediz imunogenicidade usando IEDB Class I Immunogenicity Tool.
        
        Args:
            peptide: Sequência do peptídeo
            allele: Alelo HLA
            
        Returns:
            Score de imunogenicidade (0-1) ou None em caso de erro
        """
        try:
            url = "http://tools.iedb.org/immunogenicity/class1/"
            
            data = {
                'sequence_text': peptide,
                'allele': allele,
                'method': 'IEDB_recommended'
            }
            
            response = self.iedb_session.post(url, data=data, timeout=self.timeout)
            response.raise_for_status()
            
            # Parse da resposta
            try:
                result = response.json()
                if 'score' in result:
                    return float(result['score'])
                elif 'immunogenicity_score' in result:
                    return float(result['immunogenicity_score'])
            except:
                # Fallback: tenta extrair número da resposta HTML/texto
                import re
                text = response.text
                matches = re.findall(r'[\d.]+', text)
                if matches:
                    score = float(matches[0])
                    if 0 <= score <= 1:
                        return score
            
            return None
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"Erro na requisição IEDB para {peptide}: {e}")
            return None
        except Exception as e:
            logger.warning(f"Erro ao processar resposta IEDB para {peptide}: {e}")
            return None
    
    def _predict_affinity_iedb(self, peptide: str, allele: str, method: str = "netmhccons") -> Optional[float]:
        """
        Prediz afinidade de ligação MHC usando IEDB MHC-I Binding Prediction.
        
        Args:
            peptide: Sequência do peptídeo
            allele: Alelo HLA
            method: Método de predição
            
        Returns:
            Afinidade em nM (IC50) ou None em caso de erro
        """
        try:
            url = "http://tools.iedb.org/mhci/"
            
            data = {
                'sequence_text': peptide,
                'allele': allele,
                'method': method,
                'length': len(peptide)
            }
            
            response = self.iedb_session.post(url, data=data, timeout=self.timeout)
            response.raise_for_status()
            
            # Parse da resposta
            try:
                result = response.json()
                if 'ic50' in result:
                    return float(result['ic50'])
                elif 'affinity' in result:
                    return float(result['affinity'])
            except:
                # Fallback: tenta extrair IC50 do texto
                import re
                text = response.text
                matches = re.findall(r'(\d+\.?\d*)\s*nM', text, re.IGNORECASE)
                if matches:
                    return float(matches[0])
            
            return None
            
        except requests.exceptions.RequestException as e:
            logger.warning(f"Erro na requisição IEDB Affinity para {peptide}: {e}")
            return None
        except Exception as e:
            logger.warning(f"Erro ao processar resposta IEDB Affinity para {peptide}: {e}")
            return None
    
    def _search_uniprot(self, peptide: str) -> Optional[Dict]:
        """
        Busca informações sobre uma sequência peptídica no UniProt.
        
        Args:
            peptide: Sequência do peptídeo
            
        Returns:
            Dicionário com informações encontradas ou None
        """
        try:
            url = "https://rest.uniprot.org/uniprotkb/search"
            
            params = {
                'query': f'sequence:"{peptide}"',
                'format': 'json',
                'size': 1
            }
            
            response = self.uniprot_session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            
            data = response.json()
            
            if 'results' in data and len(data['results']) > 0:
                result = data['results'][0]
                return {
                    'uniprot_id': result.get('primaryAccession', ''),
                    'protein_name': result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                    'organism': result.get('organism', {}).get('scientificName', ''),
                    'found': True
                }
            else:
                return {'found': False}
                
        except requests.exceptions.RequestException as e:
            logger.warning(f"Erro na requisição UniProt para {peptide}: {e}")
            return None
        except Exception as e:
            logger.warning(f"Erro ao processar resposta UniProt para {peptide}: {e}")
            return None
    
    def _predict_immunogenicity_wrapper(self, args: Tuple[str, str]) -> Tuple[str, Optional[float]]:
        """Wrapper para predição de imunogenicidade."""
        peptide, allele = args
        time.sleep(self.request_delay)
        score = self._predict_immunogenicity_iedb(peptide, allele)
        return peptide, score
    
    def _predict_affinity_wrapper(self, args: Tuple[str, str, str]) -> Tuple[str, Optional[float]]:
        """Wrapper para predição de afinidade."""
        peptide, allele, method = args
        time.sleep(self.request_delay)
        affinity = self._predict_affinity_iedb(peptide, allele, method)
        return peptide, affinity
    
    def _search_uniprot_wrapper(self, peptide: str) -> Tuple[str, Optional[Dict]]:
        """Wrapper para busca UniProt."""
        time.sleep(self.request_delay)
        result = self._search_uniprot(peptide)
        return peptide, result
    
    def predict_immunogenicity_batch(self, peptides: List[str], allele: str) -> Dict[str, Optional[float]]:
        """
        Prediz imunogenicidade para múltiplos peptídeos em paralelo.
        
        Args:
            peptides: Lista de sequências de peptídeos
            allele: Alelo HLA
            
        Returns:
            Dicionário mapeando peptídeo -> score de imunogenicidade
        """
        results = {}
        
        try:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_peptide = {
                    executor.submit(self._predict_immunogenicity_wrapper, (pep, allele)): pep
                    for pep in peptides
                }
                
                for future in as_completed(future_to_peptide):
                    try:
                        peptide, score = future.result()
                        results[peptide] = score
                    except Exception as e:
                        peptide = future_to_peptide[future]
                        logger.error(f"Erro ao processar {peptide}: {e}")
                        results[peptide] = None
            
            logger.info(f"Predições de imunogenicidade concluídas para {len(peptides)} peptídeos")
            return results
            
        except Exception as e:
            logger.error(f"Erro no batch de imunogenicidade: {e}")
            return {pep: None for pep in peptides}
    
    def predict_affinity_batch(self, peptides: List[str], allele: str, 
                               method: str = "netmhccons") -> Dict[str, Optional[float]]:
        """
        Prediz afinidade para múltiplos peptídeos em paralelo.
        
        Args:
            peptides: Lista de sequências de peptídeos
            allele: Alelo HLA
            method: Método de predição
            
        Returns:
            Dicionário mapeando peptídeo -> afinidade (nM)
        """
        results = {}
        
        try:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_peptide = {
                    executor.submit(self._predict_affinity_wrapper, (pep, allele, method)): pep
                    for pep in peptides
                }
                
                for future in as_completed(future_to_peptide):
                    try:
                        peptide, affinity = future.result()
                        results[peptide] = affinity
                    except Exception as e:
                        peptide = future_to_peptide[future]
                        logger.error(f"Erro ao processar {peptide}: {e}")
                        results[peptide] = None
            
            logger.info(f"Predições de afinidade concluídas para {len(peptides)} peptídeos")
            return results
            
        except Exception as e:
            logger.error(f"Erro no batch de afinidade: {e}")
            return {pep: None for pep in peptides}
    
    def search_sequences_batch(self, peptides: List[str]) -> Dict[str, Optional[Dict]]:
        """
        Busca informações no UniProt para múltiplos peptídeos em paralelo.
        
        Args:
            peptides: Lista de sequências de peptídeos
            
        Returns:
            Dicionário mapeando peptídeo -> informações UniProt
        """
        results = {}
        
        try:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                future_to_peptide = {
                    executor.submit(self._search_uniprot_wrapper, pep): pep
                    for pep in peptides
                }
                
                for future in as_completed(future_to_peptide):
                    try:
                        peptide, info = future.result()
                        results[peptide] = info
                    except Exception as e:
                        peptide = future_to_peptide[future]
                        logger.error(f"Erro ao processar {peptide}: {e}")
                        results[peptide] = None
            
            logger.info(f"Buscas UniProt concluídas para {len(peptides)} peptídeos")
            return results
            
        except Exception as e:
            logger.error(f"Erro no batch UniProt: {e}")
            return {pep: None for pep in peptides}
    
    def enrich_dataframe(self, df: pd.DataFrame, allele: str, 
                        include_uniprot: bool = False) -> pd.DataFrame:
        """
        Enriquece DataFrame com dados de APIs externas.
        
        Args:
            df: DataFrame com coluna 'peptide'
            allele: Alelo HLA para predições
            include_uniprot: Se True, inclui busca no UniProt
            
        Returns:
            DataFrame enriquecido com colunas adicionais
        """
        if 'peptide' not in df.columns:
            raise ValueError("DataFrame deve conter coluna 'peptide'")
        
        df_result = df.copy()
        peptides = df_result['peptide'].tolist()
        
        logger.info(f"Enriquecendo dados para {len(peptides)} peptídeos...")
        
        # Predição de imunogenicidade (IEDB)
        try:
            logger.info("Consultando IEDB para imunogenicidade...")
            immunogenicity_scores = self.predict_immunogenicity_batch(peptides, allele)
            df_result['iedb_immunogenicity'] = df_result['peptide'].map(immunogenicity_scores)
        except Exception as e:
            logger.error(f"Erro ao enriquecer com imunogenicidade: {e}")
            df_result['iedb_immunogenicity'] = None
        
        # Predição de afinidade (IEDB) - apenas se não existir
        if 'affinity_nm' not in df_result.columns or df_result['affinity_nm'].isna().all():
            try:
                logger.info("Consultando IEDB para afinidade...")
                affinity_scores = self.predict_affinity_batch(peptides, allele)
                df_result['iedb_affinity_nm'] = df_result['peptide'].map(affinity_scores)
            except Exception as e:
                logger.error(f"Erro ao enriquecer com afinidade: {e}")
                df_result['iedb_affinity_nm'] = None
        
        # Busca UniProt (opcional)
        if include_uniprot:
            try:
                logger.info("Consultando UniProt...")
                uniprot_info = self.search_sequences_batch(peptides)
                df_result['uniprot_found'] = df_result['peptide'].map(
                    lambda p: uniprot_info.get(p, {}).get('found', False) if uniprot_info.get(p) else False
                )
                df_result['uniprot_id'] = df_result['peptide'].map(
                    lambda p: uniprot_info.get(p, {}).get('uniprot_id', '') if uniprot_info.get(p) else ''
                )
            except Exception as e:
                logger.error(f"Erro ao enriquecer com UniProt: {e}")
                df_result['uniprot_found'] = False
                df_result['uniprot_id'] = ''
        
        logger.info("Enriquecimento concluído")
        return df_result
