# -*- coding: utf-8 -*-
"""
Módulo cliente para APIs externas (IEDB, UniProt).
Implementa processamento paralelo usando ThreadPoolExecutor para otimizar performance.
"""

import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Optional, Tuple
import pandas as pd
import numpy as np

# Tenta importar Streamlit para cache (opcional)
try:
    import streamlit as st
    STREAMLIT_AVAILABLE = True
    
    def cached_function(func):
        """Wrapper para cache do Streamlit."""
        return st.cache_data(ttl=3600)(func)
except ImportError:
    STREAMLIT_AVAILABLE = False
    def cached_function(func):
        """Wrapper dummy quando Streamlit não está disponível."""
        return func


class IEDBClient:
    """Cliente para API do IEDB (Immune Epitope Database)."""
    
    BASE_URL = "http://tools.iedb.org"
    
    def __init__(self, timeout: int = 30):
        """
        Inicializa cliente IEDB.
        
        Args:
            timeout: Timeout em segundos para requisições HTTP
        """
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Bioinformatics-Pipeline/1.0'
        })
    
    def predict_immunogenicity(self, peptide: str, allele: str = "HLA-A*02:01") -> Optional[float]:
        """
        Prediz imunogenicidade usando IEDB Class I Immunogenicity Tool.
        
        Args:
            peptide: Sequência do peptídeo
            allele: Alelo HLA (padrão: HLA-A*02:01)
            
        Returns:
            Score de imunogenicidade (0-1) ou None em caso de erro
        """
        try:
            # Endpoint do IEDB Immunogenicity Tool
            url = f"{self.BASE_URL}/immunogenicity/class1/"
            
            data = {
                'sequence_text': peptide,
                'allele': allele,
                'method': 'IEDB_recommended'
            }
            
            response = self.session.post(url, data=data, timeout=self.timeout)
            response.raise_for_status()
            
            # Parse da resposta (formato pode variar)
            # Em produção, ajustar conforme documentação atual da API
            try:
                result = response.json()
                if 'score' in result:
                    return float(result['score'])
                elif 'immunogenicity_score' in result:
                    return float(result['immunogenicity_score'])
            except:
                # Fallback: tenta extrair número da resposta HTML/texto
                text = response.text
                # Procura por padrões numéricos (score entre 0-1)
                import re
                matches = re.findall(r'[\d.]+', text)
                if matches:
                    score = float(matches[0])
                    if 0 <= score <= 1:
                        return score
            
            # Se não conseguir parsear, retorna None
            return None
            
        except requests.exceptions.RequestException as e:
            print(f"Erro na requisição IEDB para {peptide}: {e}")
            return None
        except Exception as e:
            print(f"Erro ao processar resposta IEDB para {peptide}: {e}")
            return None
    
    def predict_affinity(self, peptide: str, allele: str = "HLA-A*02:01", method: str = "netmhccons") -> Optional[float]:
        """
        Prediz afinidade de ligação MHC usando IEDB MHC-I Binding Prediction.
        
        Args:
            peptide: Sequência do peptídeo
            allele: Alelo HLA
            method: Método de predição (netmhccons, netmhcpan, etc.)
            
        Returns:
            Afinidade em nM (IC50) ou None em caso de erro
        """
        try:
            url = f"{self.BASE_URL}/mhci/"
            
            data = {
                'sequence_text': peptide,
                'allele': allele,
                'method': method,
                'length': len(peptide)
            }
            
            response = self.session.post(url, data=data, timeout=self.timeout)
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
                # Procura por valores numéricos seguidos de "nM" ou "IC50"
                matches = re.findall(r'(\d+\.?\d*)\s*nM', text, re.IGNORECASE)
                if matches:
                    return float(matches[0])
            
            return None
            
        except requests.exceptions.RequestException as e:
            print(f"Erro na requisição IEDB Affinity para {peptide}: {e}")
            return None
        except Exception as e:
            print(f"Erro ao processar resposta IEDB Affinity para {peptide}: {e}")
            return None


class UniProtClient:
    """Cliente para API do UniProt."""
    
    BASE_URL = "https://rest.uniprot.org"
    
    def __init__(self, timeout: int = 30):
        """
        Inicializa cliente UniProt.
        
        Args:
            timeout: Timeout em segundos para requisições HTTP
        """
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Bioinformatics-Pipeline/1.0'
        })
    
    def search_sequence(self, peptide: str) -> Optional[Dict]:
        """
        Busca informações sobre uma sequência peptídica no UniProt.
        
        Args:
            peptide: Sequência do peptídeo
            
        Returns:
            Dicionário com informações encontradas ou None
        """
        try:
            url = f"{self.BASE_URL}/uniprotkb/search"
            
            params = {
                'query': f'sequence:"{peptide}"',
                'format': 'json',
                'size': 1
            }
            
            response = self.session.get(url, params=params, timeout=self.timeout)
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
            print(f"Erro na requisição UniProt para {peptide}: {e}")
            return None
        except Exception as e:
            print(f"Erro ao processar resposta UniProt para {peptide}: {e}")
            return None
    
    def get_protein_info(self, uniprot_id: str) -> Optional[Dict]:
        """
        Obtém informações detalhadas de uma proteína pelo ID UniProt.
        
        Args:
            uniprot_id: ID UniProt (ex: P12345)
            
        Returns:
            Dicionário com informações da proteína ou None
        """
        try:
            url = f"{self.BASE_URL}/uniprotkb/{uniprot_id}"
            
            params = {'format': 'json'}
            
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            
            return response.json()
            
        except requests.exceptions.RequestException as e:
            print(f"Erro ao buscar proteína {uniprot_id}: {e}")
            return None
        except Exception as e:
            print(f"Erro ao processar resposta UniProt para {uniprot_id}: {e}")
            return None


class APIClientManager:
    """Gerenciador de clientes de API com processamento paralelo."""
    
    def __init__(self, max_workers: int = 5, request_delay: float = 0.1):
        """
        Inicializa gerenciador de APIs.
        
        Args:
            max_workers: Número máximo de threads para processamento paralelo
            request_delay: Delay entre requisições (segundos) para evitar rate limiting
        """
        self.max_workers = max_workers
        self.request_delay = request_delay
        self.iedb_client = IEDBClient()
        self.uniprot_client = UniProtClient()
    
    def _predict_immunogenicity_wrapper(self, args: Tuple[str, str]) -> Tuple[str, Optional[float]]:
        """Wrapper para predição de imunogenicidade."""
        peptide, allele = args
        time.sleep(self.request_delay)  # Rate limiting
        score = self.iedb_client.predict_immunogenicity(peptide, allele)
        return peptide, score
    
    def _predict_affinity_wrapper(self, args: Tuple[str, str, str]) -> Tuple[str, Optional[float]]:
        """Wrapper para predição de afinidade."""
        peptide, allele, method = args
        time.sleep(self.request_delay)
        affinity = self.iedb_client.predict_affinity(peptide, allele, method)
        return peptide, affinity
    
    def _search_sequence_wrapper(self, peptide: str) -> Tuple[str, Optional[Dict]]:
        """Wrapper para busca de sequência."""
        time.sleep(self.request_delay)
        result = self.uniprot_client.search_sequence(peptide)
        return peptide, result
    
    def predict_immunogenicity_batch(self, peptides: List[str], allele: str = "HLA-A*02:01") -> Dict[str, Optional[float]]:
        """
        Prediz imunogenicidade para múltiplos peptídeos em paralelo.
        
        Args:
            peptides: Lista de sequências de peptídeos
            allele: Alelo HLA
            
        Returns:
            Dicionário mapeando peptídeo -> score de imunogenicidade
        """
        results = {}
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submete todas as tarefas
            future_to_peptide = {
                executor.submit(self._predict_immunogenicity_wrapper, (pep, allele)): pep
                for pep in peptides
            }
            
            # Coleta resultados conforme completam
            for future in as_completed(future_to_peptide):
                try:
                    peptide, score = future.result()
                    results[peptide] = score
                except Exception as e:
                    peptide = future_to_peptide[future]
                    print(f"Erro ao processar {peptide}: {e}")
                    results[peptide] = None
        
        return results
    
    def predict_affinity_batch(self, peptides: List[str], allele: str = "HLA-A*02:01", 
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
                    print(f"Erro ao processar {peptide}: {e}")
                    results[peptide] = None
        
        return results
    
    def search_sequences_batch(self, peptides: List[str]) -> Dict[str, Optional[Dict]]:
        """
        Busca informações no UniProt para múltiplos peptídeos em paralelo.
        
        Args:
            peptides: Lista de sequências de peptídeos
            
        Returns:
            Dicionário mapeando peptídeo -> informações UniProt
        """
        results = {}
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_peptide = {
                executor.submit(self._search_sequence_wrapper, pep): pep
                for pep in peptides
            }
            
            for future in as_completed(future_to_peptide):
                try:
                    peptide, info = future.result()
                    results[peptide] = info
                except Exception as e:
                    peptide = future_to_peptide[future]
                    print(f"Erro ao processar {peptide}: {e}")
                    results[peptide] = None
        
        return results
    
    def enrich_dataframe(self, df: pd.DataFrame, allele: str = "HLA-A*02:01", 
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
        
        print(f"Enriquecendo dados para {len(peptides)} peptídeos...")
        
        # Predição de imunogenicidade (IEDB)
        print("  → Predizendo imunogenicidade (IEDB)...")
        immunogenicity_scores = self.predict_immunogenicity_batch(peptides, allele)
        df_result['iedb_immunogenicity'] = df_result['peptide'].map(immunogenicity_scores)
        
        # Predição de afinidade (IEDB) - apenas se não existir
        if 'affinity_nm' not in df_result.columns or df_result['affinity_nm'].isna().all():
            print("  → Predizendo afinidade (IEDB)...")
            affinity_scores = self.predict_affinity_batch(peptides, allele)
            df_result['iedb_affinity_nm'] = df_result['peptide'].map(affinity_scores)
        
        # Busca UniProt (opcional)
        if include_uniprot:
            print("  → Buscando informações no UniProt...")
            uniprot_info = self.search_sequences_batch(peptides)
            df_result['uniprot_found'] = df_result['peptide'].map(
                lambda p: uniprot_info.get(p, {}).get('found', False) if uniprot_info.get(p) else False
            )
            df_result['uniprot_id'] = df_result['peptide'].map(
                lambda p: uniprot_info.get(p, {}).get('uniprot_id', '') if uniprot_info.get(p) else ''
            )
        
        print("✅ Enriquecimento concluído.")
        
        return df_result
