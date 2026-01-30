import requests
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)

class APIClient:
    def __init__(self):
        self.iedb_url = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
        self.uniprot_url = "https://www.ebi.ac.uk/proteins/api/proteins/"

    def predict_mhcii_iedb(self, peptides, allele, method='nn_align'):
        """Faz predição de ligação MHC-II usando a IEDB API."""
        predictions = []
        try:
            # Formatar os dados para a requisição
            data = {
                'method': method,
                'sequence_text': "\n".join(peptides),
                'allele': allele,
                'length': len(peptides[0]) if peptides else 0
            }
            response = requests.post(self.iedb_url, data=data)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                # O formato da resposta é TSV, então vamos parsear
                for line in lines[1:]:  # Pular cabeçalho
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        predictions.append({
                            'peptide': parts[0],
                            'allele': parts[1],
                            'method': parts[2],
                            'ic50': float(parts[3]) if parts[3] != 'NA' else None
                        })
            else:
                logger.error(f"Erro na requisição à IEDB API: {response.status_code}")
        except Exception as e:
            logger.error(f"Exceção durante a predição MHC-II: {e}")
        return predictions

    def get_uniprot_annotation(self, protein_id):
        """Busca anotações do UniProt para uma proteína."""
        try:
            response = requests.get(f"{self.uniprot_url}{protein_id}")
            if response.status_code == 200:
                return response.json()
            else:
                logger.warning(f"Erro ao buscar anotações para {protein_id}: {response.status_code}")
                return None
        except Exception as e:
            logger.error(f"Exceção ao buscar anotações do UniProt: {e}")
            return None
