import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)

class MHCAnalyzer:
    def __init__(self, predictor):
        self.predictor = predictor

    def predict_mhci_epitopes(self, sequences, alleles, peptide_lengths):
        """Prediz epítopos MHC-I para uma lista de sequências."""
        all_results = []
        try:
            for header, seq in sequences:
                # Gerar peptídeos para cada comprimento
                peptides = []
                for length in peptide_lengths:
                    for i in range(len(seq) - length + 1):
                        peptides.append(str(seq[i:i+length]))
                # Fazer a predição em lotes para melhor performance
                if peptides:
                    # Usar ThreadPoolExecutor para paralelizar as predições
                    with ThreadPoolExecutor(max_workers=4) as executor:
                        futures = []
                        for allele in alleles:
                            futures.append(executor.submit(self.predictor.predict, peptides, alleles=[allele]))
                        for future in as_completed(futures):
                            result = future.result()
                            if not result.empty:
                                all_results.append({
                                    'header': header,
                                    'allele': allele,
                                    'predictions': result
                                })
            logger.info(f"Total de resultados MHC-I: {len(all_results)}")
        except Exception as e:
            logger.error(f"Erro durante a predição MHC-I: {e}")
        return all_results
