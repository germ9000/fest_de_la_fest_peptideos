from fpdf import FPDF
import pandas as pd
import logging
from datetime import datetime
import tempfile
import os

logger = logging.getLogger(__name__)

class PDFGenerator(FPDF):
    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=15)

    def header(self):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, 'Relatório de Análise de Epítopos Imunológicos', 0, 1, 'C')

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, f'Página {self.page_no()}', 0, 0, 'C')

    def add_title(self, title):
        self.set_font('Arial', 'B', 16)
        self.cell(0, 10, title, 0, 1, 'L')
        self.ln(5)

    def add_section_title(self, title):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, title, 0, 1, 'L')
        self.ln(2)

    def add_text(self, text):
        self.set_font('Arial', '', 10)
        self.multi_cell(0, 5, text)
        self.ln(2)

    def add_table(self, df):
        # Adiciona uma tabela ao PDF
        self.set_font('Arial', '', 8)
        # Cabeçalho
        for col in df.columns:
            self.cell(40, 6, col, 1, 0, 'C')
        self.ln()
        # Dados
        for _, row in df.iterrows():
            for col in df.columns:
                self.cell(40, 6, str(row[col]), 1, 0, 'C')
            self.ln()
        self.ln(5)

    def generate_report(self, fasta_data, mhci_results, mhcii_results, alleles_mhci, alleles_mhcii,
                        peptide_lengths_mhci, peptide_lengths_mhcii, ic50_threshold, percentile_threshold):
        # Criar um arquivo temporário para o PDF
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdf')
        pdf_path = temp_file.name

        self.add_page()

        # Título
        self.add_title("Relatório de Análise de Epítopos Imunológicos")
        self.add_text(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        # Seção: Parâmetros
        self.add_section_title("Parâmetros da Análise")
        params_text = f"""
        Alelos MHC-I: {', '.join(alleles_mhci)}
        Alelos MHC-II: {', '.join(alleles_mhcii)}
        Comprimentos dos peptídeos MHC-I: {', '.join(map(str, peptide_lengths_mhci))}
        Comprimentos dos peptídeos MHC-II: {', '.join(map(str, peptide_lengths_mhcii))}
        Threshold de IC50: {ic50_threshold} nM
        Threshold de Percentil: {percentile_threshold}
        """
        self.add_text(params_text)

        # Seção: Sequências Analisadas
        self.add_section_title("Sequências Analisadas")
        if fasta_data:
            fasta_df = pd.DataFrame(fasta_data)
            self.add_table(fasta_df[['id', 'length', 'isoelectric_point', 'gravy']].head())  # Mostrar apenas as primeiras
        else:
            self.add_text("Nenhuma sequência analisada.")

        # Seção: Resultados MHC-I
        self.add_section_title("Resultados MHC-I")
        if mhci_results:
            # Combinar todos os resultados em um DataFrame
            df_list = []
            for result in mhci_results:
                df_list.append(result['predictions'])
            combined_df = pd.concat(df_list, ignore_index=True)
            # Filtrar por threshold
            filtered_df = combined_df[(combined_df['ic50'] <= ic50_threshold) & (combined_df['percentile'] <= percentile_threshold)]
            if not filtered_df.empty:
                self.add_text(f"Total de epítopos preditos (com thresholds): {len(filtered_df)}")
                # Mostrar os top 20
                top_df = filtered_df.sort_values(by='ic50').head(20)
                self.add_table(top_df[['peptide', 'allele', 'ic50', 'percentile']])
            else:
                self.add_text("Nenhum epítopo MHC-I encontrado com os thresholds fornecidos.")
        else:
            self.add_text("Nenhum resultado MHC-I disponível.")

        # Seção: Resultados MHC-II
        self.add_section_title("Resultados MHC-II")
        if mhcii_results is not None and not mhcii_results.empty:
            self.add_text(f"Total de epítopos MHC-II preditos: {len(mhcii_results)}")
            # Mostrar os top 20
            top_df = mhcii_results.sort_values(by='ic50').head(20)
            self.add_table(top_df[['peptide', 'allele', 'ic50']])
        else:
            self.add_text("Nenhum resultado MHC-II disponível.")

        # Gerar o PDF
        self.output(pdf_path)
        logger.info(f"Relatório gerado: {pdf_path}")
        return pdf_path
