# -*- coding: utf-8 -*-
"""
Módulo de geração de relatórios PDF.
Usa fpdf2 para criar relatórios profissionais com resumo e tabelas de afinidade.
"""

from fpdf import FPDF
from typing import List, Dict, Optional
import pandas as pd
from datetime import datetime
import os
import tempfile
import matplotlib
matplotlib.use('Agg')  # Backend não-interativo
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from PIL import Image
import io


class PeptideReportPDF(FPDF):
    """Classe customizada para geração de relatórios de peptídeos."""
    
    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=15)
        self.set_margins(left=15, top=15, right=15)
        # Usa fontes padrão que suportam UTF-8
        # Para melhor suporte a caracteres especiais, instale fontes DejaVu separadamente
    
    def header(self):
        """Cabeçalho das páginas."""
        self.set_font('Helvetica', 'B', 16)
        self.set_text_color(44, 62, 80)  # #2C3E50
        self.cell(0, 10, 'Relatorio de Analise de Peptideos', 0, 1, 'L')
        self.set_font('Helvetica', '', 10)
        self.set_text_color(128, 128, 128)
        self.cell(0, 5, f'Gerado em: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}', 0, 1, 'L')
        self.ln(5)
    
    def footer(self):
        """Rodapé das páginas."""
        self.set_y(-15)
        self.set_font('Helvetica', '', 8)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, f'Pagina {self.page_no()}', 0, 0, 'C')
    
    def chapter_title(self, title: str):
        """Adiciona título de capítulo."""
        self.set_font('Helvetica', 'B', 14)
        self.set_text_color(44, 62, 80)
        self.cell(0, 10, title, 0, 1, 'L')
        self.ln(2)
        # Linha decorativa
        self.set_draw_color(44, 62, 80)
        self.set_line_width(0.5)
        self.line(15, self.get_y(), 195, self.get_y())
        self.ln(5)
    
    def add_text(self, text: str, font_size: int = 10, style: str = '', align: str = 'L'):
        """Adiciona texto formatado."""
        self.set_font('Helvetica', style, font_size)
        self.set_text_color(0, 0, 0)
        # Quebra texto em múltiplas linhas
        self.multi_cell(0, 5, text, align=align)
        self.ln(2)
    
    def add_summary_section(self, total_peptides: int, valid_peptides: int, 
                           top_affinity: float, avg_affinity: float):
        """Adiciona seção de resumo executivo."""
        self.chapter_title("Resumo Executivo")
        
        summary_text = f"""
Total de peptídeos analisados: {total_peptides}
Peptídeos válidos após validação: {valid_peptides}
Taxa de validação: {(valid_peptides/total_peptides*100):.1f}%

Melhor afinidade encontrada: {top_affinity:.2f} nM
Afinidade média: {avg_affinity:.2f} nM

Critérios de validação aplicados:
• Sequências contendo apenas aminoácidos canônicos (20 AA padrão)
• Tamanho entre 8-14 resíduos (compatível com MHC-I)
• Remoção de duplicatas
        """
        
        self.add_text(summary_text.strip(), font_size=10)
        self.ln(5)
    
    def add_affinity_table(self, df: pd.DataFrame, max_rows: int = 50):
        """Adiciona tabela de afinidade."""
        self.chapter_title("Tabela de Afinidade - Top Candidatos")
        
        # Prepara dados para tabela
        if 'affinity_nm' in df.columns:
            affinity_col = 'affinity_nm'
        elif 'iedb_affinity_nm' in df.columns:
            affinity_col = 'iedb_affinity_nm'
        elif 'affinity' in df.columns:
            affinity_col = 'affinity'
        else:
            self.add_text("AVISO: Coluna de afinidade nao encontrada no DataFrame.", font_size=10)
            return
        
        # Seleciona top N e colunas relevantes
        df_top = df.head(max_rows).copy()
        
        # Define colunas para exibir
        display_cols = ['peptide']
        if affinity_col in df_top.columns:
            display_cols.append(affinity_col)
        if 'mhc_score' in df_top.columns:
            display_cols.append('mhc_score')
        if 'presentation_score' in df_top.columns:
            display_cols.append('presentation_score')
        if 'iedb_immunogenicity' in df_top.columns:
            display_cols.append('iedb_immunogenicity')
        if 'percentile_rank' in df_top.columns:
            display_cols.append('percentile_rank')
        if 'final_rank_score' in df_top.columns:
            display_cols.append('final_rank_score')
        
        df_display = df_top[display_cols].copy()
        
        # Cabeçalho da tabela
        self.set_font('Helvetica', 'B', 9)
        self.set_fill_color(44, 62, 80)
        self.set_text_color(255, 255, 255)
        
        col_widths = [50, 30, 30, 30, 25, 25]  # Ajustar conforme necessário
        headers = list(df_display.columns)
        
        # Ajusta larguras conforme número de colunas
        if len(headers) <= 3:
            col_widths = [60, 50, 50]
        elif len(headers) == 4:
            col_widths = [50, 40, 40, 40]
        elif len(headers) == 5:
            col_widths = [45, 35, 35, 35, 30]
        else:
            col_widths = [40, 30, 30, 30, 25, 25]
        
        col_widths = col_widths[:len(headers)]
        total_width = sum(col_widths)
        
        # Desenha cabeçalho
        x_start = (210 - total_width) / 2  # Centraliza tabela
        x = x_start
        
        for i, header in enumerate(headers):
            # Traduz nomes de colunas para português
            header_pt = {
                'peptide': 'Peptideo',
                'affinity_nm': 'Afinidade (nM)',
                'iedb_affinity_nm': 'Afinidade IEDB (nM)',
                'mhc_score': 'Score MHC',
                'presentation_score': 'Score Apresentacao',
                'iedb_immunogenicity': 'Imunogenicidade',
                'percentile_rank': 'Rank %',
                'final_rank_score': 'Score Final'
            }.get(header, header)
            
            self.set_xy(x, self.get_y())
            self.cell(col_widths[i], 7, header_pt[:20], 1, 0, 'C', True)  # Limita tamanho do texto
            x += col_widths[i]
        
        self.ln(7)
        
        # Dados da tabela
        self.set_font('Helvetica', '', 8)
        self.set_text_color(0, 0, 0)
        self.set_fill_color(245, 245, 245)
        
        fill = False
        for idx, row in df_display.iterrows():
            if self.get_y() > 270:  # Nova página se necessário
                self.add_page()
                # Redesenha cabeçalho
                x = x_start
                self.set_font('Helvetica', 'B', 9)
                self.set_fill_color(44, 62, 80)
                self.set_text_color(255, 255, 255)
                for i, header in enumerate(headers):
                    header_pt = {
                        'peptide': 'Peptideo',
                        'affinity_nm': 'Afinidade (nM)',
                        'iedb_affinity_nm': 'Afinidade IEDB (nM)',
                        'mhc_score': 'Score MHC',
                        'presentation_score': 'Score Apresentacao',
                        'iedb_immunogenicity': 'Imunogenicidade',
                        'percentile_rank': 'Rank %',
                        'final_rank_score': 'Score Final'
                    }.get(header, header)
                    self.set_xy(x, self.get_y())
                    self.cell(col_widths[i], 7, header_pt[:20], 1, 0, 'C', True)
                    x += col_widths[i]
                self.ln(7)
                self.set_font('Helvetica', '', 8)
                self.set_text_color(0, 0, 0)
                self.set_fill_color(245, 245, 245)
            
            x = x_start
            for i, col in enumerate(headers):
                value = row[col]
                
                # Formata valores
                if pd.isna(value):
                    text = '-'
                elif isinstance(value, (int, float)):
                    if 'affinity' in col.lower():
                        text = f"{value:.2f}"
                    elif 'score' in col.lower() or 'rank' in col.lower():
                        text = f"{value:.3f}"
                    else:
                        text = f"{value:.2f}"
                else:
                    text = str(value)[:15]  # Limita tamanho
                
                self.set_xy(x, self.get_y())
                self.cell(col_widths[i], 6, text, 1, 0, 'C', fill)
                x += col_widths[i]
            
            self.ln(6)
            fill = not fill
        
        self.ln(5)
    
    def add_statistics_section(self, df: pd.DataFrame):
        """Adiciona seção de estatísticas."""
        self.chapter_title("Estatísticas Descritivas")
        
        stats_text = []
        
        # Estatísticas de afinidade
        if 'affinity_nm' in df.columns:
            aff_col = 'affinity_nm'
        elif 'iedb_affinity_nm' in df.columns:
            aff_col = 'iedb_affinity_nm'
        else:
            aff_col = None
        
        if aff_col and not df[aff_col].isna().all():
            aff_data = df[aff_col].dropna()
        stats_text.append(f"Afinidade (nM):")
        stats_text.append(f"  - Minima: {aff_data.min():.2f}")
        stats_text.append(f"  - Maxima: {aff_data.max():.2f}")
        stats_text.append(f"  - Media: {aff_data.mean():.2f}")
        stats_text.append(f"  - Mediana: {aff_data.median():.2f}")
        stats_text.append(f"  - Desvio Padrao: {aff_data.std():.2f}")
        stats_text.append("")
        
        # Estatísticas de tamanho
        if 'length' in df.columns:
            len_data = df['length']
            stats_text.append(f"Tamanho dos peptideos:")
            stats_text.append(f"  - Minimo: {len_data.min()} residuos")
            stats_text.append(f"  - Maximo: {len_data.max()} residuos")
            stats_text.append(f"  - Media: {len_data.mean():.1f} residuos")
            stats_text.append("")
        
        # Estatísticas de GRAVY
        if 'gravy' in df.columns:
            gravy_data = df['gravy'].dropna()
            stats_text.append(f"Indice GRAVY:")
            stats_text.append(f"  - Minimo: {gravy_data.min():.3f}")
            stats_text.append(f"  - Maximo: {gravy_data.max():.3f}")
            stats_text.append(f"  - Media: {gravy_data.mean():.3f}")
        
        if stats_text:
            self.add_text("\n".join(stats_text), font_size=10)
        else:
            self.add_text("Estatisticas nao disponiveis para os dados fornecidos.", font_size=10)


def generate_report(df: pd.DataFrame, output_path: str, 
                   total_peptides: Optional[int] = None,
                   include_statistics: bool = True,
                   max_table_rows: int = 50) -> str:
    """
    Gera relatório PDF completo.
    
    Args:
        df: DataFrame com resultados da análise
        output_path: Caminho para salvar o PDF
        total_peptides: Número total de peptídeos antes da validação
        include_statistics: Se True, inclui seção de estatísticas
        max_table_rows: Número máximo de linhas na tabela
        
    Returns:
        Caminho do arquivo gerado
    """
    if 'peptide' not in df.columns:
        raise ValueError("DataFrame deve conter coluna 'peptide'")
    
    # Cria PDF
    pdf = PeptideReportPDF()
    pdf.add_page()
    
    # Calcula estatísticas para resumo
    valid_peptides = len(df)
    if total_peptides is None:
        total_peptides = valid_peptides
    
    # Determina coluna de afinidade
    affinity_col = None
    for col in ['affinity_nm', 'iedb_affinity_nm', 'affinity']:
        if col in df.columns and not df[col].isna().all():
            affinity_col = col
            break
    
    if affinity_col:
        top_affinity = df[affinity_col].min()
        avg_affinity = df[affinity_col].mean()
    else:
        top_affinity = 0.0
        avg_affinity = 0.0
    
    # Adiciona seções
    pdf.add_summary_section(total_peptides, valid_peptides, top_affinity, avg_affinity)
    pdf.add_page()
    pdf.add_affinity_table(df, max_rows=max_table_rows)
    
    if include_statistics:
        pdf.add_page()
        pdf.add_statistics_section(df)
    
    # Salva PDF
    pdf.output(output_path)
    
    return output_path


def generate_simple_report(df: pd.DataFrame, output_path: str) -> str:
    """
    Gera relatório PDF simplificado (apenas resumo e tabela).
    
    Args:
        df: DataFrame com resultados
        output_path: Caminho para salvar o PDF
        
    Returns:
        Caminho do arquivo gerado
    """
    return generate_report(df, output_path, include_statistics=False, max_table_rows=30)


def create_affinity_chart_image(df: pd.DataFrame, max_peptides: int = 20) -> Optional[str]:
    """
    Cria gráfico de afinidade e salva como imagem temporária.
    
    Args:
        df: DataFrame com resultados
        max_peptides: Número máximo de peptídeos para exibir
        
    Returns:
        Caminho da imagem temporária ou None
    """
    if 'affinity_nm' not in df.columns or df['affinity_nm'].isna().all():
        return None
    
    try:
        df_plot = df.head(max_peptides).copy()
        df_plot = df_plot.sort_values('affinity_nm', ascending=True)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        colors = ['#27AE60' if x < 50 else '#F39C12' if x < 500 else '#E74C3C' 
                  for x in df_plot['affinity_nm']]
        
        bars = ax.barh(range(len(df_plot)), df_plot['affinity_nm'], color=colors)
        ax.set_yticks(range(len(df_plot)))
        ax.set_yticklabels(df_plot['peptide'], fontsize=8)
        ax.set_xlabel('Afinidade (nM) - Menor é Melhor', fontweight='bold')
        ax.set_title(f'Top {max_peptides} Peptídeos por Afinidade (IC50)', fontweight='bold', fontsize=12)
        ax.axvline(x=50, color='green', linestyle='--', linewidth=1.5, label='Limite Forte (50nM)')
        ax.legend()
        ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        
        # Salva em arquivo temporário
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.png')
        plt.savefig(temp_file.name, dpi=150, bbox_inches='tight')
        plt.close()
        
        return temp_file.name
    except Exception as e:
        print(f"Erro ao criar gráfico: {e}")
        return None


def add_image_to_pdf(pdf: PeptideReportPDF, image_path: str, width: float = 180):
    """
    Adiciona imagem ao PDF.
    
    Args:
        pdf: Instância do PDF
        image_path: Caminho da imagem
        width: Largura da imagem (mm)
    """
    try:
        if os.path.exists(image_path):
            pdf.image(image_path, x=(210 - width) / 2, w=width)
            pdf.ln(5)
    except Exception as e:
        print(f"Erro ao adicionar imagem ao PDF: {e}")


def generate_report_with_plots(df: pd.DataFrame, output_path: str,
                               total_peptides: Optional[int] = None,
                               include_statistics: bool = True,
                               max_table_rows: int = 50,
                               allele: str = "HLA-A*02:01") -> str:
    """
    Gera relatório PDF completo com gráficos.
    
    Args:
        df: DataFrame com resultados da análise
        output_path: Caminho para salvar o PDF
        total_peptides: Número total de peptídeos antes da validação
        include_statistics: Se True, inclui seção de estatísticas
        max_table_rows: Número máximo de linhas na tabela
        allele: Alelo HLA utilizado
        
    Returns:
        Caminho do arquivo gerado
    """
    if 'peptide' not in df.columns:
        raise ValueError("DataFrame deve conter coluna 'peptide'")
    
    # Cria PDF
    pdf = PeptideReportPDF()
    pdf.add_page()
    
    # Calcula estatísticas para resumo
    valid_peptides = len(df)
    if total_peptides is None:
        total_peptides = valid_peptides
    
    # Determina coluna de afinidade
    affinity_col = None
    for col in ['affinity_nm', 'iedb_affinity_nm', 'affinity']:
        if col in df.columns and not df[col].isna().all():
            affinity_col = col
            break
    
    if affinity_col:
        top_affinity = df[affinity_col].min()
        avg_affinity = df[affinity_col].mean()
    else:
        top_affinity = 0.0
        avg_affinity = 0.0
    
    # Adiciona seções
    pdf.add_summary_section(total_peptides, valid_peptides, top_affinity, avg_affinity)
    
    # Adiciona metadados
    pdf.chapter_title("Parametros da Analise")
    metadata_text = f"""
Alelo HLA utilizado: {allele}
Data de geracao: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
Total de peptideos processados: {total_peptides}
Peptideos validos: {valid_peptides}
    """
    pdf.add_text(metadata_text.strip(), font_size=10)
    
    # Adiciona gráfico de afinidade
    pdf.add_page()
    pdf.chapter_title("Grafico de Afinidade")
    chart_path = create_affinity_chart_image(df, max_peptides=20)
    if chart_path:
        add_image_to_pdf(pdf, chart_path, width=170)
        # Remove arquivo temporário após uso
        try:
            os.unlink(chart_path)
        except:
            pass
    
    # Tabela de afinidade
    pdf.add_page()
    pdf.add_affinity_table(df, max_rows=max_table_rows)
    
    # Estatísticas (se solicitado)
    if include_statistics:
        pdf.add_page()
        pdf.add_statistics_section(df)
    
    # Salva PDF
    pdf.output(output_path)
    
    return output_path
