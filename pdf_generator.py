# -*- coding: utf-8 -*-
"""
Módulo de geração de relatórios PDF profissionais.
Usa fpdf2 para criar relatórios com gráficos e tabelas formatadas.
"""

import logging
import os
import tempfile
from typing import Optional
import pandas as pd
from datetime import datetime
from fpdf import FPDF

# Configuração de logging
logger = logging.getLogger(__name__)

# Tenta importar matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')  # Backend não-interativo
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("Matplotlib não disponível. Gráficos não serão incluídos no PDF.")


class PDFGenerator:
    """Classe para geração de relatórios PDF profissionais."""
    
    def __init__(self):
        """Inicializa gerador de PDF."""
        logger.info("PDFGenerator inicializado")
    
    def _create_affinity_chart_image(self, df: pd.DataFrame, max_peptides: int = 20) -> Optional[str]:
        """
        Cria gráfico de afinidade e salva como imagem temporária.
        
        Args:
            df: DataFrame com resultados
            max_peptides: Número máximo de peptídeos para exibir
            
        Returns:
            Caminho da imagem temporária ou None
        """
        if not MATPLOTLIB_AVAILABLE:
            return None
        
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
            
            temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.png')
            plt.savefig(temp_file.name, dpi=150, bbox_inches='tight')
            plt.close()
            
            return temp_file.name
        except Exception as e:
            logger.warning(f"Erro ao criar gráfico: {e}")
            return None
    
    def generate_report(self, df: pd.DataFrame, output_path: str,
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
        
        try:
            logger.info(f"Gerando relatório PDF: {output_path}")
            
            # Cria PDF
            pdf = PeptideReportPDF()
            pdf.add_page()
            
            # Calcula estatísticas
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
            
            # Metadados
            pdf.chapter_title("Parametros da Analise")
            metadata_text = f"""
Alelo HLA utilizado: {allele}
Data de geracao: {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}
Total de peptideos processados: {total_peptides}
Peptideos validos: {valid_peptides}
            """
            pdf.add_text(metadata_text.strip(), font_size=10)
            
            # Gráfico de afinidade
            pdf.add_page()
            pdf.chapter_title("Grafico de Afinidade")
            chart_path = self._create_affinity_chart_image(df, max_peptides=20)
            if chart_path:
                try:
                    pdf.image(chart_path, x=20, w=170)
                    pdf.ln(5)
                    # Remove arquivo temporário
                    os.unlink(chart_path)
                except Exception as e:
                    logger.warning(f"Erro ao adicionar gráfico ao PDF: {e}")
            
            # Tabela de afinidade
            pdf.add_page()
            pdf.add_affinity_table(df, max_rows=max_table_rows)
            
            # Estatísticas (se solicitado)
            if include_statistics:
                pdf.add_page()
                pdf.add_statistics_section(df)
            
            # Salva PDF
            pdf.output(output_path)
            logger.info(f"Relatório PDF gerado com sucesso: {output_path}")
            return output_path
            
        except Exception as e:
            logger.error(f"Erro ao gerar PDF: {e}")
            raise


class PeptideReportPDF(FPDF):
    """Classe customizada para geração de relatórios de peptídeos."""
    
    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=15)
        self.set_margins(left=15, top=15, right=15)
    
    def header(self):
        """Cabeçalho das páginas."""
        self.set_font('Helvetica', 'B', 16)
        self.set_text_color(44, 62, 80)
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
        self.set_draw_color(44, 62, 80)
        self.set_line_width(0.5)
        self.line(15, self.get_y(), 195, self.get_y())
        self.ln(5)
    
    def add_text(self, text: str, font_size: int = 10, style: str = '', align: str = 'L'):
        """Adiciona texto formatado."""
        self.set_font('Helvetica', style, font_size)
        self.set_text_color(0, 0, 0)
        self.multi_cell(0, 5, text, align=align)
        self.ln(2)
    
    def add_summary_section(self, total_peptides: int, valid_peptides: int,
                           top_affinity: float, avg_affinity: float):
        """Adiciona seção de resumo executivo."""
        self.chapter_title("Resumo Executivo")
        
        summary_text = f"""
Total de peptideos analisados: {total_peptides}
Peptideos validos apos validacao: {valid_peptides}
Taxa de validacao: {(valid_peptides/total_peptides*100):.1f}%

Melhor afinidade encontrada: {top_affinity:.2f} nM
Afinidade media: {avg_affinity:.2f} nM

Criterios de validacao aplicados:
• Sequencias contendo apenas aminoacidos canonicos (20 AA padrao)
• Tamanho entre 8-14 residuos (compativel com MHC-I)
• Remocao de duplicatas
        """
        
        self.add_text(summary_text.strip(), font_size=10)
        self.ln(5)
    
    def add_affinity_table(self, df: pd.DataFrame, max_rows: int = 50):
        """Adiciona tabela de afinidade."""
        self.chapter_title("Tabela de Afinidade - Top Candidatos")
        
        # Determina coluna de afinidade
        if 'affinity_nm' in df.columns:
            affinity_col = 'affinity_nm'
        elif 'iedb_affinity_nm' in df.columns:
            affinity_col = 'iedb_affinity_nm'
        elif 'affinity' in df.columns:
            affinity_col = 'affinity'
        else:
            self.add_text("AVISO: Coluna de afinidade nao encontrada.", font_size=10)
            return
        
        df_top = df.head(max_rows).copy()
        
        # Colunas para exibir
        display_cols = ['peptide']
        if affinity_col in df_top.columns:
            display_cols.append(affinity_col)
        if 'mhc_score' in df_top.columns:
            display_cols.append('mhc_score')
        if 'iedb_immunogenicity' in df_top.columns:
            display_cols.append('iedb_immunogenicity')
        if 'percentile_rank' in df_top.columns:
            display_cols.append('percentile_rank')
        
        df_display = df_top[display_cols].copy()
        
        # Cabeçalho
        self.set_font('Helvetica', 'B', 9)
        self.set_fill_color(44, 62, 80)
        self.set_text_color(255, 255, 255)
        
        col_widths = [50, 30, 30, 25, 25]
        headers = list(df_display.columns)
        col_widths = col_widths[:len(headers)]
        total_width = sum(col_widths)
        x_start = (210 - total_width) / 2
        x = x_start
        
        header_map = {
            'peptide': 'Peptideo',
            'affinity_nm': 'Afinidade (nM)',
            'iedb_affinity_nm': 'Afinidade IEDB (nM)',
            'mhc_score': 'Score MHC',
            'iedb_immunogenicity': 'Imunogenicidade',
            'percentile_rank': 'Rank %'
        }
        
        for i, header in enumerate(headers):
            header_pt = header_map.get(header, header)
            self.set_xy(x, self.get_y())
            self.cell(col_widths[i], 7, header_pt[:20], 1, 0, 'C', True)
            x += col_widths[i]
        
        self.ln(7)
        
        # Dados
        self.set_font('Helvetica', '', 8)
        self.set_text_color(0, 0, 0)
        self.set_fill_color(245, 245, 245)
        
        fill = False
        for idx, row in df_display.iterrows():
            if self.get_y() > 270:
                self.add_page()
                x = x_start
                self.set_font('Helvetica', 'B', 9)
                self.set_fill_color(44, 62, 80)
                self.set_text_color(255, 255, 255)
                for i, header in enumerate(headers):
                    header_pt = header_map.get(header, header)
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
                    text = str(value)[:15]
                
                self.set_xy(x, self.get_y())
                self.cell(col_widths[i], 6, text, 1, 0, 'C', fill)
                x += col_widths[i]
            
            self.ln(6)
            fill = not fill
        
        self.ln(5)
    
    def add_statistics_section(self, df: pd.DataFrame):
        """Adiciona seção de estatísticas."""
        self.chapter_title("Estatisticas Descritivas")
        
        stats_text = []
        
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
        
        if 'length' in df.columns:
            len_data = df['length']
            stats_text.append(f"Tamanho dos peptideos:")
            stats_text.append(f"  - Minimo: {len_data.min()} residuos")
            stats_text.append(f"  - Maximo: {len_data.max()} residuos")
            stats_text.append(f"  - Media: {len_data.mean():.1f} residuos")
        
        if stats_text:
            self.add_text("\n".join(stats_text), font_size=10)
        else:
            self.add_text("Estatisticas nao disponiveis.", font_size=10)
