# -*- coding: utf-8 -*-
"""
Aplicação Streamlit - Pipeline de Análise de Peptídeos
=======================================================
Interface web completa para análise de epítopos e geração de dossiês.

Autor: Engenheiro de Bioinformática Sênior
Versão: 2.0.0
"""

import logging
import sys
from typing import Optional
import streamlit as st
import pandas as pd
import numpy as np
import io
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime

# Configuração de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Importa módulos da src/
try:
    from src import (
        FastaProcessor,
        MHCAnalyzer,
        APIManager,
        PDFGenerator,
        ConservationAnalyzer
    )
except ImportError as e:
    st.error(f"❌ Erro ao importar módulos: {e}")
    st.error("Certifique-se de que todos os arquivos estão na pasta src/")
    st.stop()

# Configuração da página
st.set_page_config(
    page_title="Pipeline de Análise de Peptídeos",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# CSS customizado
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #2C3E50;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)


@st.cache_resource
def get_mhc_analyzer():
    """
    Carrega analisador MHC uma única vez na memória.
    Usa @st.cache_resource para evitar recarregar modelos pesados.
    """
    logger.info("Carregando MHCAnalyzer (cacheado)")
    return MHCAnalyzer()


@st.cache_resource
def get_fasta_processor(min_length: int = 8, max_length: int = 14):
    """
    Carrega processador FASTA com parâmetros específicos.
    """
    return FastaProcessor(min_length=min_length, max_length=max_length)


def create_affinity_chart(df: pd.DataFrame) -> Optional[go.Figure]:
    """Cria gráfico de afinidade interativo."""
    if 'affinity_nm' not in df.columns or df['affinity_nm'].isna().all():
        return None
    
    df_plot = df.head(20).copy()
    df_plot = df_plot.sort_values('affinity_nm', ascending=True)
    
    fig = go.Figure()
    
    colors = ['#27AE60' if x < 50 else '#F39C12' if x < 500 else '#E74C3C' 
              for x in df_plot['affinity_nm']]
    
    fig.add_trace(go.Bar(
        x=df_plot['affinity_nm'],
        y=df_plot['peptide'],
        orientation='h',
        marker_color=colors,
        text=[f"{x:.2f} nM" for x in df_plot['affinity_nm']],
        textposition='outside',
        name='Afinidade'
    ))
    
    fig.add_vline(x=50, line_dash="dash", line_color="green", 
                  annotation_text="Limite Forte (50nM)")
    
    fig.update_layout(
        title="Top 20 Peptídeos por Afinidade (IC50)",
        xaxis_title="Afinidade (nM) - Menor é Melhor",
        yaxis_title="Peptídeo",
        height=600,
        showlegend=False,
        template="plotly_white"
    )
    
    return fig


def main():
    """Função principal da aplicação."""
    
    # Header
    st.markdown('<p class="main-header">🧬 Pipeline de Análise de Peptídeos</p>', 
                unsafe_allow_html=True)
    st.markdown("---")
    
    # Sidebar - Configurações
    with st.sidebar:
        st.header("⚙️ Configurações")
        
        # Upload de arquivo
        st.subheader("📁 Entrada de Dados")
        uploaded_file = st.file_uploader(
            "Carregue arquivo FASTA, TSV, CSV ou Excel",
            type=['fasta', 'fa', 'faa', 'tsv', 'csv', 'txt', 'xlsx', 'xls'],
            help="Formatos suportados: FASTA, TSV, CSV, Excel"
        )
        
        # Parâmetros de análise
        st.subheader("🔬 Parâmetros de Análise")
        
        # Tipo de MHC
        mhc_type = st.radio(
            "Tipo de MHC",
            ["MHC-I", "MHC-II", "Ambos"],
            index=0,
            help="Selecione o tipo de MHC para análise"
        )
        
        # Alelo HLA-I
        allele_class1_options = [
            "HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02", 
            "HLA-B*08:01", "HLA-C*07:01", "HLA-A*01:01"
        ]
        selected_allele_class1 = st.selectbox(
            "Alelo HLA-I (MHC-I)",
            allele_class1_options,
            index=0,
            disabled=(mhc_type == "MHC-II")
        )
        
        # Alelo HLA-II
        allele_class2_options = [
            "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01",
            "HLA-DRB1*07:01", "HLA-DRB1*11:01", "HLA-DRB1*15:01"
        ]
        selected_allele_class2 = st.selectbox(
            "Alelo HLA-II (MHC-II)",
            allele_class2_options,
            index=0,
            disabled=(mhc_type == "MHC-I")
        )
        
        # Tamanho do k-mer
        min_length = st.slider("Tamanho mínimo (resíduos)", 8, 14, 8)
        max_length = st.slider("Tamanho máximo (resíduos)", 8, 14, 14)
        
        # Opções avançadas
        st.subheader("🔧 Opções Avançadas")
        
        use_api_enrichment = st.checkbox(
            "Enriquecer com APIs externas (IEDB/UniProt)",
            value=False,
            help="Consulta APIs externas (pode ser mais lento)"
        )
        
        include_conservation = st.checkbox(
            "Análise de conservação",
            value=False
        )
        
        max_workers = st.slider("Threads paralelas", 1, 10, 5)
        
        include_uniprot = st.checkbox("Incluir busca UniProt", value=False)
        
        # Botão de processamento
        st.markdown("---")
        process_button = st.button(
            "🚀 Processar Análise",
            type="primary",
            use_container_width=True
        )
    
    # Área principal
    if not uploaded_file and not process_button:
        st.info("👈 Por favor, carregue um arquivo e configure os parâmetros na barra lateral.")
        return
    
    if not uploaded_file:
        st.warning("⚠️ Por favor, carregue um arquivo antes de processar.")
        return
    
    # Processamento
    if process_button:
        with st.spinner("🔄 Processando análise... Isso pode levar alguns minutos."):
            try:
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                # Etapa 1: Carregamento e Validação
                status_text.text("📥 Carregando arquivo...")
                progress_bar.progress(10)
                
                file_bytes = uploaded_file.read()
                processor = get_fasta_processor(min_length, max_length)
                df_input = processor.load_peptides_from_bytes(file_bytes, uploaded_file.name)
                total_peptides = len(df_input)
                
                status_text.text(f"✅ {total_peptides} peptídeos carregados. Validando...")
                progress_bar.progress(30)
                
                df_valid = processor.validate_peptides(df_input)
                valid_count = len(df_valid)
                
                status_text.text(f"✅ {valid_count} peptídeos válidos. Calculando propriedades...")
                progress_bar.progress(50)
                
                df_result = processor.add_physchem_properties(df_valid)
                progress_bar.progress(60)
                
                # Etapa 2: Predições MHC
                status_text.text("🧬 Executando predições MHC...")
                progress_bar.progress(70)
                
                mhc_analyzer = get_mhc_analyzer()
                
                if mhc_type in ["MHC-I", "Ambos"]:
                    df_result = mhc_analyzer.predict_class1(df_result, [selected_allele_class1])
                
                if mhc_type in ["MHC-II", "Ambos"]:
                    df_result = mhc_analyzer.predict_class2(df_result, [selected_allele_class2])
                
                progress_bar.progress(75)
                
                # Etapa 3: Conservação (opcional)
                if include_conservation:
                    status_text.text("🔬 Calculando conservação...")
                    progress_bar.progress(80)
                    conservation_analyzer = ConservationAnalyzer()
                    df_result = conservation_analyzer.add_conservation_to_dataframe(df_result)
                
                # Etapa 4: APIs (opcional)
                if use_api_enrichment:
                    status_text.text("🌐 Consultando APIs externas...")
                    progress_bar.progress(85)
                    api_manager = APIManager(max_workers=max_workers, request_delay=0.2)
                    selected_allele = selected_allele_class1 if mhc_type != "MHC-II" else selected_allele_class2
                    df_result = api_manager.enrich_dataframe(
                        df_result,
                        allele=selected_allele,
                        include_uniprot=include_uniprot
                    )
                
                # Etapa 5: Scores Finais
                status_text.text("📊 Calculando scores finais...")
                progress_bar.progress(90)
                
                df_result = mhc_analyzer.calculate_final_scores(df_result)
                
                progress_bar.progress(100)
                status_text.text("✅ Análise concluída!")
                
                # Salva no session state
                st.session_state['df_result'] = df_result
                st.session_state['total_peptides'] = total_peptides
                st.session_state['valid_peptides'] = valid_count
                st.session_state['selected_allele'] = selected_allele_class1 if mhc_type != "MHC-II" else selected_allele_class2
                st.session_state['mhc_type'] = mhc_type
                
                st.success("✅ Análise processada com sucesso!")
                
            except Exception as e:
                logger.error(f"Erro durante processamento: {e}", exc_info=True)
                st.error(f"❌ Erro durante o processamento: {str(e)}")
                st.exception(e)
                return
    
    # Verifica se há resultados
    if 'df_result' not in st.session_state:
        return
    
    df_result = st.session_state['df_result']
    total_peptides = st.session_state.get('total_peptides', len(df_result))
    valid_peptides = st.session_state.get('valid_peptides', len(df_result))
    
    # Tabs - Visualização
    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Resumo", 
        "🧬 Sequências", 
        "🎯 Análise de Epítopos", 
        "📄 Relatório Final"
    ])
    
    with tab1:
        st.header("📊 Resumo da Análise")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total de Peptídeos", total_peptides)
        with col2:
            st.metric("Peptídeos Válidos", valid_peptides)
        with col3:
            if 'affinity_nm' in df_result.columns:
                best_affinity = df_result['affinity_nm'].min()
                st.metric("Melhor Afinidade", f"{best_affinity:.2f} nM")
        with col4:
            validation_rate = (valid_peptides / total_peptides * 100) if total_peptides > 0 else 0
            st.metric("Taxa de Validação", f"{validation_rate:.1f}%")
        
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Top 20 por Afinidade")
            affinity_chart = create_affinity_chart(df_result)
            if affinity_chart:
                st.plotly_chart(affinity_chart, use_container_width=True)
        
        with col2:
            if 'gravy' in df_result.columns and 'pI' in df_result.columns:
                scatter_chart = px.scatter(
                    df_result.head(50),
                    x='gravy',
                    y='pI',
                    size='affinity_nm' if 'affinity_nm' in df_result.columns else None,
                    color='affinity_nm' if 'affinity_nm' in df_result.columns else None,
                    hover_data=['peptide'],
                    title="Propriedades Físico-Químicas"
                )
                st.plotly_chart(scatter_chart, use_container_width=True)
    
    with tab2:
        st.header("🧬 Sequências Analisadas")
        
        search_term = st.text_input("🔍 Buscar peptídeo", "")
        if search_term:
            df_display = df_result[df_result['peptide'].str.contains(search_term, case=False, na=False)]
        else:
            df_display = df_result
        
        display_cols = ['peptide']
        if 'affinity_nm' in df_display.columns:
            display_cols.append('affinity_nm')
        if 'mhc_score' in df_display.columns:
            display_cols.append('mhc_score')
        
        st.dataframe(df_display[display_cols], use_container_width=True, height=600)
    
    with tab3:
        st.header("🎯 Análise Detalhada de Epítopos")
        top10 = df_result.head(10)
        st.dataframe(top10, use_container_width=True)
    
    with tab4:
        st.header("📄 Relatório Final")
        
        report_name = st.text_input(
            "Nome do arquivo PDF",
            value=f"relatorio_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
        )
        
        if st.button("📥 Gerar e Baixar Relatório PDF", type="primary"):
            with st.spinner("Gerando relatório PDF..."):
                try:
                    pdf_generator = PDFGenerator()
                    pdf_path = pdf_generator.generate_report(
                        df_result,
                        report_name,
                        total_peptides=total_peptides,
                        include_statistics=True,
                        max_table_rows=50,
                        allele=st.session_state.get('selected_allele', selected_allele_class1)
                    )
                    
                    with open(pdf_path, "rb") as pdf_file:
                        st.download_button(
                            label="⬇️ Baixar Relatório PDF",
                            data=pdf_file.read(),
                            file_name=report_name,
                            mime="application/pdf"
                        )
                    
                    st.success("✅ Relatório gerado com sucesso!")
                except Exception as e:
                    logger.error(f"Erro ao gerar PDF: {e}", exc_info=True)
                    st.error(f"❌ Erro ao gerar relatório: {str(e)}")

# Etapa 1: Carregamento e Validação
                status_text.text("📥 Carregando arquivo...")
                progress_bar.progress(10)
                
                file_bytes = uploaded_file.read()
                processor = get_fasta_processor(min_length, max_length)
                df_input = processor.load_peptides_from_bytes(file_bytes, uploaded_file.name)
                total_peptides = len(df_input)
                
                # --- NOVA MENSAGEM DE FEEDBACK ---
                st.info(f"📄 Arquivo lido: {total_peptides} sequências brutas encontradas.")
                
                if total_peptides == 0:
                    st.error("⚠️ Nenhuma sequência foi encontrada no arquivo. Verifique a formatação.")
                    st.stop()
                # ----------------------------------

                status_text.text(f"✅ {total_peptides} peptídeos carregados. Validando...")
                progress_bar.progress(30)
                
                df_valid = processor.validate_peptides(df_input)
                valid_count = len(df_valid)
                
                # --- NOVA MENSAGEM DE FEEDBACK ---
                if valid_count < total_peptides:
                    st.warning(f"🧹 Filtragem: {total_peptides - valid_count} sequências removidas (inválidas ou tamanho incorreto).")
                st.success(f"🧬 Processando {valid_count} sequências válidas.")
                
                if valid_count == 0:
                    st.error("⚠️ Nenhuma sequência válida restou após os filtros de tamanho/caracteres.")
                    st.stop()
           #-----------------------------
                # ... (código anterior de carregamento) ...
                
                status_text.text(f"✅ {total_peptides} peptídeos carregados. Validando...")
                progress_bar.progress(30)
                
                df_valid = processor.validate_peptides(df_input)
                valid_count = len(df_valid)
                
                # --- BLOCO DE SEGURANÇA NOVO ---
                if valid_count == 0:
                    st.error(f"❌ Erro de Validação: Das {total_peptides} sequências carregadas, 0 restaram.")
                    st.warning("Dica: Verifique se o 'Tamanho mínimo' e 'Máximo' nos filtros (barra lateral) condizem com seus dados. Se você carregou proteínas inteiras, aumente o tamanho máximo para 1000 ou mais.")
                    st.stop() # PARA AQUI E NÃO TENTA CONTINUAR
                # -------------------------------
                
                st.success(f"🧬 {valid_count} sequências válidas mantidas.")
                
                status_text.text(f"✅ {valid_count} peptídeos válidos. Calculando propriedades...")
                progress_bar.progress(50)
                
                # Agora é seguro chamar, pois garantimos que não está vazio
                df_result = processor.add_physchem_properties(df_valid)
                # ----------------------------------
                
                status_text.text(f"✅ {valid_count} peptídeos válidos. Calculando propriedades...")
                progress_bar.progress(50)
                
                # Aqui chama a função que estava dando erro antes
                df_result = processor.add_physchem_properties(df_valid)
                
if __name__ == "__main__":
    main()
