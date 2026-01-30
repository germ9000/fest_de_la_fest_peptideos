"""
Aplicação principal Streamlit para análise de epítopos imunológicos.
"""
import streamlit as st
import pandas as pd
import logging
from typing import List, Optional, Dict, Any
import tempfile
import os

from src.data_handler import FastaProcessor, ProteinAnalyzer
from src.mhc_analyzer import MHCAnalyzer
from src.api_client import IEDBClient, UniProtClient
from src.report_gen import PDFReportGenerator

# Configuração de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuração da página
st.set_page_config(
    page_title="ImmunoEpitope Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@st.cache_resource
def load_mhc_analyzer() -> MHCAnalyzer:
    """Carrega o analisador MHCflurry uma vez na memória."""
    logger.info("Carregando MHCflurry analyzer...")
    return MHCAnalyzer()

class ImmunoEpitopeApp:
    """Classe principal da aplicação Streamlit."""
    
    def __init__(self):
        self.initialize_session_state()
        self.setup_ui()
        
    def initialize_session_state(self):
        """Inicializa o estado da sessão."""
        if 'analysis_results' not in st.session_state:
            st.session_state.analysis_results = None
        if 'fasta_data' not in st.session_state:
            st.session_state.fasta_data = None
        if 'mhc_predictions' not in st.session_state:
            st.session_state.mhc_predictions = {'mhc_i': None, 'mhc_ii': None}
        if 'uniprot_data' not in st.session_state:
            st.session_state.uniprot_data = None
    
    def setup_ui(self):
        """Configura a interface do usuário."""
        st.title("🧬 ImmunoEpitope Pipeline")
        st.markdown("""
        Aplicação para análise de epítopos imunológicos a partir de sequências proteicas.
        Faça upload de um arquivo FASTA para iniciar a análise.
        """)
        
        self.render_sidebar()
        self.render_main_content()
    
    def render_sidebar(self):
        """Renderiza a barra lateral com controles."""
        with st.sidebar:
            st.header("⚙️ Configurações de Análise")
            
            # Upload do arquivo
            uploaded_file = st.file_uploader(
                "📁 Upload de arquivo FASTA",
                type=['fasta', 'fa', 'txt'],
                help="Faça upload de um arquivo FASTA contendo sequências proteicas"
            )
            
            if uploaded_file:
                self.process_uploaded_file(uploaded_file)
            
            # Parâmetros MHC-I
            st.subheader("🎯 Parâmetros MHC-I")
            self.mhci_alleles = st.text_input(
                "Alelos MHC-I (separados por vírgula)",
                value="HLA-A*02:01,HLA-A*03:01,HLA-B*07:02",
                help="Ex: HLA-A*02:01,HLA-B*07:02"
            )
            self.mhci_lengths = st.multiselect(
                "Comprimentos de peptídeos",
                options=[8, 9, 10, 11, 12],
                default=[9],
                help="Selecione os comprimentos para predição"
            )
            
            # Parâmetros MHC-II
            st.subheader("🎯 Parâmetros MHC-II")
            self.mhcii_alleles = st.text_input(
                "Alelos MHC-II (separados por vírgula)",
                value="HLA-DRB1*01:01,HLA-DRB1*04:01",
                help="Ex: HLA-DRB1*01:01,HLA-DRB1*04:01"
            )
            
            # Thresholds
            st.subheader("📊 Thresholds")
            self.ic50_threshold = st.number_input(
                "IC50 Threshold (nM)",
                min_value=0.0,
                max_value=10000.0,
                value=50.0,
                step=1.0,
                help="Peptídeos com IC50 <= valor são considerados ligantes"
            )
            self.percentile_threshold = st.number_input(
                "Percentile Threshold",
                min_value=0.0,
                max_value=100.0,
                value=2.0,
                step=0.1,
                help="Peptídeos com percentile <= valor são considerados ligantes"
            )
            
            # Botão de análise
            if st.button("🚀 Executar Análise", type="primary", use_container_width=True):
                if st.session_state.fasta_data:
                    self.run_analysis()
                else:
                    st.error("Por favor, faça upload de um arquivo FASTA primeiro.")
    
    def process_uploaded_file(self, uploaded_file):
        """Processa o arquivo FASTA carregado."""
        try:
            # Salva temporariamente o arquivo
            with tempfile.NamedTemporaryFile(delete=False, suffix='.fasta') as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_path = tmp_file.name
            
            # Processa o FASTA
            processor = FastaProcessor(tmp_path)
            fasta_data = processor.process()
            
            if fasta_data:
                st.session_state.fasta_data = fasta_data
                st.success(f"✅ {len(fasta_data)} sequência(s) carregada(s)")
                
                # Mostra prévia
                with st.expander("👁️ Prévia das sequências"):
                    preview_df = pd.DataFrame([
                        {
                            'ID': seq.id,
                            'Descrição': seq.description[:50] + "..." if len(seq.description) > 50 else seq.description,
                            'Tamanho': len(seq.seq)
                        }
                        for seq in fasta_data[:5]
                    ])
                    st.dataframe(preview_df, use_container_width=True)
            
            # Limpa arquivo temporário
            os.unlink(tmp_path)
            
        except Exception as e:
            st.error(f"Erro ao processar arquivo: {str(e)}")
            logger.error(f"Erro no processamento de FASTA: {e}", exc_info=True)
    
    def run_analysis(self):
        """Executa a análise completa."""
        with st.spinner("🔬 Executando análise... Isso pode levar alguns minutos."):
            try:
                # 1. Análise físico-química
                protein_analyzer = ProteinAnalyzer()
                physchem_results = protein_analyzer.analyze_proteins(
                    st.session_state.fasta_data
                )
                
                # 2. Predição MHC-I
                mhc_analyzer = load_mhc_analyzer()
                mhci_predictions = mhc_analyzer.predict_mhc_i(
                    sequences=st.session_state.fasta_data,
                    alleles=self.mhci_alleles.split(','),
                    lengths=self.mhci_lengths
                )
                
                # 3. Predição MHC-II (IEDB)
                iedb_client = IEDBClient()
                mhcii_predictions = iedb_client.predict_mhc_ii(
                    sequences=st.session_state.fasta_data,
                    alleles=self.mhcii_alleles.split(',')
                )
                
                # 4. Anotações UniProt
                uniprot_client = UniProtClient()
                uniprot_data = uniprot_client.fetch_annotations(
                    [seq.id for seq in st.session_state.fasta_data]
                )
                
                # Armazena resultados
                st.session_state.analysis_results = {
                    'physchem': physchem_results,
                    'mhci': mhci_predictions,
                    'mhcii': mhcii_predictions,
                    'uniprot': uniprot_data,
                    'parameters': {
                        'mhci_alleles': self.mhci_alleles,
                        'mhcii_alleles': self.mhcii_alleles,
                        'ic50_threshold': self.ic50_threshold,
                        'percentile_threshold': self.percentile_threshold
                    }
                }
                
                st.success("✅ Análise concluída!")
                
            except Exception as e:
                st.error(f"Erro na análise: {str(e)}")
                logger.error(f"Erro na execução da análise: {e}", exc_info=True)
    
    def render_main_content(self):
        """Renderiza o conteúdo principal com abas."""
        if not st.session_state.analysis_results:
            self.render_welcome()
            return
        
        tab1, tab2, tab3, tab4 = st.tabs([
            "📊 Sequências",
            "🎯 Análise de Epítopos",
            "📈 Visualizações",
            "📄 Relatório Final"
        ])
        
        with tab1:
            self.render_sequences_tab()
        
        with tab2:
            self.render_epitope_tab()
        
        with tab3:
            self.render_visualizations_tab()
        
        with tab4:
            self.render_report_tab()
    
    def render_welcome(self):
        """Renderiza tela de boas-vindas."""
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            st.image("https://cdn-icons-png.flaticon.com/512/1998/1998673.png", width=200)
            st.markdown("""
            ## Bem-vindo ao ImmunoEpitope Pipeline!
            
            1. 📁 **Faça upload** de um arquivo FASTA na barra lateral
            2. ⚙️ **Configure** os parâmetros de análise
            3. 🚀 **Clique em 'Executar Análise'**
            4. 📊 **Explore** os resultados nas abas
            
            ### 📚 Formato FASTA Esperado:
            ```fasta
            >sp|P04637|P53_HUMAN Cellular tumor antigen p53
            MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
            >sp|P62993|GRB2_HUMAN Growth factor receptor-bound protein 2
            MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
            ```
            """)
    
    def render_sequences_tab(self):
        """Renderiza aba de sequências."""
        st.header("📊 Análise das Sequências")
        
        if st.session_state.fasta_data:
            # Tabela de sequências
            seq_data = []
            for seq in st.session_state.fasta_data:
                seq_data.append({
                    'ID': seq.id,
                    'Descrição': seq.description,
                    'Tamanho': len(seq.seq),
                    'GC %': self.calculate_gc_content(str(seq.seq)),
                    'Peso Molecular': f"{len(seq.seq) * 110:.0f} Da"
                })
            
            df = pd.DataFrame(seq_data)
            st.dataframe(df, use_container_width=True)
            
            # Propriedades físico-químicas
            if 'physchem' in st.session_state.analysis_results:
                st.subheader("🧪 Propriedades Físico-Químicas")
                physchem_df = pd.DataFrame(st.session_state.analysis_results['physchem'])
                st.dataframe(physchem_df.style.background_gradient(subset=['pI', 'GRAVY']), use_container_width=True)
    
    def render_epitope_tab(self):
        """Renderiza aba de análise de epítopos."""
        st.header("🎯 Predição de Epítopos")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("MHC-I - Melhores Epítopos")
            if st.session_state.analysis_results['mhci'] is not None:
                mhci_df = pd.DataFrame(st.session_state.analysis_results['mhci'])
                
                # Aplica threshold
                filtered_df = mhci_df[
                    (mhci_df['ic50'] <= self.ic50_threshold) & 
                    (mhci_df['percentile'] <= self.percentile_threshold)
                ]
                
                if not filtered_df.empty:
                    # Estiliza a tabela
                    def highlight_binders(row):
                        if row['ic50'] < 50:
                            return ['background-color: #90EE90'] * len(row)
                        elif row['ic50'] < 500:
                            return ['background-color: #FFFACD'] * len(row)
                        return [''] * len(row)
                    
                    st.dataframe(
                        filtered_df.sort_values('ic50').head(20).style.apply(
                            highlight_binders, axis=1
                        ),
                        use_container_width=True
                    )
                    
                    # Métricas
                    col1a, col2a, col3a = st.columns(3)
                    with col1a:
                        st.metric("Total de Ligantes", len(filtered_df))
                    with col2a:
                        st.metric("IC50 Médio", f"{filtered_df['ic50'].mean():.1f} nM")
                    with col3a:
                        st.metric("Melhor IC50", f"{filtered_df['ic50'].min():.1f} nM")
                else:
                    st.warning("Nenhum epítopo forte encontrado com os thresholds atuais.")
        
        with col2:
            st.subheader("MHC-II - Melhores Epítopos")
            if st.session_state.analysis_results['mhcii'] is not None:
                mhcii_df = pd.DataFrame(st.session_state.analysis_results['mhcii'])
                
                if not mhcii_df.empty:
                    st.dataframe(
                        mhcii_df.sort_values('score').head(20),
                        use_container_width=True
                    )
    
    def render_visualizations_tab(self):
        """Renderiza aba de visualizações."""
        st.header("📈 Visualizações")
        
        if st.session_state.analysis_results['mhci'] is not None:
            import plotly.express as px
            import plotly.graph_objects as go
            
            mhci_df = pd.DataFrame(st.session_state.analysis_results['mhci'])
            
            # Gráfico 1: Distribuição de IC50
            fig1 = px.histogram(
                mhci_df, 
                x='ic50',
                nbins=50,
                title='Distribuição de IC50',
                labels={'ic50': 'IC50 (nM)'}
            )
            fig1.add_vline(
                x=self.ic50_threshold,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Threshold: {self.ic50_threshold}nM"
            )
            st.plotly_chart(fig1, use_container_width=True)
            
            # Gráfico 2: IC50 vs Percentile
            fig2 = px.scatter(
                mhci_df,
                x='ic50',
                y='percentile',
                color='allele',
                hover_data=['peptide'],
                title='IC50 vs Percentile por Alelo'
            )
            st.plotly_chart(fig2, use_container_width=True)
    
    def render_report_tab(self):
        """Renderiza aba de relatório."""
        st.header("📄 Relatório Final")
        
        if st.session_state.analysis_results:
            # Botão para gerar PDF
            if st.button("📥 Gerar Relatório PDF", type="primary"):
                with st.spinner("Gerando relatório PDF..."):
                    try:
                        report_gen = PDFReportGenerator()
                        
                        # Gera o relatório
                        pdf_path = report_gen.generate_report(
                            fasta_data=st.session_state.fasta_data,
                            analysis_results=st.session_state.analysis_results,
                            parameters=st.session_state.analysis_results['parameters']
                        )
                        
                        # Disponibiliza para download
                        with open(pdf_path, "rb") as pdf_file:
                            pdf_bytes = pdf_file.read()
                        
                        st.download_button(
                            label="⬇️ Baixar Relatório PDF",
                            data=pdf_bytes,
                            file_name="immunoepitope_report.pdf",
                            mime="application/pdf"
                        )
                        
                        st.success("Relatório gerado com sucesso!")
                        
                    except Exception as e:
                        st.error(f"Erro ao gerar relatório: {str(e)}")
            
            # Prévia do relatório
            with st.expander("👁️ Prévia do Relatório"):
                st.markdown(self.generate_report_preview())
    
    def calculate_gc_content(self, sequence: str) -> float:
        """Calcula conteúdo de GC."""
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0
    
    def generate_report_preview(self) -> str:
        """Gera prévia do relatório em markdown."""
        results = st.session_state.analysis_results
        
        preview = f"""
        # Relatório de Análise ImmunoEpitope
        
        ## 📅 Data da Análise
        {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
        
        ## ⚙️ Parâmetros Utilizados
        - **Alelos MHC-I**: {results['parameters']['mhci_alleles']}
        - **Alelos MHC-II**: {results['parameters']['mhcii_alleles']}
        - **Threshold IC50**: {results['parameters']['ic50_threshold']} nM
        - **Threshold Percentile**: {results['parameters']['percentile_threshold']} %
        
        ## 📊 Estatísticas
        - **Sequências analisadas**: {len(st.session_state.fasta_data)}
        - **Epítopos MHC-I preditos**: {len(results['mhci']) if results['mhci'] else 0}
        - **Epítopos MHC-II preditos**: {len(results['mhcii']) if results['mhcii'] else 0}
        
        ## 🎯 Top 5 Epítopos MHC-I
        """
        
        if results['mhci']:
            mhci_df = pd.DataFrame(results['mhci'])
            top_5 = mhci_df.sort_values('ic50').head(5)
            
            preview += "\n| Peptídeo | Alelo | IC50 (nM) | Percentile |\n"
            preview += "|----------|-------|-----------|------------|\n"
            for _, row in top_5.iterrows():
                preview += f"| {row['peptide']} | {row['allele']} | {row['ic50']:.2f} | {row['percentile']:.2f} |\n"
        
        return preview

def main():
    """Função principal."""
    try:
        app = ImmunoEpitopeApp()
    except Exception as e:
        st.error(f"Erro crítico na aplicação: {str(e)}")
        logger.error(f"Erro crítico: {e}", exc_info=True)

if __name__ == "__main__":
    main()
