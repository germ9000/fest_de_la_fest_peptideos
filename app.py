import streamlit as st
import pandas as pd
from src.data_handler import FastaProcessor
from src.mhc_analyzer import MHCAnalyzer
from src.api_client import APIClient
from src.report_gen import PDFGenerator
import logging
import time

# Configuração básica de logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configuração da página
st.set_page_config(
    page_title="ImmunoEpitope Pipeline",
    page_icon="🧬",
    layout="wide"
)

# Título da aplicação
st.title("🧬 ImmunoEpitope Pipeline")
st.markdown("""
Esta aplicação processa arquivos FASTA para predição de epítopos imunológicos, 
utilizando **MHCflurry** para predição de ligação MHC-I e **IEDB API** para predição MHC-II.
""")

# Sidebar para configurações
st.sidebar.header("Configurações de Análise")

# Upload do arquivo FASTA
uploaded_file = st.sidebar.file_uploader("Faça upload do arquivo FASTA", type=['fasta', 'fa', 'txt'])

# Parâmetros de predição MHC-I
st.sidebar.subheader("Parâmetros de Predição MHC-I")
alleles_mhci = st.sidebar.text_input("Alelos MHC-I (separados por vírgula)", "HLA-A*02:01")
alleles_mhci = [allele.strip() for allele in alleles_mhci.split(",")]

# Parâmetros de predição MHC-II
st.sidebar.subheader("Parâmetros de Predição MHC-II")
alleles_mhcii = st.sidebar.text_input("Alelos MHC-II (separados por vírgula)", "HLA-DRB1*01:01")
alleles_mhcii = [allele.strip() for allele in alleles_mhcii.split(",")]

# Tamanhos dos peptídeos (k-mer)
peptide_lengths_mhci = st.sidebar.multiselect("Comprimentos dos peptídeos para MHC-I", [8, 9, 10, 11], default=[8, 9])
peptide_lengths_mhcii = st.sidebar.multiselect("Comprimentos dos peptídeos para MHC-II", [13, 14, 15], default=[15])

# Thresholds
st.sidebar.subheader("Thresholds de Afinidade")
ic50_threshold = st.sidebar.number_input("IC50 Threshold (nM) para MHC-I", min_value=0.0, value=50.0, step=1.0)
percentile_threshold = st.sidebar.number_input("Percentil Threshold para MHC-I", min_value=0.0, value=2.0, step=0.1)

# Botão para iniciar a análise
run_analysis = st.sidebar.button("Iniciar Análise")

# Inicialização de variáveis de estado
if 'fasta_processor' not in st.session_state:
    st.session_state.fasta_processor = None
if 'mhc_analyzer' not in st.session_state:
    st.session_state.mhc_analyzer = None
if 'api_client' not in st.session_state:
    st.session_state.api_client = None
if 'results_mhci' not in st.session_state:
    st.session_state.results_mhci = None
if 'results_mhcii' not in st.session_state:
    st.session_state.results_mhcii = None
if 'uniprot_data' not in st.session_state:
    st.session_state.uniprot_data = None

# Cache para o modelo MHCflurry
@st.cache_resource
def load_mhcflurry_model():
    logger.info("Carregando modelo MHCflurry...")
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    return predictor

# Função principal de análise
def analyze_fasta(uploaded_file, alleles_mhci, alleles_mhcii, peptide_lengths_mhci, peptide_lengths_mhcii):
    # Processar o arquivo FASTA
    fasta_processor = FastaProcessor(uploaded_file)
    sequences = fasta_processor.parse_fasta()
    st.session_state.fasta_processor = fasta_processor

    # Inicializar o analisador MHC com o modelo carregado
    predictor = load_mhcflurry_model()
    mhc_analyzer = MHCAnalyzer(predictor)
    st.session_state.mhc_analyzer = mhc_analyzer

    # Inicializar o cliente de API
    api_client = APIClient()
    st.session_state.api_client = api_client

    # Abas para resultados
    tab1, tab2, tab3, tab4 = st.tabs(["Sequências", "Análise de Epítopos MHC-I", "Análise de Epítopos MHC-II", "Relatório Final"])

    with tab1:
        st.header("Sequências Proteicas Carregadas")
        seq_df = pd.DataFrame([(header, str(seq)) for header, seq in sequences], columns=["ID", "Sequência"])
        st.dataframe(seq_df)

    # Predições MHC-I
    with tab2:
        st.header("Predição de Epítopos MHC-I")
        if sequences:
            with st.spinner("Realizando predições MHC-I. Isso pode levar alguns instantes..."):
                results_mhci = mhc_analyzer.predict_mhci_epitopes(sequences, alleles_mhci, peptide_lengths_mhci)
                st.session_state.results_mhci = results_mhci
                if results_mhci:
                    # Combinar todos os resultados em um DataFrame
                    df_list = []
                    for result in results_mhci:
                        df_list.append(result['predictions'])
                    combined_df = pd.concat(df_list, ignore_index=True)
                    # Filtrar por threshold
                    filtered_df = combined_df[(combined_df['ic50'] <= ic50_threshold) & (combined_df['percentile'] <= percentile_threshold)]
                    st.subheader(f"Epítopos Preditos (IC50 <= {ic50_threshold} nM e Percentil <= {percentile_threshold})")
                    st.dataframe(filtered_df)
                    # Opção de download
                    csv = filtered_df.to_csv(index=False)
                    st.download_button(
                        label="Baixar resultados MHC-I como CSV",
                        data=csv,
                        file_name="mhci_predictions.csv",
                        mime="text/csv",
                    )
                else:
                    st.warning("Nenhuma predição MHC-I foi retornada.")

    # Predições MHC-II via IEDB API
    with tab3:
        st.header("Predição de Epítopos MHC-II (via IEDB API)")
        if sequences:
            with st.spinner("Realizando predições MHC-II. Isso pode levar alguns minutos devido às chamadas de API..."):
                # Vamos limitar o número de sequências e peptídeos para não sobrecarregar a API
                # Por exemplo, pegar apenas a primeira sequência
                if sequences:
                    first_seq = sequences[0]
                    # Gerar peptídeos da primeira sequência
                    peptides = fasta_processor.generate_peptides(str(first_seq[1]), peptide_lengths_mhcii)
                    # Limitar a 50 peptídeos para demonstração
                    peptides = peptides[:50]
                    # Predizer para cada alelo MHC-II
                    all_results = []
                    for allele in alleles_mhcii:
                        results = api_client.predict_mhcii_iedb(peptides, allele, method='nn_align')
                        if results:
                            all_results.extend(results)
                    if all_results:
                        results_mhcii = pd.DataFrame(all_results)
                        st.session_state.results_mhcii = results_mhcii
                        st.dataframe(results_mhcii)
                        # Opção de download
                        csv = results_mhcii.to_csv(index=False)
                        st.download_button(
                            label="Baixar resultados MHC-II como CSV",
                            data=csv,
                            file_name="mhcii_predictions.csv",
                            mime="text/csv",
                        )
                    else:
                        st.warning("Nenhuma predição MHC-II foi retornada.")

    # Relatório Final
    with tab4:
        st.header("Relatório Final")
        if st.button("Gerar Relatório PDF"):
            if (st.session_state.results_mhci is not None) or (st.session_state.results_mhcii is not None):
                with st.spinner("Gerando relatório PDF..."):
                    # Coletar dados para o relatório
                    fasta_data = fasta_processor.get_fasta_data()
                    # Gerar o PDF
                    pdf_generator = PDFGenerator()
                    pdf_file = pdf_generator.generate_report(
                        fasta_data=fasta_data,
                        mhci_results=st.session_state.results_mhci,
                        mhcii_results=st.session_state.results_mhcii,
                        alleles_mhci=alleles_mhci,
                        alleles_mhcii=alleles_mhcii,
                        peptide_lengths_mhci=peptide_lengths_mhci,
                        peptide_lengths_mhcii=peptide_lengths_mhcii,
                        ic50_threshold=ic50_threshold,
                        percentile_threshold=percentile_threshold
                    )
                    # Disponibilizar o download
                    with open(pdf_file, "rb") as f:
                        st.download_button(
                            label="Baixar Relatório PDF",
                            data=f,
                            file_name="immunoepitope_report.pdf",
                            mime="application/pdf"
                        )
            else:
                st.error("Nenhum resultado disponível para gerar o relatório. Execute as análises primeiro.")

# Executar a análise quando o botão for pressionado
if run_analysis and uploaded_file is not None:
    analyze_fasta(uploaded_file, alleles_mhci, alleles_mhcii, peptide_lengths_mhci, peptide_lengths_mhcii)
elif run_analysis and uploaded_file is None:
    st.error("Por favor, faça upload de um arquivo FASTA.")
