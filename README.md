# 🧬 Pipeline de Análise de Peptídeos

Aplicação web profissional para análise de epítopos e predição de apresentação MHC-I/II usando MHCflurry.

[![Streamlit](https://img.shields.io/badge/Streamlit-FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white)](https://streamlit.io/)
[![Python](https://img.shields.io/badge/Python-3.8+-blue?style=for-the-badge&logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)](LICENSE)

## 📋 Índice

- [Características](#-características)
- [Instalação](#-instalação)
- [Uso](#-uso)
- [Estrutura do Projeto](#-estrutura-do-projeto)
- [Citação](#-citação)
- [Contribuição](#-contribuição)
- [Licença](#-licença)

## ✨ Características

- ✅ **Interface Web Intuitiva**: Aplicação Streamlit moderna e responsiva
- ✅ **Predições MHC-I/II**: Suporte completo para ambos os tipos usando MHCflurry
- ✅ **Download Automático**: Modelos MHCflurry baixados automaticamente se não presentes
- ✅ **Análise de Conservação**: Cálculo de conservação posicional e em relação a referências
- ✅ **Propriedades Físico-Químicas**: Peso molecular, pI, GRAVY, hidrofobicidade Kyte-Doolittle
- ✅ **Enriquecimento com APIs**: Integração com IEDB e UniProt (opcional)
- ✅ **Relatórios PDF Profissionais**: Geração automática de dossiês com gráficos
- ✅ **Processamento Paralelo**: ThreadPoolExecutor para otimização de performance
- ✅ **Logging Profissional**: Sistema de logs completo para debugging
- ✅ **Código Modular**: Arquitetura em classes para fácil manutenção

## 🚀 Instalação

### Pré-requisitos

- Python 3.8 ou superior
- pip (gerenciador de pacotes Python)

### Instalação Local

1. **Clone o repositório**:
```bash
git clone https://github.com/SEU_USUARIO/SEU_REPO.git
cd SEU_REPO
```

2. **Instale as dependências**:
```bash
pip install -r requirements.txt
```

3. **Baixe os modelos MHCflurry** (primeira vez):
```bash
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads fetch models_class2_presentation
```

**Nota**: Os modelos serão baixados automaticamente na primeira execução se não estiverem presentes.

4. **Execute a aplicação**:
```bash
streamlit run app.py
```

A aplicação abrirá automaticamente em `http://localhost:8501`

### Deploy no Streamlit Cloud

1. **Faça push do código para GitHub**
2. **Acesse [share.streamlit.io](https://share.streamlit.io)**
3. **Faça login com GitHub**
4. **Clique em "New app"**
5. **Configure**:
   - Repository: seu repositório
   - Branch: `main`
   - Main file: `app.py`
6. **Clique em "Deploy"**

O Streamlit Cloud detecta automaticamente:
- `requirements.txt` - Instala dependências Python
- `packages.txt` - Instala pacotes do sistema (se necessário)
- `.streamlit/config.toml` - Configurações da aplicação

## 📖 Uso

### Interface Web

1. **Carregue seu arquivo**: FASTA, TSV, CSV ou Excel na sidebar
2. **Configure parâmetros**:
   - Tipo de MHC (I, II ou Ambos)
   - Alelo HLA
   - Tamanho do k-mer (8-14 resíduos)
   - Opções avançadas (APIs, conservação)
3. **Clique em "Processar Análise"**
4. **Visualize resultados** nas tabs:
   - **Resumo**: Métricas principais e gráficos
   - **Sequências**: Tabela completa com filtros
   - **Análise de Epítopos**: Top candidatos detalhados
   - **Relatório Final**: Geração de PDF

### Uso Programático

```python
from src import FastaProcessor, MHCAnalyzer, PDFGenerator

# Processa arquivo FASTA
processor = FastaProcessor(min_length=8, max_length=14)
df = processor.load_peptides_from_bytes(file_bytes, "peptides.fasta")
df_valid = processor.validate_peptides(df)
df_props = processor.add_physchem_properties(df_valid)

# Predições MHC
mhc_analyzer = MHCAnalyzer()
df_mhc = mhc_analyzer.predict_class1(df_props, ["HLA-A*02:01"])
df_final = mhc_analyzer.calculate_final_scores(df_mhc)

# Gera relatório PDF
pdf_gen = PDFGenerator()
pdf_gen.generate_report(df_final, "relatorio.pdf")
```

## 📁 Estrutura do Projeto

```
.
├── app.py                    # Aplicação Streamlit principal
├── src/                      # Módulos do projeto
│   ├── __init__.py
│   ├── data_processor.py    # Processamento de FASTA
│   ├── mhc_analyzer.py      # Predições MHC-I/II
│   ├── api_manager.py       # Gerenciamento de APIs
│   ├── pdf_generator.py     # Geração de PDFs
│   └── conservation.py      # Análise de conservação
├── tests/                    # Testes unitários
│   ├── __init__.py
│   └── test_fasta.py        # Testes de leitura FASTA
├── data/                     # Arquivos de exemplo
│   └── .gitkeep
├── .streamlit/              # Configurações Streamlit
│   └── config.toml
├── requirements.txt         # Dependências Python
├── packages.txt            # Pacotes do sistema (Streamlit Cloud)
└── README.md               # Este arquivo
```

## 📚 Citação

Se você usar este pipeline em sua pesquisa, por favor cite:

### Ferramentas Utilizadas

**MHCflurry**:
```
O'Donnell, T. J., Rubinsteyn, A., Bonsack, M., Riemer, A. B., Laserson, U., & Hammerbacher, J. (2018).
MHCflurry: Open-Source Class I MHC Binding Affinity Prediction.
Cell Systems, 7(1), 129-132.e4.
```

**IEDB (Immune Epitope Database)**:
```
Vita, R., Mahajan, S., Overton, J. A., Dhanda, S. K., Martini, S., Cantrell, J. R., ... & Peters, B. (2019).
The Immune Epitope Database (IEDB): 2018 update.
Nucleic Acids Research, 47(D1), D339-D343.
```

**UniProt**:
```
UniProt Consortium. (2023).
UniProt: the Universal Protein Knowledgebase in 2023.
Nucleic Acids Research, 51(D1), D523-D531.
```

**Biopython**:
```
Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & de Hoon, M. J. (2009).
Biopython: freely available Python tools for computational molecular biology and bioinformatics.
Bioinformatics, 25(11), 1422-1423.
```

### Este Projeto

```
Pipeline de Análise de Peptídeos v2.0.0
Engenheiro de Bioinformática Sênior
2026
```

## 🧪 Testes

Execute os testes unitários:

```bash
python -m pytest tests/
```

Ou usando unittest:

```bash
python -m unittest tests.test_fasta
```

## 🔧 Desenvolvimento

### Estrutura de Classes

- **FastaProcessor**: Processamento e validação de sequências
- **MHCAnalyzer**: Predições MHC com download automático de modelos
- **APIManager**: Gerenciamento de APIs externas com processamento paralelo
- **PDFGenerator**: Geração de relatórios PDF profissionais
- **ConservationAnalyzer**: Análise de conservação de sequências

### Logging

O projeto usa logging profissional. Configure o nível:

```python
import logging
logging.basicConfig(level=logging.DEBUG)  # Para debug detalhado
```

### Contribuindo

1. Fork o projeto
2. Crie uma branch para sua feature (`git checkout -b feature/AmazingFeature`)
3. Commit suas mudanças (`git commit -m 'Add some AmazingFeature'`)
4. Push para a branch (`git push origin feature/AmazingFeature`)
5. Abra um Pull Request

## 📝 Licença

Este projeto está sob a licença MIT. Veja o arquivo `LICENSE` para mais detalhes.

## 🙏 Agradecimentos

- MHCflurry por fornecer modelos de predição de alta qualidade
- IEDB e UniProt por disponibilizarem APIs públicas
- Comunidade Streamlit por uma excelente plataforma de desenvolvimento web

## 📧 Contato

Para questões, problemas ou sugestões, abra uma issue no GitHub.

---

**Desenvolvido com ❤️ por Engenheiros de Bioinformática**
