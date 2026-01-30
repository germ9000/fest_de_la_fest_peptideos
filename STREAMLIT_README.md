# 🧬 Aplicação Streamlit - Pipeline de Análise de Peptídeos

Aplicação web completa e interativa para análise de epítopos e geração de dossiês em PDF.

## 🚀 Instalação Rápida

### 1. Instale as dependências:
```bash
pip install -r requirements.txt
```

### 2. Baixe os modelos MHCflurry (primeira vez):
```bash
mhcflurry-downloads fetch models_class1_presentation
```

### 3. Execute a aplicação:
```bash
streamlit run app.py
```

A aplicação abrirá automaticamente no navegador em `http://localhost:8501`

## 📋 Funcionalidades

### ✅ Interface Intuitiva
- **Sidebar de Configurações**: Parâmetros de análise facilmente ajustáveis
- **Upload de Arquivos**: Suporte a FASTA, TSV, CSV e Excel
- **Tabs Organizadas**: Resultados divididos em seções lógicas

### 🔬 Análises Disponíveis
1. **Validação Automática**: Filtra peptídeos inválidos
2. **Propriedades Físico-Químicas**:
   - Peso molecular
   - Ponto isoelétrico (pI)
   - Índice GRAVY
   - **Hidrofobicidade Kyte-Doolittle** (novo!)
   - Índice de instabilidade
   - Índice alifático

3. **Predições MHC-I**:
   - Afinidade de ligação (MHCflurry)
   - Score de apresentação
   - Percentil de rank

4. **Enriquecimento com APIs** (opcional):
   - **IEDB**: Imunogenicidade e afinidade
   - **UniProt**: Anotações funcionais

### 📊 Visualizações Interativas
- Gráficos Plotly interativos
- Tabelas com cores dinâmicas (afinidade < 50nM destacada)
- Distribuições estatísticas
- Scatter plots de propriedades

### 📄 Relatórios PDF Profissionais
- Resumo executivo
- Gráficos de afinidade
- Tabelas formatadas
- Metadados da análise
- Estatísticas descritivas

## 🎯 Como Usar

### Passo 1: Carregue seu arquivo
1. Na barra lateral, clique em "Browse files"
2. Selecione arquivo FASTA, TSV, CSV ou Excel
3. O arquivo será carregado automaticamente

### Passo 2: Configure os parâmetros
- **Alelo HLA**: Selecione o alelo para predição (ex: HLA-A*02:01)
- **Tamanho do k-mer**: Ajuste min/max (padrão: 8-14 para MHC-I)
- **APIs Externas**: Ative/desative enriquecimento (pode ser mais lento)
- **Threads**: Ajuste paralelismo (padrão: 5)

### Passo 3: Processe a análise
1. Clique em "🚀 Processar Análise"
2. Aguarde o processamento (barra de progresso)
3. Resultados aparecerão automaticamente nas tabs

### Passo 4: Visualize os resultados

#### Tab "📊 Resumo"
- Métricas principais
- Gráficos de afinidade
- Estatísticas descritivas

#### Tab "🧬 Sequências"
- Tabela completa de peptídeos
- Filtros por afinidade
- Busca de sequências
- Cores dinâmicas (verde = forte, amarelo = moderado, vermelho = fraco)

#### Tab "🎯 Análise de Epítopos"
- Top 10 candidatos
- Visualizações detalhadas
- Distribuições

#### Tab "📄 Relatório Final"
- Geração de PDF profissional
- Exportação de dados (Excel, CSV, TSV)

## ⚙️ Configurações Avançadas

### Performance
- **Threads Paralelas**: Aumente para processar mais rápido (cuidado com rate limiting)
- **Cache**: APIs são cacheadas automaticamente (TTL: 1 hora)

### APIs Externas
- **IEDB**: Predições de imunogenicidade e afinidade
- **UniProt**: Busca informações sobre proteínas (mais lento)

### Exportação
- **PDF**: Relatório completo com gráficos
- **Excel**: Dados completos em planilha
- **CSV/TSV**: Formato texto para análise externa

## 🐛 Troubleshooting

### Erro: "Modelos MHCflurry não encontrados"
```bash
mhcflurry-downloads fetch models_class1_presentation
```

### Erro: "Módulo não encontrado"
```bash
pip install -r requirements.txt
```

### APIs retornando None
- Verifique conexão com internet
- APIs podem estar temporariamente indisponíveis
- Aumente `request_delay` na configuração

### Aplicação lenta
- Reduza número de threads
- Desative APIs externas se não necessário
- Processe menos peptídeos por vez

## 📊 Exemplo de Saída

### Métricas Principais
- Total de peptídeos analisados
- Peptídeos válidos após validação
- Melhor afinidade encontrada (nM)
- Taxa de validação (%)

### Gráficos
- Top 20 por afinidade (barra horizontal)
- Propriedades físico-químicas (scatter plot)
- Distribuições (histogramas)

### Tabelas
- Colunas: Peptídeo, Afinidade (nM), Score MHC, Propriedades
- Cores: Verde (<50nM), Amarelo (50-500nM), Vermelho (>500nM)

## 🔧 Arquitetura

```
app.py                 # Interface Streamlit principal
├── data_handler.py   # Manipulação de dados
├── api_client.py     # Clientes de API (IEDB, UniProt)
└── report_gen.py     # Geração de PDFs
```

## 📝 Notas Técnicas

- **Cache**: Streamlit cache automático para APIs (evita chamadas repetidas)
- **Threading**: Processamento paralelo com ThreadPoolExecutor
- **Rate Limiting**: Delays automáticos entre requisições
- **Validação**: Filtros biológicos rigorosos (aminoácidos canônicos, tamanho)

## 🎓 Referências

- **MHCflurry**: Predição de apresentação MHC-I
- **IEDB**: Immune Epitope Database
- **UniProt**: Universal Protein Resource
- **Kyte-Doolittle**: Escala de hidrofobicidade (1982)

## 📧 Suporte

Para problemas ou dúvidas, consulte a documentação dos módulos individuais ou abra uma issue.

---

**Desenvolvido por Engenheiro de Bioinformática Sênior** 🧬
