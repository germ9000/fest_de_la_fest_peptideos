# 🚀 Guia Rápido de Início

## ⚡ Execução Local (3 Passos)

### 1️⃣ Instalar Dependências

**Windows:**
```powershell
# Opção 1: Script automático
.\setup.bat

# Opção 2: Manual
pip install -r requirements.txt
mhcflurry-downloads fetch models_class1_presentation
```

**Linux/Mac:**
```bash
# Opção 1: Script automático
chmod +x setup.sh
./setup.sh

# Opção 2: Manual
pip install -r requirements.txt
mhcflurry-downloads fetch models_class1_presentation
```

**Python:**
```bash
python setup.py
```

### 2️⃣ Testar Instalação (Opcional)

```bash
python test_installation.py
```

Deve mostrar: ✅ Todos os módulos instalados corretamente!

### 3️⃣ Executar Aplicação

```bash
streamlit run app.py
```

A aplicação abrirá automaticamente em: **http://localhost:8501**

---

## 📋 Como Usar

1. **Carregue arquivo**: FASTA, TSV, CSV ou Excel na sidebar
2. **Configure parâmetros**:
   - Tipo de MHC (I, II ou Ambos)
   - Alelo HLA
   - Tamanho do k-mer
3. **Clique em "Processar Análise"**
4. **Visualize resultados** nas tabs

---

## 🌐 Deploy no Streamlit Cloud

### Pré-requisitos
- Conta GitHub
- Repositório com código

### Passos

1. **Crie repositório no GitHub**
```bash
git init
git add .
git commit -m "Streamlit app"
git remote add origin https://github.com/SEU_USUARIO/SEU_REPO.git
git push -u origin main
```

2. **Deploy no Streamlit Cloud**
   - Acesse: https://share.streamlit.io
   - Login com GitHub
   - Clique "New app"
   - Configure:
     - Repository: seu repositório
     - Main file: `app.py`
   - Clique "Deploy"

3. **Aguarde** (primeira vez pode levar 5-10 minutos)

---

## 🆕 Novas Funcionalidades

### ✅ Implementadas

1. **Predições MHC-II**
   - Suporte completo para HLA-DR, DQ, DP
   - Análise combinada MHC-I + MHC-II

2. **Análise de Conservação**
   - Score de conservação posicional
   - Conservação em relação a sequências de referência
   - Entropia de Shannon

3. **Integração NetMHCpan**
   - Predições adicionais via API
   - Comparação com MHCflurry

---

## 📁 Estrutura de Arquivos

```
.
├── app.py                    # 🎯 Aplicação principal
├── data_handler.py          # 📊 Manipulação de dados
├── api_client.py            # 🌐 APIs (IEDB, UniProt)
├── report_gen.py           # 📄 Geração de PDFs
├── mhc_predictions.py      # 🧬 Predições MHC-I/II
├── conservation_analysis.py # 🔬 Análise de conservação
├── netmhcpan_client.py     # 🧪 Cliente NetMHCpan
├── requirements.txt         # 📦 Dependências
├── setup.py                # ⚙️ Setup automático
├── test_installation.py    # ✅ Teste de instalação
└── README_DEPLOY.md        # 📖 Guia de deploy
```

---

## ❓ Problemas?

Consulte: **TROUBLESHOOTING.md**

Problemas comuns:
- `ModuleNotFoundError` → `pip install -r requirements.txt`
- `MHCflurry models not found` → `mhcflurry-downloads fetch models_class1_presentation`
- `Port 8501 in use` → Feche outras instâncias do Streamlit

---

## 🎯 Exemplo de Uso

1. Abra `app.py` no Streamlit
2. Faça upload de arquivo FASTA com peptídeos
3. Selecione:
   - Tipo: MHC-I
   - Alelo: HLA-A*02:01
   - Tamanho: 8-14 resíduos
4. Clique "Processar"
5. Veja resultados nas tabs
6. Gere PDF na tab "Relatório Final"

---

**Pronto para usar! 🚀**
