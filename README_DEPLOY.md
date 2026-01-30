# 🚀 Guia de Deploy - Streamlit Cloud

## Deploy no Streamlit Cloud (Gratuito)

### Pré-requisitos
1. Conta no GitHub
2. Repositório com o código
3. Conta no [Streamlit Cloud](https://streamlit.io/cloud)

### Passo a Passo

#### 1. Preparar Repositório GitHub

```bash
# Inicialize git (se ainda não fez)
git init
git add .
git commit -m "Initial commit - Streamlit app"

# Crie repositório no GitHub e conecte
git remote add origin https://github.com/SEU_USUARIO/SEU_REPO.git
git branch -M main
git push -u origin main
```

#### 2. Deploy no Streamlit Cloud

1. Acesse [share.streamlit.io](https://share.streamlit.io)
2. Faça login com GitHub
3. Clique em "New app"
4. Configure:
   - **Repository**: Seu repositório GitHub
   - **Branch**: `main` (ou sua branch principal)
   - **Main file path**: `app.py`
   - **App URL**: Escolha um nome único

5. Clique em "Deploy"

#### 3. Configurações Importantes

O Streamlit Cloud detecta automaticamente:
- `requirements.txt` - Instala dependências Python
- `packages.txt` - Instala pacotes do sistema (se necessário)
- `.streamlit/config.toml` - Configurações da aplicação

### ⚠️ Problemas Comuns e Soluções

#### Problema: "Module not found"
**Solução**: Verifique se todas as dependências estão em `requirements.txt`

#### Problema: "MHCflurry models not found"
**Solução**: Adicione ao `setup.sh`:
```bash
#!/bin/bash
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads fetch models_class2_presentation
```

#### Problema: "Timeout during installation"
**Solução**: 
- Reduza dependências pesadas
- Use `packages.txt` para dependências do sistema
- Considere usar cache do Streamlit

### 📝 Arquivos Necessários para Deploy

```
seu-repo/
├── app.py                    # ✅ Aplicação principal
├── requirements.txt          # ✅ Dependências Python
├── packages.txt              # ✅ Pacotes do sistema (opcional)
├── .streamlit/
│   └── config.toml          # ✅ Configurações
├── data_handler.py          # ✅ Módulos
├── api_client.py            # ✅ Módulos
├── report_gen.py            # ✅ Módulos
├── mhc_predictions.py       # ✅ Módulos
├── conservation_analysis.py # ✅ Módulos
└── netmhcpan_client.py      # ✅ Módulos
```

### 🔧 Configuração Avançada

#### Variáveis de Ambiente (se necessário)

No Streamlit Cloud, vá em "Settings" > "Secrets" e adicione:

```toml
[api_keys]
iedb_key = "sua_chave_aqui"
uniprot_key = "sua_chave_aqui"
```

Acesse no código com:
```python
import streamlit as st
api_key = st.secrets["api_keys"]["iedb_key"]
```

### 📊 Limites do Streamlit Cloud

- **RAM**: 1GB (gratuito)
- **CPU**: Compartilhado
- **Storage**: 1GB
- **Timeout**: 30 segundos por requisição

### 💡 Dicas de Performance

1. **Use cache**: `@st.cache_data` para dados pesados
2. **Limite tamanho**: Processe arquivos menores (<10MB)
3. **Otimize imports**: Importe apenas o necessário
4. **Use progress bars**: Mostre progresso para usuário

### 🐛 Troubleshooting

#### App não inicia
- Verifique logs em "Manage app" > "Logs"
- Confirme que `app.py` está na raiz
- Verifique sintaxe Python

#### Erro de memória
- Reduza `max_workers` no código
- Processe menos peptídeos por vez
- Desative APIs externas se não necessário

#### Erro de timeout
- Reduza timeout de APIs
- Processe em lotes menores
- Use cache agressivamente

### 📞 Suporte

- [Documentação Streamlit Cloud](https://docs.streamlit.io/streamlit-cloud)
- [Fórum Streamlit](https://discuss.streamlit.io/)
- [GitHub Issues](https://github.com/streamlit/streamlit/issues)

---

## 🖥️ Execução Local

### Windows

```powershell
# Instale dependências
pip install -r requirements.txt

# Baixe modelos MHCflurry
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads fetch models_class2_presentation

# Execute
streamlit run app.py
```

### Linux/Mac

```bash
# Instale dependências
pip install -r requirements.txt

# Baixe modelos MHCflurry
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads fetch models_class2_presentation

# Execute
streamlit run app.py
```

### Usando setup.py

```bash
python setup.py
streamlit run app.py
```

---

**Boa sorte com seu deploy! 🚀**
