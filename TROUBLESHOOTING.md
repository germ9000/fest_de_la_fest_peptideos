# 🔧 Guia de Troubleshooting

## Problemas Comuns e Soluções

### ❌ Erro: "ModuleNotFoundError: No module named 'X'"

**Causa**: Módulo não instalado

**Solução**:
```bash
pip install -r requirements.txt
```

Se persistir, instale manualmente:
```bash
pip install streamlit pandas numpy biopython mhcflurry plotly fpdf2 Pillow scikit-learn requests
```

---

### ❌ Erro: "MHCflurry models not found"

**Causa**: Modelos não foram baixados

**Solução**:
```bash
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads fetch models_class2_presentation
```

**Verificação**:
```python
from mhcflurry import Class1PresentationPredictor
predictor = Class1PresentationPredictor.load()  # Deve funcionar sem erro
```

---

### ❌ Erro: "Streamlit não encontrado"

**Causa**: Streamlit não instalado ou ambiente virtual não ativado

**Solução**:
```bash
# Instale Streamlit
pip install streamlit

# Verifique instalação
streamlit --version

# Execute
streamlit run app.py
```

---

### ❌ Erro ao importar módulos locais

**Causa**: Caminho incorreto ou módulos não encontrados

**Solução**:
1. Certifique-se de estar no diretório correto:
```bash
cd "C:\Users\Gabriel Eduardo\Desktop\analises finais peptideos"
```

2. Verifique se todos os arquivos estão presentes:
   - `app.py`
   - `data_handler.py`
   - `api_client.py`
   - `report_gen.py`
   - `mhc_predictions.py`
   - `conservation_analysis.py`
   - `netmhcpan_client.py`

3. Execute Python para testar imports:
```python
python -c "from data_handler import load_peptides_from_bytes; print('OK')"
```

---

### ❌ Erro: "Port 8501 already in use"

**Causa**: Outra instância do Streamlit está rodando

**Solução**:
```bash
# Windows
netstat -ano | findstr :8501
taskkill /PID <PID> /F

# Linux/Mac
lsof -ti:8501 | xargs kill -9

# Ou use outra porta
streamlit run app.py --server.port 8502
```

---

### ❌ Erro: "Permission denied" ao executar scripts

**Windows**:
```powershell
# Execute como Administrador ou
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

**Linux/Mac**:
```bash
chmod +x setup.sh
./setup.sh
```

---

### ❌ Erro: "FileNotFoundError" ao ler arquivo

**Causa**: Caminho do arquivo incorreto ou arquivo não existe

**Solução**:
1. Use caminhos absolutos ou relativos corretos
2. Verifique se o arquivo existe:
```python
import os
print(os.path.exists("seu_arquivo.fasta"))
```

---

### ❌ Aplicação muito lenta

**Causas e Soluções**:

1. **Muitos peptídeos**: Processe em lotes menores
2. **APIs externas**: Desative se não necessário
3. **Threads**: Reduza `max_workers` na sidebar
4. **Cache**: Use `@st.cache_data` (já implementado)

**Otimizações**:
- Reduza número de peptídeos por análise
- Desative NetMHCpan se não necessário
- Desative análise de conservação se não necessário
- Use apenas MHC-I ou MHC-II, não ambos

---

### ❌ Erro: "API timeout" ou "Connection error"

**Causa**: Problemas de rede ou APIs indisponíveis

**Solução**:
1. Verifique conexão com internet
2. Aumente `request_delay` na sidebar
3. Desative APIs externas temporariamente
4. APIs podem estar temporariamente indisponíveis

---

### ❌ Erro ao gerar PDF

**Causa**: Problemas com fpdf2 ou Pillow

**Solução**:
```bash
pip install --upgrade fpdf2 Pillow
```

Se persistir, verifique permissões de escrita no diretório.

---

### ❌ Erro: "Bio.Align not found"

**Causa**: Versão antiga do Biopython

**Solução**:
```bash
pip install --upgrade biopython
```

---

### ❌ Erro no Windows: "mhcflurry-downloads não é reconhecido"

**Causa**: Script não está no PATH

**Solução**:
```bash
# Use Python diretamente
python -m mhcflurry.downloads fetch models_class1_presentation
```

Ou adicione ao PATH:
```powershell
$env:Path += ";C:\Users\<SEU_USUARIO>\AppData\Local\Programs\Python\Python<versao>\Scripts"
```

---

## 🧪 Teste Rápido

Execute este script para verificar instalação:

```python
# test_installation.py
import sys

def test_imports():
    modules = [
        'streamlit',
        'pandas',
        'numpy',
        'Bio',
        'mhcflurry',
        'plotly',
        'fpdf',
        'PIL',
        'sklearn',
        'requests'
    ]
    
    failed = []
    for module in modules:
        try:
            __import__(module)
            print(f"✅ {module}")
        except ImportError:
            print(f"❌ {module}")
            failed.append(module)
    
    if failed:
        print(f"\n❌ Módulos faltando: {', '.join(failed)}")
        print("Execute: pip install -r requirements.txt")
        return False
    else:
        print("\n✅ Todos os módulos instalados corretamente!")
        return True

if __name__ == "__main__":
    test_imports()
```

Execute:
```bash
python test_installation.py
```

---

## 📞 Ainda com Problemas?

1. **Verifique logs**: Streamlit mostra erros no terminal
2. **Modo debug**: Execute com `streamlit run app.py --logger.level=debug`
3. **Ambiente virtual**: Use venv para isolar dependências:
```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
venv\Scripts\activate     # Windows
pip install -r requirements.txt
```

4. **Versão Python**: Requer Python 3.8+
```bash
python --version
```

---

## ✅ Checklist de Instalação

- [ ] Python 3.8+ instalado
- [ ] Dependências instaladas (`pip install -r requirements.txt`)
- [ ] Modelos MHCflurry baixados
- [ ] Todos os arquivos .py presentes
- [ ] Porta 8501 disponível
- [ ] Conexão com internet (para APIs)

---

**Boa sorte! 🚀**
