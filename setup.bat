@echo off
REM Script de setup para Windows

echo ==========================================
echo   Setup - Pipeline de Analise de Peptideos
echo ==========================================
echo.

REM Instala dependências Python
echo 📦 Instalando dependencias Python...
python -m pip install -r requirements.txt

REM Baixa modelos MHCflurry
echo.
echo 🧬 Baixando modelos MHCflurry...
echo ⚠️  Isso pode levar alguns minutos...
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads fetch models_class2_presentation

echo.
echo ✅ Setup concluido!
echo.
echo Para executar a aplicacao:
echo   streamlit run app.py
echo.

pause
