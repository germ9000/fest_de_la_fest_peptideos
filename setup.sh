#!/bin/bash
# Script de setup para Linux/Mac

echo "=========================================="
echo "  Setup - Pipeline de Análise de Peptídeos"
echo "=========================================="

# Instala dependências Python
echo ""
echo "📦 Instalando dependências Python..."
pip install -r requirements.txt

# Baixa modelos MHCflurry
echo ""
echo "🧬 Baixando modelos MHCflurry..."
echo "⚠️  Isso pode levar alguns minutos..."
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads fetch models_class2_presentation

echo ""
echo "✅ Setup concluído!"
echo ""
echo "Para executar a aplicação:"
echo "  streamlit run app.py"
echo ""
