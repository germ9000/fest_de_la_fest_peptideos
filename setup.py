# -*- coding: utf-8 -*-
"""
Script de setup para instalação e configuração do ambiente.
"""

import subprocess
import sys
import os


def run_command(command, description):
    """Executa comando e trata erros."""
    print(f"\n{'='*60}")
    print(f"📦 {description}")
    print(f"{'='*60}")
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro: {e.stderr}")
        return False


def main():
    """Função principal de setup."""
    print("""
    ╔══════════════════════════════════════════════════════════╗
    ║   Setup - Pipeline de Análise de Peptídeos              ║
    ║   Aplicação Streamlit                                    ║
    ╚══════════════════════════════════════════════════════════╝
    """)
    
    # 1. Instalar dependências Python
    print("\n🔧 Passo 1: Instalando dependências Python...")
    if not run_command(
        f"{sys.executable} -m pip install -r requirements.txt",
        "Instalando pacotes do requirements.txt"
    ):
        print("❌ Falha na instalação de dependências.")
        return False
    
    # 2. Baixar modelos MHCflurry
    print("\n🧬 Passo 2: Baixando modelos MHCflurry...")
    print("⚠️  Isso pode levar alguns minutos na primeira vez...")
    if not run_command(
        "mhcflurry-downloads fetch models_class1_presentation",
        "Baixando modelos MHC-I"
    ):
        print("⚠️  Aviso: Falha ao baixar modelos MHC-I. Você pode fazer isso manualmente depois.")
    
    # Tenta baixar modelos MHC-II também
    run_command(
        "mhcflurry-downloads fetch models_class2_presentation",
        "Baixando modelos MHC-II (opcional)"
    )
    
    # 3. Verificar instalação
    print("\n✅ Passo 3: Verificando instalação...")
    try:
        import streamlit
        import pandas
        import numpy
        import biopython
        print("✅ Streamlit instalado")
        print("✅ Pandas instalado")
        print("✅ NumPy instalado")
        print("✅ Biopython instalado")
        
        try:
            from mhcflurry import Class1PresentationPredictor
            print("✅ MHCflurry instalado")
        except:
            print("⚠️  MHCflurry não encontrado. Execute: pip install mhcflurry")
        
        try:
            import plotly
            print("✅ Plotly instalado")
        except:
            print("⚠️  Plotly não encontrado. Execute: pip install plotly")
        
    except ImportError as e:
        print(f"❌ Erro na verificação: {e}")
        return False
    
    print("""
    ╔══════════════════════════════════════════════════════════╗
    ║   ✅ Setup concluído com sucesso!                       ║
    ╚══════════════════════════════════════════════════════════╝
    
    Para executar a aplicação:
    
        streamlit run app.py
    
    Ou use:
    
        python -m streamlit run app.py
    
    A aplicação abrirá automaticamente em http://localhost:8501
    """)
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
