# -*- coding: utf-8 -*-
"""
Script de teste para verificar instalação.
Execute: python test_installation.py
"""

import sys

def test_imports():
    """Testa imports de todos os módulos necessários."""
    modules = {
        'streamlit': 'Streamlit',
        'pandas': 'Pandas',
        'numpy': 'NumPy',
        'Bio': 'Biopython',
        'mhcflurry': 'MHCflurry',
        'plotly': 'Plotly',
        'fpdf': 'fpdf2',
        'PIL': 'Pillow',
        'sklearn': 'scikit-learn',
        'requests': 'Requests'
    }
    
    print("=" * 60)
    print("TESTE DE INSTALAÇÃO - Pipeline de Análise de Peptídeos")
    print("=" * 60)
    print()
    
    failed = []
    for module, name in modules.items():
        try:
            __import__(module)
            print(f"✅ {name:20} - OK")
        except ImportError as e:
            print(f"❌ {name:20} - FALTANDO")
            failed.append(name)
    
    print()
    print("=" * 60)
    
    if failed:
        print(f"❌ Módulos faltando: {', '.join(failed)}")
        print()
        print("Execute para instalar:")
        print("  pip install -r requirements.txt")
        print()
        return False
    else:
        print("✅ Todos os módulos instalados corretamente!")
        print()
        
        # Testa módulos locais
        print("Testando módulos locais...")
        local_modules = [
            'data_handler',
            'api_client',
            'report_gen',
            'mhc_predictions',
            'conservation_analysis',
            'netmhcpan_client'
        ]
        
        local_failed = []
        for module in local_modules:
            try:
                __import__(module)
                print(f"✅ {module:25} - OK")
            except ImportError as e:
                print(f"❌ {module:25} - ERRO: {e}")
                local_failed.append(module)
        
        if local_failed:
            print()
            print(f"❌ Módulos locais com problemas: {', '.join(local_failed)}")
            print("Verifique se os arquivos estão no diretório correto.")
            return False
        
        print()
        print("=" * 60)
        print("✅ INSTALAÇÃO COMPLETA E FUNCIONAL!")
        print("=" * 60)
        print()
        print("Para executar a aplicação:")
        print("  streamlit run app.py")
        print()
        return True

if __name__ == "__main__":
    success = test_imports()
    sys.exit(0 if success else 1)
