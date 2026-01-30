# -*- coding: utf-8 -*-
"""
Pipeline Modular de Análise de Peptídeos - Versão Refatorada
=============================================================

Este script demonstra o uso dos módulos modulares:
- data_handler.py: Manipulação de dados e validação
- api_client.py: Chamadas paralelas às APIs (IEDB/UniProt)
- report_gen.py: Geração de relatórios PDF

Autor: Engenheiro de Bioinformática Sênior
Data: 2026-01-30
"""

import sys
import os
from pathlib import Path

# Importa módulos locais
from data_handler import (
    load_peptides, 
    validate_peptides, 
    add_physchem_properties
)
from api_client import APIClientManager
from report_gen import generate_report

# Imports para pipeline MHCflurry (mantido para compatibilidade)
import pandas as pd
import numpy as np
from mhcflurry import Class1PresentationPredictor


def run_mhcflurry_predictions(df: pd.DataFrame, alleles: list) -> pd.DataFrame:
    """
    Executa predições usando MHCflurry (local, não requer API).
    
    Args:
        df: DataFrame com coluna 'peptide'
        alleles: Lista de alelos HLA (ex: ['HLA-A*02:01'])
        
    Returns:
        DataFrame enriquecido com predições MHCflurry
    """
    print(f"Rodando predição MHCflurry para {alleles}...")
    
    predictor = Class1PresentationPredictor.load()
    
    # Predição para todos os alelos
    preds = []
    for allele in alleles:
        results = predictor.predict(
            peptides=df['peptide'].values,
            alleles=[allele],
            include_affinity_percentile=True,
            verbose=0
        )
        results['allele'] = allele
        preds.append(results)
    
    # Combina resultados e pega o melhor alelo para cada peptídeo
    full_preds = pd.concat(preds)
    best_preds = full_preds.loc[full_preds.groupby('peptide')['presentation_score'].idxmax()]
    
    # Merge com DataFrame original
    df_result = pd.merge(
        df, 
        best_preds[['peptide', 'affinity', 'presentation_score', 'allele']], 
        on='peptide', 
        how='left'
    )
    
    # Renomeia colunas para consistência
    if 'affinity' in df_result.columns:
        df_result['affinity_nm'] = df_result['affinity']
    
    if 'presentation_score' in df_result.columns:
        df_result['mhc_score'] = df_result['presentation_score']
    
    # Calcula percentile rank se disponível
    if 'presentation_percentile' in best_preds.columns:
        df_result['percentile_rank'] = df_result['peptide'].map(
            best_preds.set_index('peptide')['presentation_percentile']
        )
    else:
        # Calcula rank relativo
        df_result['percentile_rank'] = df_result['mhc_score'].rank(pct=True, ascending=False)
    
    return df_result


def calculate_final_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calcula scores finais normalizados para ranking.
    
    Args:
        df: DataFrame com predições
        
    Returns:
        DataFrame com scores finais calculados
    """
    df = df.copy()
    
    # Normaliza scores
    from sklearn.preprocessing import MinMaxScaler
    scaler = MinMaxScaler()
    
    # Normaliza score de apresentação MHC
    if 'mhc_score' in df.columns:
        try:
            df['norm_mhc'] = scaler.fit_transform(df[['mhc_score']].fillna(0))
        except:
            df['norm_mhc'] = 0.5
    
    # Score final (pode ser ajustado conforme necessidade)
    if 'norm_mhc' in df.columns:
        df['final_rank_score'] = df['norm_mhc']
    else:
        df['final_rank_score'] = 0.0
    
    # Ordena por score final
    df = df.sort_values(by='final_rank_score', ascending=False).reset_index(drop=True)
    
    return df


def main():
    """Função principal do pipeline."""
    
    print("=" * 70)
    print("PIPELINE MODULAR DE ANÁLISE DE PEPTÍDEOS")
    print("=" * 70)
    print()
    
    # ========================================================================
    # CONFIGURAÇÃO
    # ========================================================================
    
    # Caminho do arquivo de entrada (ajustar conforme necessário)
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        # Solicita arquivo interativamente
        input_file = input("Digite o caminho do arquivo de peptídeos: ").strip()
    
    if not os.path.exists(input_file):
        print(f"❌ Erro: Arquivo não encontrado: {input_file}")
        return
    
    # Alelos HLA para análise
    ALLELES = ["HLA-A*02:01"]  # Pode ser expandido: ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02"]
    
    # Configurações de API
    USE_API_ENRICHMENT = True  # Se False, usa apenas MHCflurry local
    MAX_WORKERS = 5  # Número de threads para processamento paralelo
    
    # ========================================================================
    # ETAPA 1: CARREGAMENTO E VALIDAÇÃO DE DADOS
    # ========================================================================
    
    print("\n[ETAPA 1] Carregando e validando dados...")
    print("-" * 70)
    
    try:
        # Carrega peptídeos usando data_handler
        df_input = load_peptides(input_file)
        total_peptides = len(df_input)
        print(f"✓ {total_peptides} peptídeos carregados do arquivo")
        
        # Valida peptídeos
        df_valid = validate_peptides(df_input, min_length=8, max_length=14)
        valid_count = len(df_valid)
        print(f"✓ {valid_count} peptídeos válidos após validação")
        print(f"  Taxa de validação: {(valid_count/total_peptides*100):.1f}%")
        
        if df_valid.empty:
            print("❌ Nenhum peptídeo válido encontrado após validação.")
            return
        
    except Exception as e:
        print(f"❌ Erro ao carregar/validar dados: {e}")
        return
    
    # ========================================================================
    # ETAPA 2: PROPRIEDADES FÍSICO-QUÍMICAS
    # ========================================================================
    
    print("\n[ETAPA 2] Calculando propriedades físico-químicas...")
    print("-" * 70)
    
    try:
        df_result = add_physchem_properties(df_valid)
        print(f"✓ Propriedades calculadas para {len(df_result)} peptídeos")
    except Exception as e:
        print(f"❌ Erro ao calcular propriedades: {e}")
        return
    
    # ========================================================================
    # ETAPA 3: PREDIÇÕES MHC (MHCflurry Local)
    # ========================================================================
    
    print("\n[ETAPA 3] Executando predições MHCflurry...")
    print("-" * 70)
    
    try:
        df_result = run_mhcflurry_predictions(df_result, ALLELES)
        print(f"✓ Predições MHCflurry concluídas")
        
        if 'affinity_nm' in df_result.columns:
            best_affinity = df_result['affinity_nm'].min()
            avg_affinity = df_result['affinity_nm'].mean()
            print(f"  Melhor afinidade: {best_affinity:.2f} nM")
            print(f"  Afinidade média: {avg_affinity:.2f} nM")
    except Exception as e:
        print(f"❌ Erro nas predições MHCflurry: {e}")
        return
    
    # ========================================================================
    # ETAPA 4: ENRIQUECIMENTO COM APIs EXTERNAS (Opcional)
    # ========================================================================
    
    if USE_API_ENRICHMENT:
        print("\n[ETAPA 4] Enriquecendo dados com APIs externas (IEDB/UniProt)...")
        print("-" * 70)
        print("⚠️  Nota: Esta etapa pode demorar devido a rate limiting das APIs")
        
        try:
            api_manager = APIClientManager(max_workers=MAX_WORKERS, request_delay=0.2)
            df_result = api_manager.enrich_dataframe(
                df_result, 
                allele=ALLELES[0],
                include_uniprot=False  # Pode ser True para incluir busca UniProt
            )
            print("✓ Enriquecimento concluído")
        except Exception as e:
            print(f"⚠️  Aviso: Erro no enriquecimento de APIs (continuando sem dados de API): {e}")
            # Continua sem dados de API
    
    # ========================================================================
    # ETAPA 5: CÁLCULO DE SCORES FINAIS
    # ========================================================================
    
    print("\n[ETAPA 5] Calculando scores finais...")
    print("-" * 70)
    
    try:
        df_result = calculate_final_scores(df_result)
        print(f"✓ Scores finais calculados")
    except Exception as e:
        print(f"❌ Erro ao calcular scores finais: {e}")
        return
    
    # ========================================================================
    # ETAPA 6: EXPORTAÇÃO DE RESULTADOS
    # ========================================================================
    
    print("\n[ETAPA 6] Exportando resultados...")
    print("-" * 70)
    
    # Gera nome do arquivo de saída
    input_path = Path(input_file)
    output_base = input_path.stem
    
    # Exporta para Excel
    excel_output = f"{output_base}_RESULTADO.xlsx"
    try:
        df_result.to_excel(excel_output, index=False)
        print(f"✓ Dados exportados para Excel: {excel_output}")
    except Exception as e:
        print(f"⚠️  Erro ao exportar Excel: {e}")
    
    # Gera relatório PDF
    pdf_output = f"{output_base}_RELATORIO.pdf"
    try:
        generate_report(
            df_result,
            pdf_output,
            total_peptides=total_peptides,
            include_statistics=True,
            max_table_rows=50
        )
        print(f"✓ Relatório PDF gerado: {pdf_output}")
    except Exception as e:
        print(f"⚠️  Erro ao gerar PDF: {e}")
    
    # ========================================================================
    # RESUMO FINAL
    # ========================================================================
    
    print("\n" + "=" * 70)
    print("RESUMO FINAL")
    print("=" * 70)
    print(f"Total de peptídeos processados: {len(df_result)}")
    
    if 'affinity_nm' in df_result.columns:
        print(f"Melhor afinidade: {df_result['affinity_nm'].min():.2f} nM")
        print(f"Afinidade média: {df_result['affinity_nm'].mean():.2f} nM")
    
    print("\n🏆 TOP 5 CANDIDATOS:")
    top5 = df_result.head(5)
    display_cols = ['peptide']
    if 'affinity_nm' in df_result.columns:
        display_cols.append('affinity_nm')
    if 'mhc_score' in df_result.columns:
        display_cols.append('mhc_score')
    if 'final_rank_score' in df_result.columns:
        display_cols.append('final_rank_score')
    
    print(top5[display_cols].to_string(index=False))
    print("\n✅ Pipeline concluído com sucesso!")
    print("=" * 70)


if __name__ == "__main__":
    main()
