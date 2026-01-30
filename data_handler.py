# -*- coding: utf-8 -*-
"""
Módulo de manipulação de dados para análise de peptídeos.
Responsável por leitura de arquivos FASTA, validação e cálculo de propriedades físico-químicas.
"""

import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from typing import List, Dict, Optional, Union
import io


def read_fasta(filepath: str) -> pd.DataFrame:
    """
    Lê arquivo FASTA e extrai sequências de peptídeos.
    
    Args:
        filepath: Caminho para o arquivo FASTA
        
    Returns:
        DataFrame com coluna 'peptide' contendo as sequências
    """
    peptides = []
    
    try:
        with open(filepath, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()
                # Remove caracteres não-aminoácidos
                seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq)
                if seq:
                    peptides.append(seq)
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo FASTA: {e}")
    
    if not peptides:
        raise ValueError("Nenhuma sequência válida encontrada no arquivo FASTA")
    
    df = pd.DataFrame({'peptide': peptides})
    return df


def read_tsv_csv(filepath: str, sep: Optional[str] = None) -> pd.DataFrame:
    """
    Lê arquivo TSV/CSV e extrai peptídeos.
    
    Args:
        filepath: Caminho para o arquivo
        sep: Separador (None para detecção automática)
        
    Returns:
        DataFrame com coluna 'peptide'
    """
    try:
        if filepath.endswith('.tsv'):
            df = pd.read_csv(filepath, sep='\t', header=None)
        elif filepath.endswith('.csv'):
            df = pd.read_csv(filepath, sep=sep if sep else ',', header=None)
        else:
            # Tenta detectar automaticamente
            df = pd.read_csv(filepath, sep=sep, header=None)
        
        # Lineariza o DataFrame
        peptides_flat = df.values.flatten()
        df_clean = pd.DataFrame(peptides_flat, columns=['peptide'])
        df_clean = df_clean.dropna()
        df_clean['peptide'] = df_clean['peptide'].astype(str).str.strip()
        df_clean = df_clean[df_clean['peptide'] != '']
        
        return df_clean
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo TSV/CSV: {e}")


def read_excel(filepath: str, sheet_name: Optional[Union[str, int]] = 0) -> pd.DataFrame:
    """
    Lê arquivo Excel e extrai peptídeos.
    
    Args:
        filepath: Caminho para o arquivo Excel
        sheet_name: Nome ou índice da planilha
        
    Returns:
        DataFrame com coluna 'peptide'
    """
    try:
        df = pd.read_excel(filepath, sheet_name=sheet_name)
        
        # Tenta encontrar coluna de peptídeos
        if 'peptide' in df.columns:
            df_clean = df[['peptide']].copy()
        else:
            # Procura coluna com sequências válidas
            for col in df.columns:
                if df[col].astype(str).str.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$').any():
                    df_clean = pd.DataFrame({'peptide': df[col]})
                    break
            else:
                # Se não encontrar, usa primeira coluna
                df_clean = pd.DataFrame({'peptide': df.iloc[:, 0]})
        
        df_clean = df_clean.dropna()
        df_clean['peptide'] = df_clean['peptide'].astype(str).str.strip()
        df_clean = df_clean[df_clean['peptide'] != '']
        
        return df_clean
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo Excel: {e}")


def load_peptides(filepath: str) -> pd.DataFrame:
    """
    Carrega peptídeos de diversos formatos (FASTA, TSV, CSV, Excel).
    
    Args:
        filepath: Caminho para o arquivo
        
    Returns:
        DataFrame com coluna 'peptide'
    """
    filepath_lower = filepath.lower()
    
    if filepath_lower.endswith(('.fasta', '.fa', '.faa')):
        return read_fasta(filepath)
    elif filepath_lower.endswith(('.tsv', '.txt')):
        return read_tsv_csv(filepath, sep='\t')
    elif filepath_lower.endswith('.csv'):
        return read_tsv_csv(filepath, sep=',')
    elif filepath_lower.endswith(('.xlsx', '.xls')):
        return read_excel(filepath)
    else:
        # Tenta como TSV por padrão
        return read_tsv_csv(filepath)


def read_fasta_from_bytes(file_bytes: bytes) -> pd.DataFrame:
    """
    Lê arquivo FASTA a partir de bytes (útil para Streamlit file_uploader).
    
    Args:
        file_bytes: Conteúdo do arquivo em bytes
        
    Returns:
        DataFrame com coluna 'peptide'
    """
    peptides = []
    
    try:
        file_string = file_bytes.decode('utf-8')
        file_io = io.StringIO(file_string)
        
        for record in SeqIO.parse(file_io, "fasta"):
            seq = str(record.seq).upper()
            seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq)
            if seq:
                peptides.append(seq)
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo FASTA: {e}")
    
    if not peptides:
        raise ValueError("Nenhuma sequência válida encontrada no arquivo FASTA")
    
    return pd.DataFrame({'peptide': peptides})


def read_tsv_csv_from_bytes(file_bytes: bytes, filename: str, sep: Optional[str] = None) -> pd.DataFrame:
    """
    Lê arquivo TSV/CSV a partir de bytes.
    
    Args:
        file_bytes: Conteúdo do arquivo em bytes
        filename: Nome do arquivo (para detectar formato)
        sep: Separador (None para detecção automática)
        
    Returns:
        DataFrame com coluna 'peptide'
    """
    try:
        if filename.lower().endswith('.tsv') or (sep is None and '\t' in file_bytes.decode('utf-8', errors='ignore')[:1000]):
            df = pd.read_csv(io.BytesIO(file_bytes), sep='\t', header=None)
        else:
            df = pd.read_csv(io.BytesIO(file_bytes), sep=sep if sep else ',', header=None)
        
        peptides_flat = df.values.flatten()
        df_clean = pd.DataFrame(peptides_flat, columns=['peptide'])
        df_clean = df_clean.dropna()
        df_clean['peptide'] = df_clean['peptide'].astype(str).str.strip()
        df_clean = df_clean[df_clean['peptide'] != '']
        
        return df_clean
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo TSV/CSV: {e}")


def read_excel_from_bytes(file_bytes: bytes, sheet_name: Optional[Union[str, int]] = 0) -> pd.DataFrame:
    """
    Lê arquivo Excel a partir de bytes.
    
    Args:
        file_bytes: Conteúdo do arquivo em bytes
        sheet_name: Nome ou índice da planilha
        
    Returns:
        DataFrame com coluna 'peptide'
    """
    try:
        df = pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet_name)
        
        if 'peptide' in df.columns:
            df_clean = df[['peptide']].copy()
        else:
            for col in df.columns:
                if df[col].astype(str).str.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$').any():
                    df_clean = pd.DataFrame({'peptide': df[col]})
                    break
            else:
                df_clean = pd.DataFrame({'peptide': df.iloc[:, 0]})
        
        df_clean = df_clean.dropna()
        df_clean['peptide'] = df_clean['peptide'].astype(str).str.strip()
        df_clean = df_clean[df_clean['peptide'] != '']
        
        return df_clean
    except Exception as e:
        raise ValueError(f"Erro ao ler arquivo Excel: {e}")


def load_peptides_from_bytes(file_bytes: bytes, filename: str) -> pd.DataFrame:
    """
    Carrega peptídeos a partir de bytes (para uso com Streamlit file_uploader).
    
    Args:
        file_bytes: Conteúdo do arquivo em bytes
        filename: Nome do arquivo (para detectar formato)
        
    Returns:
        DataFrame com coluna 'peptide'
    """
    filename_lower = filename.lower()
    
    if filename_lower.endswith(('.fasta', '.fa', '.faa')):
        return read_fasta_from_bytes(file_bytes)
    elif filename_lower.endswith(('.tsv', '.txt')):
        return read_tsv_csv_from_bytes(file_bytes, filename, sep='\t')
    elif filename_lower.endswith('.csv'):
        return read_tsv_csv_from_bytes(file_bytes, filename, sep=',')
    elif filename_lower.endswith(('.xlsx', '.xls')):
        return read_excel_from_bytes(file_bytes)
    else:
        # Tenta como TSV por padrão
        return read_tsv_csv_from_bytes(file_bytes, filename, sep='\t')


def validate_peptides(df: pd.DataFrame, min_length: int = 8, max_length: int = 14) -> pd.DataFrame:
    """
    Valida e filtra peptídeos baseado em critérios biológicos.
    
    Args:
        df: DataFrame com coluna 'peptide'
        min_length: Tamanho mínimo do peptídeo (padrão: 8 para MHC-I)
        max_length: Tamanho máximo do peptídeo (padrão: 14 para MHC-I)
        
    Returns:
        DataFrame filtrado com apenas peptídeos válidos
    """
    if 'peptide' not in df.columns:
        raise ValueError("DataFrame deve conter coluna 'peptide'")
    
    df = df.copy()
    
    # Padroniza para maiúsculas
    df['peptide'] = df['peptide'].astype(str).str.upper().str.strip()
    
    # Regex para aminoácidos canônicos
    aa_regex = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")
    
    # Valida sequências
    df["valid"] = df['peptide'].apply(lambda x: bool(aa_regex.match(x)))
    
    # Remove inválidos e duplicatas
    df_valid = df[df["valid"]].drop(columns=["valid"]).drop_duplicates(subset=['peptide']).copy()
    
    # Filtro de tamanho (MHC-I liga em 8-14 meros)
    df_valid = df_valid[df_valid['peptide'].str.len().between(min_length, max_length)]
    
    return df_valid.reset_index(drop=True)


def calculate_kyte_doolittle_hydrophobicity(peptide: str) -> float:
    """
    Calcula hidrofobicidade média usando a escala de Kyte-Doolittle.
    
    Args:
        peptide: Sequência do peptídeo
        
    Returns:
        Hidrofobicidade média (valores positivos = hidrofóbico)
    """
    # Escala de hidrofobicidade de Kyte-Doolittle (1982)
    kd_scale = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5,
        'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
        'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
    
    if not peptide:
        return 0.0
    
    try:
        hydrophobicity_values = [kd_scale.get(aa, 0.0) for aa in peptide.upper()]
        return round(sum(hydrophobicity_values) / len(hydrophobicity_values), 3)
    except:
        return 0.0


def calculate_physchem_properties(peptide: str) -> Dict[str, float]:
    """
    Calcula propriedades físico-químicas de um peptídeo usando Biopython.
    Inclui hidrofobicidade de Kyte-Doolittle.
    
    Args:
        peptide: Sequência do peptídeo
        
    Returns:
        Dicionário com propriedades calculadas
    """
    try:
        pa = ProteinAnalysis(peptide)
        
        # Calcula propriedades básicas
        properties = {
            "length": len(peptide),
            "mw": round(pa.molecular_weight(), 2),
            "pI": round(pa.isoelectric_point(), 2),
            "gravy": round(pa.gravy(), 3),
        }
        
        # Hidrofobicidade Kyte-Doolittle
        properties["kd_hydrophobicity"] = calculate_kyte_doolittle_hydrophobicity(peptide)
        
        # Calcula índice de instabilidade (Guruprasad)
        try:
            properties["instability_index"] = round(pa.instability_index(), 2)
        except:
            properties["instability_index"] = 0.0
        
        # Calcula índice alifático
        try:
            aa_percent = pa.get_amino_acids_percent()
            properties["aliphatic_index"] = round(
                (aa_percent.get('A', 0) + 2.9 * aa_percent.get('V', 0) + 
                 3.9 * (aa_percent.get('I', 0) + aa_percent.get('L', 0))) * 100, 2
            )
        except:
            properties["aliphatic_index"] = 0.0
        
        return properties
    except Exception as e:
        # Retorna valores padrão em caso de erro
        return {
            "length": len(peptide),
            "mw": 0.0,
            "pI": 0.0,
            "gravy": 0.0,
            "kd_hydrophobicity": 0.0,
            "instability_index": 0.0,
            "aliphatic_index": 0.0
        }


def add_physchem_properties(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adiciona propriedades físico-químicas ao DataFrame de peptídeos.
    
    Args:
        df: DataFrame com coluna 'peptide'
        
    Returns:
        DataFrame com colunas adicionais de propriedades
    """
    if 'peptide' not in df.columns:
        raise ValueError("DataFrame deve conter coluna 'peptide'")
    
    # Calcula propriedades para cada peptídeo
    props = df['peptide'].apply(calculate_physchem_properties).apply(pd.Series)
    
    # Remove colunas duplicadas antes de concatenar
    cols_to_remove = [c for c in props.columns if c in df.columns]
    if cols_to_remove:
        df = df.drop(columns=cols_to_remove)
    
    # Concatena propriedades
    df_result = pd.concat([df.reset_index(drop=True), props.reset_index(drop=True)], axis=1)
    
    return df_result
