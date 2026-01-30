# -*- coding: utf-8 -*-
"""
Pipeline de Análise de Peptídeos - Módulos Principais
"""

__version__ = "2.0.0"
__author__ = "Engenheiro de Bioinformática Sênior"

from .data_processor import FastaProcessor
from .mhc_analyzer import MHCAnalyzer
from .api_manager import APIManager
from .pdf_generator import PDFGenerator
from .conservation import ConservationAnalyzer

__all__ = [
    'FastaProcessor',
    'MHCAnalyzer',
    'APIManager',
    'PDFGenerator',
    'ConservationAnalyzer',
]
