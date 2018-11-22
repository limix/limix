""" Quantitative trait locus.

Functions
---------
- scan

Classes
-------
- QTLModel
- GWAS_LMM 
- GWAS_MTLMM 
"""


from ._scan import scan
from ._iscan import iscan
from ._model import QTLModel
from .gwas import GWAS_LMM, GWAS_MTLMM 

__all__ = ["scan", "iscan", "QTLModel", "GWAS_LMM", "GWAS_MTLMM"]
