""" Quantitative trait locus.

Functions
---------
- st_scan

Classes
-------
- QTLResult
- GWAS_LMM
- GWAS_MTLMM
"""


from ._result import QTLResult
from ._mt_scan import mt_scan
from ._st_iscan import st_iscan
from ._st_scan import st_scan
from ._st_sscan import st_sscan

__all__ = ["st_scan", "st_iscan", "st_sscan", "QTLResult", "mt_scan"]
