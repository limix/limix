"""
Quantitative trait locus analysis.
"""

from ._iscan import iscan
from ._scan import scan
from ._st_sscan import st_sscan

__all__ = ["st_sscan", "scan", "iscan"]
