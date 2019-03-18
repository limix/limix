"""
Quantitative trait locus analysis.
"""


from ._mt_scan import mt_scan
from ._st_iscan import st_iscan
from ._st_scan import st_scan
from ._st_sscan import st_sscan

__all__ = ["st_scan", "st_iscan", "st_sscan", "mt_scan"]
