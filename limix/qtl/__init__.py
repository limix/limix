""" Quantitative trait locus.

Functions
---------
- scan

Classes
-------
- QTLModel
"""


from ._scan import scan
from ._iscan import iscan
from ._model import QTLModel

__all__ = ["scan", "iscan", "QTLModel"]
