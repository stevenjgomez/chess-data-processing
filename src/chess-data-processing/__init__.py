'''
Stack, orient, and index data from CHESS ID4B.
'''

from _meta import __author__, __copyright__, __license__, __version__
from ._indexing import *
from ._stacking import *

# What to import when running "from package import *"
__all__ = ['stack_em_all', 'stack_tempdep', 'orient_index_tempdep']
