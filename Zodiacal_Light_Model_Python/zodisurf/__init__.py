"""
ZODI Model - A Python package for modeling and analyzing Zodiacal Light data.
"""

from .get_zmod import get_zmod, read_zpars, mk_zdata, set_default_engine, get_default_engine
from .zkernelpy import zkernel
from .solar_sp import solar_sp
from .skysurf_params import get_albedo, get_hong_params, get_mult, get_emiss, put_zpar

__version__ = "0.3.0"  # Bumped for new vectorization/threading!
__author__ = "Rosalia O'Brien, Tejovrash Acharya, Claude 5.5 Opus"
__email__ = "rosalia.d.obrien@gmail.com" 

__all__ = [
    'get_zmod', 'read_zpars', 'mk_zdata', 'set_default_engine', 'get_default_engine',
    'zkernel', 'solar_sp', 
    'get_albedo', 'get_hong_params', 'get_mult', 'get_emiss', 'put_zpar'
]