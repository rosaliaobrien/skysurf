"""
ZODI Model - A Python package for modeling and analyzing Zodiacal Light data.
"""

from .get_zmod import get_zmod, read_zpars, mk_zdata
from .zkernelpy import zkernel
from .solar_sp import solar_sp
from .skysurf_params import get_albedo, get_hong_params, get_mult, get_emiss, put_zpar

__version__ = "0.2.1"
__author__ = "Rosalia O'Brien, Tejovrash Acharya"
__email__ = "your.email@domain.com"

__all__ = ['get_zmod', 'read_zpars', 'mk_zdata', 'zkernel', 'solar_sp', 
           'get_albedo', 'get_hong_params', 'get_mult', 'get_emiss', 'put_zpar'] 