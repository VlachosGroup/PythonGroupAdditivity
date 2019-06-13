r"""
=====================================
RING notation for chemical structures
=====================================
This is a custom written abstract syntax tree parser for RING. Implementation
is based on the '<http://research.cems.umn.edu/bhan/#software>'.
Rangarajan, S., Kaminski, T., Van Wyk, E., Bhan A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Algorithms of RING", Computers and Chemical Engineering 64 (2014) 124,
    10.1016/j.compchemeng.2014.02.007
Rangarajan, S., Bhan, A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Descrpition of RING", Computers and Chemical Engineering 45 (2012) 114,
    10.1016/j.compchemeng.2012.06.008
"""

from .Reader import Read

__all__ = ['Read']
