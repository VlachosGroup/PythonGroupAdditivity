r"""
===============================================================
RING notation for chemical structures  (:mod:`chemtk.smiles`)
===============================================================
This is a custom written abstract syntax tree parser for RING. Implementation
is based on the '<http://research.cems.umn.edu/bhan/#software>'.

--------
Examples
--------

>>> from chemtk import smiles
>>> chem = smiles.read('O')  # Read SMILES notation for water.
>>> for atom in chem.atoms:
...     print atom
O
H
H
>>> 


---------
Reference
---------

.. autodoc puts stuff extracted from docstrings here.
"""

from Reader import Read

__all__ = ['Read']

del Reader
