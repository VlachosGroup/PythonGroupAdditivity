"""
==================
Physical constants
==================

Values of common physical constants are defined here as instances of
:class:`pgradd.Units.Quantity`.

"""

from . Units import eval_qty


#: Planck's constant
PLANCK_CONSTANT = eval_qty('6.626068*10^-34 J s')

#: Boltzman's constant
BOLTZMANN_CONSTANT = eval_qty('1.3806503*10^-23 J/K')

#: Boltzman's constant
GAS_CONSTANT = eval_qty('8.314472 J/(mol K)')

#: Avogadro's number
AVOGADRO_NUMBER = eval_qty('6.02214179*10^23 mol^-1')
