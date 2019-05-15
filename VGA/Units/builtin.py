from .db import units_db
from .qty import Quantity, FundamentalUnits, eval_qty


# Define base SI units.
base_SI_units = [
  # (<base name>, <multiple>, <primitive name>)
  ('m', 1.0, 'm'),
  ('g', 0.001, 'kg'),
  ('s', 1.0, 's'),
  ('A', 1.0, 'A'),
  ('K', 1.0, 'K'),
  ('mol', 1.0, 'mol'),
  ('cd', 1.0, 'cd')
]
for base, mult, prim in base_SI_units:
    units_db.add(base, Quantity(mult, FundamentalUnits.new(prim)))

derived_SI_units = [
    ('N', 'kg m/s^2'), ('Pa', 'N/m^2'), ('J', 'N m'), ('W', 'J/s'),
    ('C', 'A s'), ('V', 'J/C'), ('F', 'C/V'), ('Ohm', 'V/A'),
    # etc.  (may define others later)
]
for name, val in derived_SI_units:
    units_db.add(name, eval_qty(val))

other_units = [
    # Includes non-SI metric units and English units.
    # Count:
    ('molecule', '1/(6.02214179*10^23) mol'),
    # Length:
    ('in', '2.54 cm'),
    ('ft', '12 in'),
    # Volume:
    ('L', '100 cm^3'),
    # Time:
    ('min', '60 s'),
    ('h', '60 min'),
    # Mass:
    ('u', '1.660538921*10^-27 kg'),
    ('lb', '0.45359237 kg'),
    ('t', '1000 kg'),
    # Force:
    ('dyn', '10^-5 N'),
    ('lbf', '4.44822162 N'),
    # Pressure:
    ('bar', '100.0 kPa'),
    ('atm', '101.325 kPa'),
    ('torr', '1.0/760 atm'),
    ('psi', 'lbf/in^2'),
    # Energy:
    ('cal', '4.184 J'),
    ('erg', 'dyn cm'),
    ('BTU', '1.05435026444 kJ'),
    ('eV', '1.602176487*10^-19 C V'),
    # Power:
    ('hp', '33000 ft lbf/min'),
    # Viscosity:
    ('P', 'g/(cm s)'),
    ('St', 'cm^2/s'),
]
for name, val in other_units:
    units_db.add(name, eval_qty(val))
