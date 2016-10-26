from VGA.GroupAdd.Library import GroupLibrary
import VGA.ThermoChem
lib = GroupLibrary.Load('benson')
descriptors = lib.GetDescriptors('C1CO1')
print descriptors
thermochem = lib.Estimate(descriptors,'thermochem')
print thermochem.eval_ND_H(298.15)

from VGA.GroupAdd.Library import GroupLibrary
import VGA.ThermoChem
lib = GroupLibrary.Load('Salciccioli2012')
descriptors = lib.GetDescriptors('C([Pt])C[Pt]')
print descriptors
thermochem = lib.Estimate(descriptors,'thermochem')
print thermochem.eval_ND_H(298.15)