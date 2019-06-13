from pgradd.GroupAdd.Library import GroupLibrary
import pgradd.ThermoChem
lib = GroupLibrary.Load('GuSolventGA2017Aq')
descriptors = lib.GetDescriptors('C(=O)([Pt])O')
#descriptors = lib.GetDescriptors('C=C')
print(descriptors)
thermochem = lib.Estimate(descriptors,'thermochem')
print((thermochem.get_HoRT(500)))