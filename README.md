# Group Additivity
Python package for group additivity scheme use. 
Currently, Benson's GA scheme is employed.

**Required Packages**
- RDkit (http://www.rdkit.org/)

**Citations**
- Rangarajan et al. "Language-oriented rule-based reaction network generation and analysis: Algorithms of RING", Comput. Chem. Eng. 2014, 64, 124
- Rangarajan et al. "Language-oriented rule-based reaction network generation and analysis: Descrpition of RING", Comput. Chem. Eng. 2012, 45, 114
- Benson et al. "Additivity rules for the estimation of thermochemical properties." Chem. Rev., 1969, 69 (3), 279-324
- Salciccioli et al. "Density Functional Theory-Derived Group Additivity and Linear Scaling Methods for Prediction of Oxygenate Stability on Metal Catalysts: Adsorption of Open-Ring Alcohol and Polyol Dehydrogenation Intermediates on Pt-Based Metals." J. Phys. Chem. C, 2010, 114 (47) 20155-20166
- Salciccioli et al. "Adsorption of Acid, Ester, and Ether Functional Groups on Pt: Fast Prediction of Thermochemical Properties of Adsorbed Oxygenates via DFT-Based Group Additivity Methods." J. Phys. Chem. C, 2012, 116(2), 1873-1886
- Vorotnikov et al. "Group Additivity for Estimating Thermochemical Properties of Furanic Compounds on Pd(111)." Ind. Eng. Chem. Res., 2014, 53 (30), 11929-11938
- Vorotnikov et al. "Group Additivity and Modified Linear Scaling Relations for Estimating Surface Thermochemistry on Transition Metal Surfaces: Application to Furanics." J. Phys. Chem. C, 2015, 119 (19), 10417-10426
- Gu et al. "Group Additivity for Thermochemical Property Estimation of Lignin Monomers on Pt(111)." J. Phys. Chem. C, 2016, 120 (34), 19234-19241

**Example**  
Benson's Gas Group Additivity Example:
```
In:
from VGA.GroupAdd.Library import GroupLibrary
import VGA.ThermoChem
lib = GroupLibrary.Load('benson')
descriptors = lib.GetDescriptors('C1CO1')
print descriptors
thermochem = lib.Estimate(descriptors,'thermochem')
print thermochem.eval_ND_H(298.15)

Out:
defaultdict(<type 'int'>, {'O(C)2': 1, 'C(C)(H)2(O)': 2, 'C1CO1': 1})
-19.9132141547
```
Salciccioli et al. J. Phys. Chem. C, 2012, 116 (2), pp 1873–1886 Example:
```
In:
from VGA.GroupAdd.Library import GroupLibrary
import VGA.ThermoChem
lib = GroupLibrary.Load('Salciccioli2012')
descriptors = lib.GetDescriptors('C([Pt])C[Pt]')
print descriptors
thermochem = lib.Estimate(descriptors,'thermochem')
print thermochem.eval_ND_H(298.15)

Out:
defaultdict(<type 'int'>, {'surface-ring strain': 0.217, 'C(C)(H)2(Pt)': 2})
37.6249461725
```
