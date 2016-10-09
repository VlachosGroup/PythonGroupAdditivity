# Group Additivity
Python package for group additivity scheme use. 
Currently, Benson's GA scheme is employed.

**Requires**
RDkit (http://www.rdkit.org/)

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
```
from VGA.GroupAdd.Library import GroupLibrary
import VGA.ThermoChem
lib = GroupLibrary.Load('benson')
groups = lib.GetGroups('CC')
print groups
thermochem = lib.Estimate(groups,'thermochem')
print thermochem.eval_ND_H(298.15)

defaultdict(<type 'int'>, {'C(C)(H)3': 2})
-34.4280812417
```
