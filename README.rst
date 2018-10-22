========================
Vlachos Group Additivity
========================
Python package for group additivity scheme use. See Wiki page for more info. Below are implemented group additivity scheme.

- Benson's gas molecule group additivity (BensonGA)
- Salciccioli et al. (2012) adsorbate on Pt(111) group additivity scheme (SalciccioliGA2012)
- Gu et al. (2017) solvated adsorbate on Pt(111) group additivity scheme (GuSolventGA2017Aq, GuSolventGA2017Vac)
- Wittreich (2018) adsorbate on Pt(111). Subset of Gu et al. including only surface species, group values regressed with OLS/GLS (Maximum Liklihood) and DFT data processed with pMuTT (GRWSurface2018)

Developers
----------
- Geun Ho Gu

Maintainers
-----------
- Gerhard R Wittreich, P.E. (wittregr@udel.edu)

Required Packages
-----------------
- Python2/Python3
- `rdkit`_ >= 2018.03.4.0
- ipython >= 7.0.0
- `numpy`_ >= 1.15.1
- `pyyaml`_ >= 3.0
- `scipy`_ >= 1.1.0

Citations
---------
- Rangarajan et al. "Language-oriented rule-based reaction network generation and analysis: Algorithms of RING", Comput. Chem. Eng. 2014, 64, 124
- Rangarajan et al. "Language-oriented rule-based reaction network generation and analysis: Descrpition of RING", Comput. Chem. Eng. 2012, 45, 114
- Benson et al. "Additivity rules for the estimation of thermochemical properties." Chem. Rev., 1969, 69 (3), 279-324
- Salciccioli et al. "Density Functional Theory-Derived Group Additivity and Linear Scaling Methods for Prediction of Oxygenate Stability on Metal Catalysts: Adsorption of Open-Ring Alcohol and Polyol Dehydrogenation Intermediates on Pt-Based Metals." J. Phys. Chem. C, 2010, 114 (47) 20155-20166
- Salciccioli et al. "Adsorption of Acid, Ester, and Ether Functional Groups on Pt: Fast Prediction of Thermochemical Properties of Adsorbed Oxygenates via DFT-Based Group Additivity Methods." J. Phys. Chem. C, 2012, 116(2), 1873-1886
- Vorotnikov et al. "Group Additivity for Estimating Thermochemical Properties of Furanic Compounds on Pd(111)." Ind. Eng. Chem. Res., 2014, 53 (30), 11929-11938
- Vorotnikov et al. "Group Additivity and Modified Linear Scaling Relations for Estimating Surface Thermochemistry on Transition Metal Surfaces: Application to Furanics." J. Phys. Chem. C, 2015, 119 (19), 10417-10426
- Gu et al. "Group Additivity for Thermochemical Property Estimation of Lignin Monomers on Pt(111)." J. Phys. Chem. C, 2016, 120 (34), 19234-19241
- Gu et al. "Group Additivity for Aqueous Phase Thermochemical Properties of Alcohols on Pt(111)." J. Phys. Chem. C, submitted

Examples
--------

**Benson's Gas Group Additivity Example::**

    In:
    from VGA.GroupAdd.Library import GroupLibrary
    import VGA.ThermoChem
    lib = GroupLibrary.Load('BensonGA')
    descriptors = lib.GetDescriptors('C1CO1')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.eval_ND_H(298.15))

    Out:
    defaultdict(int, {'C(C)(H)2(O)': 2, 'O(C)2': 1, 'Oxirane': 1})
    -21.09467743150278

**Salciccioli et al. J. Phys. Chem. C, 2012, 116 (2), pp 1873–1886 Example::**

    In:
    from VGA.GroupAdd.Library import GroupLibrary
    import VGA.ThermoChem
    lib = GroupLibrary.Load('SalciccioliGA2012')
    descriptors = lib.GetDescriptors('C([Pt])C[Pt]')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.eval_ND_H(298.15))

    Out:
    defaultdict(<class 'int'>, {'C(C)(H)2(Pt)': 2, 'surface-ring strain': 0.217})
    37.62494617247582

**Gu et al. J. Phys. Chem. C, submitted Example** 
::

    In:
    from VGA.GroupAdd.Library import GroupLibrary
    import VGA.ThermoChem
    lib = GroupLibrary.Load('GuSolventGA2017Aq')
    descriptors = lib.GetDescriptors('C(=O)([Pt])O')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.eval_ND_H(500))

    Out:
    defaultdict(<class 'int'>, {'CO(O)(Pt)+O(CO)(H)': 1.0})
    -109.86212002776878

**Wittreich Example** 
 ::

    In:
    from VGA.GroupAdd.Library import GroupLibrary
    import VGA.ThermoChem
    lib = GroupLibrary.Load('GRWSurface2018')
    descriptors = lib.GetDescriptors('[Pt]C([Pt])C([Pt])([Pt])C=O')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.eval_ND_H(750))

    Out:
    defaultdict(<class 'int'>, {'C(C)(H)(Pt)2': 1, 'C(C)(CO)(Pt)2': 1, 'CO(C)(H)': 1,
                                'CPt2CPt2': 1, 'CCPt2': 1, 'surface-ring strain': 0.392})
    -13.42320778481884


.. _scipy: https://www.scipy.org/
.. _rdkit: https://www.rdkit.org/
.. _numpy: http://www.numpy.org/
.. _pyyaml: https://pyyaml.org/
