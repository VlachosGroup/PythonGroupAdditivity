Python Group Additivity (pgradd)
================================
  
A Python package and database, developed by the Vlachos Research Group at the University of Delaware implements the **F**\ irst-**P**\ rinciples **S**\ emi-**E**\ mpirical (FPSE) **G**\ roup **A**\ dditivity
(GA) method for estimating thermodynamic properties of molecules. First introduced by Benson et al. for gas molecules and
was later extended by Kua et al. to species adsorbed on catalytic surfaces. GA relies on graph theory defining each molecule
as a collection of groups and their frequency of occurrence. The values of GA groups are determined from DFT-calculated
thermodynamic properties of a (training) set of molecules by linear regression to minimize the difference of thermodynamic
properties of molecules predicted by the GA from those estimated via DFT. This package implements four group additivity
schemes in six databases (See below) and will convert a molecule entered as a **S**\ implified **M**\ olecular-**I**\ nput
**L**\ ine-**E**\ ntry **S**\ ystem (`SMILES`_) providing the constituent groups, their frequency of occurrence, and estimated
thermodynamic properties for that molecule. pgradd also provides a general GA framework for implementing a custom group additivity scheme from your *ab initio*\  data and regression to groups.

.. image:: https://github.com/VlachosGroup/PythonGroupAdditivity/blob/master/docs/pGrAdd_RGB_github.png
   :target: https://vlachosgroup.github.io/PythonGroupAdditivity/
   :align: center
 
-  Benson's gas molecule group additivity (BensonGA)
-  Salciccioli et al. (2012) adsorbate on Pt(111) group additivity scheme (SalciccioliGA2012)
-  Gu et al. (2017) solvated adsorbate on Pt(111) group additivity scheme (GuSolventGA2017Aq, GuSolventGA2017Vac)
-  Wittreich (2018) adsorbate on Pt(111). Subset of Gu et al. including only surface species, group values regressed with OLS/GLS (Maximum Likelihood) and DFT data processed with `pmutt`_ (GRWSurface2018)
-  Wittreich (2018) solvated adsorbate on Pt(111). Subset of Gu et al. including only surface species, group values regressed with OLS/GLS (Maximum Likelihood) and DFT data processed with `pmutt`_ (GRWAqueous2018)
-  Xie (2022) Database for hydrocarbon species on Ru(0001) (XieGA2022)

Citing this work
----------------
G.R. Wittreich, D.G. Vlachos, Python Group Additivity (pGrAdd) software for estimating species thermochemical properties Comput. Phys. Commun. 273 (2022) 108277 https://doi.org/10.1016/j.cpc.2021.108277

Developers
----------

-  Gerhard R Wittreich, P.E.
-  Geun Ho Gu
-  Michael Salciccioli
-  Stephen M. Edie

Required Packages
-----------------

-  Python2/Python3
-  `pmutt`_ >= 1.3.2
-  `rdkit`_ >= 2018.03.4.0
-  ipython >= 7.3.1
-  `numpy`_ >= 1.21.5
-  `pyyaml`_ >= 6.0
-  `scipy`_ >= 1.9.3

Getting Started
---------------

1. Install using pip::

    pip install --user pgradd

2. Run the unit tests. Navigate to the **tests**\  directory, input the command shown below, and look for an **OK**\  response. (**Note:**\  The number of tests/time may change with subsequent versions)::

    python -m unittest
    
    ........................................
    ----------------------------------------------------------------------
    Ran 40 tests in 7.849s

    OK

3. Look at examples below

License
-------

This project is licensed under the MIT License - see the `LICENSE`_ file for details.

Contributing
------------

If you have a suggestion or find a bug, please post to our `Issues page`_ with 
the |enhancement_label| or |bug_label| tag respectively.

Finally, if you would like to add to the body of code, please:

- fork the development branch
- make the desired changes
- write the appropriate unit tests
- submit a `pull request`_.

Questions
---------

If you are having issues, please post to our `Issues page`_ with the ``help wanted`` or ``question`` tag. We 
will do our best to assist.

Special Thanks
--------------

-  Dr. Jeffrey Frey (pip and conda compatibility)

Citations
---------

-  Rangarajan et al. "Language-oriented rule-based reaction network generation and analysis: Algorithms of RING", Comput. Chem. Eng. 2014, 64, 124. https://doi.org/10.1016/j.compchemeng.2014.02.007
-  Rangarajan et al. "Language-oriented rule-based reaction network generation and analysis: Descrpition of RING", Comput. Chem. Eng. 2012, 45, 114. https://doi.org/10.1016/j.compchemeng.2012.06.008
-  Benson et al. "Additivity rules for the estimation of thermochemical properties." Chem. Rev., 1969, 69 (3), 279-324
-  Salciccioli et al. "Density Functional Theory-Derived Group Additivity and Linear Scaling Methods for Prediction of Oxygenate Stability on Metal Catalysts: Adsorption of Open-Ring Alcohol and Polyol Dehydrogenation Intermediates on Pt-Based Metals." J. Phys. Chem. C, 2010, 114 (47) 20155-20166. https://doi.org/10.1021/jp107836a
-  Kua J, Goddard WA (1998) Chemisorption of Organics on Platinum. 2. Chemisorption of C 2 H x and CH x on Pt(111). J Phys Chem B 102:9492–9500. https://doi.org/10.1021/jp982527s
-  Kua J, Faglioni F, Goddard WA (2000) Thermochemistry for hydrocarbon intermediates chemisorbed on metal surfaces: CH(n-m)(CH3)(m) with n = 1, 2, 3 and m ≤ n on Pt, Ir, Os, Pd, Rh, and Ru. J Am Chem Soc 122:2309–2321. https://doi.org/10.1021/ja993336l
-  Salciccioli et al. "Adsorption of Acid, Ester, and Ether Functional Groups on Pt: Fast Prediction of Thermochemical Properties of Adsorbed Oxygenates via DFT-Based Group Additivity Methods." J. Phys. Chem. C, 2012, 116(2), 1873-1886. https://doi.org/10.1021/jp2091413
-  Vorotnikov et al. "Group Additivity for Estimating Thermochemical Properties of Furanic Compounds on Pd(111)." Ind. Eng. Chem. Res., 2014, 53 (30), 11929-11938. https://doi.org/10.1021/ie502049a
-  Vorotnikov et al. "Group Additivity and Modified Linear Scaling Relations for Estimating Surface Thermochemistry on Transition Metal Surfaces: Application to Furanics." J. Phys. Chem. C, 2015, 119 (19), 10417-10426. https://doi.org/10.1021/acs.jpcc.5b01696
-  Gu et al. "Group Additivity for Thermochemical Property Estimation of Lignin Monomers on Pt(111)." J. Phys. Chem. C, 2016, 120 (34), 19234-19241. https://doi.org/10.1021/acs.jpcc.6b06430
-  Gu GH, Schweitzer B, Michel C, et al (2017) Group additivity for aqueous phase thermochemical properties of alcohols on Pt(111). J Phys Chem C 121:21510–21519. https://doi.org/10.1021/acs.jpcc.7b07340
-  Xie, T.; Wittreich, G. R.; Vlachos, D. G. Multiscale Modeling of Hydrogenolysis of Ethane and Propane on Ru(0001): Implications for Plastics Recycling. Appl. Catal. B Environ. 2022, 316 (June), 121597. https://doi.org/10.1016/j.apcatb.2022.121597

Examples
--------

**Benson's Gas Group Additivity Example**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('BensonGA')
    descriptors = lib.GetDescriptors('C1CO1')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_HoRT(298.15))

    Out:
    defaultdict(int, {'C(C)(H)2(O)': 2, 'O(C)2': 1, 'Oxirane': 1})
    -21.09467743150278


**Salciccioli et al. J. Phys. Chem. C, 2012, 116 (2), pp 1873-1886 Example**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('SalciccioliGA2012')
    descriptors = lib.GetDescriptors('C([Pt])C[Pt]')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_H(298.15, units='kcal/mol'))

    Out:
    defaultdict(<class 'int'>, {'C(C)(H)2(Pt)': 2, 'surface-ring strain': 0.217})
    -11.307743997749277

**Gu et al. J. Phys. Chem. C, 2017, 121 pp 21510–21519 Example**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('GuSolventGA2017Aq')
    descriptors = lib.GetDescriptors('C(=O)([Pt])O')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_HoRT(500))

    Out:
    defaultdict(<class 'int'>, {'CO(O)(Pt)+O(CO)(H)': 1.0})
    -109.86212002776878


**Wittreich Surface Example**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('GRWSurface2018')
    descriptors = lib.GetDescriptors('[Pt]C([Pt])C([Pt])([Pt])C=O')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_HoRT(750), '[Dimensionless]')
    print(thermochem.get_H(750, 'kcal/mol'), '[kcal/mol]')

    Out:
    defaultdict(<class 'int'>, {'C(C)(H)(Pt)2': 1, 'C(C)(CO)(Pt)2': 1, 'CO(C)(H)': 1,
                                'CPt2CPt2': 1, 'CCPt2': 1, 'surface-ring strain': 0.392})
    -13.423119203382337 [Dimensionless]
    -20.005853103142883 [kcal/mol]

**Wittreich Solvated Surface Example**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('GRWAqueous2018')
    descriptors = lib.GetDescriptors('C(=O)([Pt])O')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_HoRT(500), '[Dimensionless]')
    print(thermochem.get_H(500, 'kJ/mol'), '[kJ/mol]')

    Out:
    defaultdict(<class 'int'>, {'CO(O)(Pt)+O(CO)(H)': 1.0})
    -107.57909464133714 [Dimensionless]
    -447.23102885789655 [kJ/mol]
    
**Xie Ru(0001) Surface Example 1**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('XieGA2022')
    descriptors = lib.GetDescriptors('[Ru]C([Ru])C([Ru])([Ru])C')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_HoRT(500), '[Dimensionless]')
    print(thermochem.get_H(500, 'kJ/mol'), '[kJ/mol]')

    Out:
    defaultdict(<class 'int'>, {'C(C)(H)(Ru)2': 1, 'C(C)2(Ru)2': 1, 'C(C)(H)3': 1, 'CRu2CRu2': 1})
     -35.040312149773726 [Dimensionless]
    -145.6706333743726   [kJ/mol]
    
**Xie Ru(0001) Surface Example 2**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('XieGA2022')
    descriptors = lib.GetDescriptors('CCC')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_HoRT(500), '[Dimensionless]')
    print(thermochem.get_H(500, 'kJ/mol'), '[kJ/mol]')

    Out:
    defaultdict(<class 'int'>, {'C(C)(H)3': 2, 'C(C)2(H)2': 1})
    -41.49969417868688 [Dimensionless]
   -172.52376948049303 [kJ/mol]

**Free Energy of Formation by including Entropy of the Elements**::

    In:
    from pgradd.GroupAdd.Library import GroupLibrary
    import pgradd.ThermoChem
    lib = GroupLibrary.Load('BensonGA')
    descriptors = lib.GetDescriptors('CCCCCC')
    print(descriptors)
    thermochem = lib.Estimate(descriptors,'thermochem')
    print(thermochem.get_GoRT(T=298.15, S_elements=True), '[Dimensionless]')
    print(thermochem.get_G(T=298.15, units='kJ/mol', S_elements=True), '[kJ/mol]')

    Out:
    defaultdict(<class 'int'>, {'C(C)(H)3': 2, 'C(C)2(H)2': 4})
    -3.1192349163716244 [Dimensionless]
    -7.732446702038452  [kJ/mol]

.. |bug_label| image:: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/images/labels/bug_small.png
   :height: 20
   :target: https://github.com/VlachosGroup/PythonGroupAdditivity/issues?utf8=%E2%9C%93&q=label%3Abug

.. |enhancement_label| image:: https://raw.githubusercontent.com/VlachosGroup/pMuTT/master/docs/source/images/labels/enhancement_small.png
   :height: 20
   :target: https://github.com/VlachosGroup/PythonGroupAdditivity/issues?utf8=%E2%9C%93&q=label%3Aenhancement

.. _`scipy`: https://www.scipy.org/
.. _`rdkit`: https://www.rdkit.org/
.. _`numpy`: http://www.numpy.org/
.. _`pyyaml`: https://pyyaml.org/
.. _`SMILES`: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
.. _`pmutt`: https://github.com/VlachosGroup/pMuTT
.. _`LICENSE`: https://github.com/VlachosGroup/VlachosGroupAdditivity/blob/master/LICENSE.md
.. _`Issues page`: https://github.com/VlachosGroup/VlachosGroupAdditivity/issues
.. _`pull request`: https://github.com/VlachosGroup/VlachosGroupAdditivity/pulls
