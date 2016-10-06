from msr.GroupAdd.Scheme import GroupAdditivityScheme
from msr.GroupAdd.FindResStruc import enumerate_res_struc
from rdkit import Chem
from collections import defaultdict
import yaml
molecules = yaml.load(file('C:\\Users\\Gu\\Documents\\PythonScripts\\GA_input\\solvation_GA\\Solvatino_GA_smiles.yaml','r'))
# load scheme
scheme = GroupAdditivityScheme('solvation2')
total_groups = defaultdict(int)
fid = open('solvation_groups.dat', 'w')
for mol in molecules:
    #print mol
    # get resonance structure
    groups = scheme.GetGroups(mol)
    line = mol + '\t'
    for group in sorted(groups):
        total_groups[group] += 1
        line += '\t%s\t%d'%(group,groups[group])
    print line
    fid.write(line+'\n')
fid.close()
#for group in sorted(total_groups):
    #print '%s\t%d'%(group,total_groups[group])
