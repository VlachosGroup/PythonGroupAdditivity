"""
=========================================
Defenition to draw RDKIT mol object (:mod:`pgradd.DrawMol`)
=========================================

Coverts a rdkit mol object to a svg image and display.

"""

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG, display


# http://rdkit.blogspot.com/2015/02/new-drawing-code.html
def moltosvg(mol, highlight=[], molSize=(400, 400), kekulize=True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except Exception:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)

    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])

    # Atom Label
    opts = drawer.drawOptions()
    # Atom name and index
    for i in range(mol.GetNumAtoms()):
        opts.atomLabels[i] = mol.GetAtomWithIdx(i).GetSymbol()+str(i)
    # radicals and charges
    for atom in mol.GetAtoms():
        nr = atom.GetNumRadicalElectrons()
        nc = atom.GetFormalCharge()
        if nr > 0:
            string = atom.GetSymbol() + ':'*divmod(nr, 2)[0] +\
                '.'*divmod(nr, 2)[1]
            opts.atomLabels[atom.GetIdx()] += string
        elif nc == 1:
            string = atom.GetSymbol() + '+'
            opts.atomLabels[atom.GetIdx()] += string
        elif nc > 1:
            string = atom.GetSymbol() + '+' + str(nc)
            opts.atomLabels[atom.GetIdx()] += string
        elif nc == -1:
            string = atom.GetSymbol() + '-'
            opts.atomLabels[atom.GetIdx()] += string
        elif nc < -1:
            string = atom.GetSymbol() + '-' + str(nc)
            opts.atomLabels[atom.GetIdx()] += string

    # highlight
    if highlight:
        drawer.DrawMolecule(mc, highlightAtoms=highlight)
    else:
        drawer.DrawMolecule(mc)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    svg.replace('svg:', '')
    display(SVG(svg))
