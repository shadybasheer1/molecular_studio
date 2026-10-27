from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

examples = {
    "Water": "O",
    "Methane": "C",
    "Ethane": "CC",
    "Ethanol": "CCO",
    "Acetone": "CC(=O)C",
    "Acetic acid": "CC(=O)O",
    "Benzene": "c1ccccc1",
    "Toluene": "Cc1ccccc1",
    "Phenol": "Oc1ccccc1",
    "Aniline": "Nc1ccccc1",
    "Glycine": "NCC(=O)O",
    "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    "Lactic acid": "CC(O)C(=O)O",
}
for name, smiles in examples.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("FAIL:", name, smiles)
    else:
        print("OK:", name, rdMolDescriptors.CalcMolFormula(mol), Chem.MolToSmiles(mol))
