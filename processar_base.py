def processamento_base(base):
    import pandas as pd
    
    base = pd.DataFrame(base)
    
    colunas = ['idx_base_original', 'peaks',
        'precursor_type',
        'modo',
        'precursor_mz',
        'formula',
        'inchi',
        'smiles',
        'compound_name',
        'parent_mass', 'base', 'exact_mass']

    for i, j in enumerate(colunas):
            # print(i, j)
            #a = j
        base.rename(columns = {i: j}, inplace=True)

    # from rdkit.ML.Descriptors import MoleculeDescriptors
    # import rdkit
    # from rdkit import Chem
    # from itertools import combinations
    # from rdkit import Chem
    # from rdkit.Chem import Draw
    # from rdkit.Chem import Descriptors
    # from rdkit.ML.Descriptors import MoleculeDescriptors
    # def calc_descriptors(df):
    #     df = df.copy()  # para não dar o SettingWithCopyWarning
    #     descriptor_list = ["ExactMolWt"]  #  molecular surface area
    #     for desc in descriptor_list:
    #         names = [desc]
    #         calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
    #         df[desc]  = df["smiles"].apply(lambda x: calc.CalcDescriptors(Chem.MolFromSmiles(x))[0])    

    #     return df

    # df = calc_descriptors(base)
    # base['exact_mass'] = df['ExactMolWt']
    
    return base