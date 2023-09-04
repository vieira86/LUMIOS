def similaridade2(df):
    import pandas as pd
    import numpy as np

    # novo_array = []

    # for i in range(len(df)):
    #     a = df['peaks'][i]
    #     new_array = []
    #     for i in a:
    #         if i != 0:
    #             new_array.append(i)
                
    #     novo_array.append(new_array)

    # df['peaks'] = novo_array


    # novo_array2 = []
    # for i in range(len(df)):
    #     a = df['peaks_base'][i]
    #     new_array2 = []
    #     for i in a:
    #         if i != 0:
    #             new_array2.append(i)
                
    #     novo_array2.append(new_array2)

    # df['peaks_base'] = novo_array2


    def completando_array (a,b):
        import math
        while len(a) != len(b):
            if len(a) < len(b):
                a = np.append(a, 0)
            if len(a) > len(b):
                b = np.append(b, 0)
        return a, b

    def cosine_similarity(x, y):
        from numpy.linalg import norm
        return np.dot(x, y) / (norm(x) * norm(y))

    results = []
    for i in range(len(df)):
        a = df['peaks'][i]
        b = df['peaks_base'][i]
        
        a, b = completando_array(a,b)
        similarity = cosine_similarity(a,b)
        if similarity > 1:
            similarity2 = cosine_similarity(b,a)
            results.append(similarity2)
        else:
            results.append(similarity)

    df['similaridade'] = results

    resultado = df.loc[df['similaridade']>0.5]

    resultado = resultado.sort_values(by=['index_base', 'similaridade']).drop_duplicates(subset = ['index_amostra'], keep = 'last')
   
    # resultado = resultado.sort_values(by=['index_amostra', 'similaridade']).drop_duplicates(subset = ['index_amostra'], keep = 'first')
    # resultado = resultado.drop_duplicates(subset = 'index_base', keep = 'last')

    anotacoes = resultado.drop_duplicates(subset = 'smiles')
    # anotacoes = resultado
    erro =[]
    for i in anotacoes['parent_mass']:
        a = float(i)
        erro.append(a)
    anotacoes['parent_mass'] = erro

    erro =[]
    for i in anotacoes['exact_mass']:
        a = float(i)
        erro.append(a)
    anotacoes['exact_mass'] = erro

    anotacoes['erro'] = abs(((anotacoes['parent_mass'] - anotacoes['exact_mass']) / anotacoes['exact_mass']) *1000000)

    return anotacoes