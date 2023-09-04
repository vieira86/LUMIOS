def similaridade(df):
    import pandas as pd
    
    parent_mass =[]
    for i in df['parent_mass']:
        x1 = float(i)
        parent_mass.append(x1)
        
    df['parent_mass'] = parent_mass
    
    exact_mass =[]
    for i in df['exact_mass']:
        x2 = float(i)
        exact_mass.append(x2)
        
    df['exact_mass'] = exact_mass
    
    commom = []
    for i in range(len(df)):
        juncao = (set(df['peaks'][i]) & set(df['peaks_base'][i]))

        #juncao.remove(0.)
        commom.append(juncao)

    combined = []
    for i in range(len(df)):
        combinar = (set(df['peaks'][i]) | set(df['peaks_base'][i]))

        #juncao.remove(0.)
        combined.append(combinar)

    comprimento_comum = []
    for i in commom:
        comprimento_comum.append(len(i))

    comprimento_combinado = []
    for i in combined:
        comprimento_combinado.append((len(i)))

    import numpy as np
    comprimento_comum_array = np.array(comprimento_comum) 
    comprimento_combinado_array = np.array(comprimento_combinado)


    similaridade = []

    for i in range(len(comprimento_comum_array)):
        a = (comprimento_comum_array[i]/comprimento_combinado_array[i])*100
        similaridade.append(a)

    df['similaridade'] = similaridade

    resultado = df.loc[df['similaridade']>10]

    # resultado = resultado.sort_values(by=['index_base', 'similaridade']).drop_duplicates(subset = ['index_amostra'], keep = 'last')
    resultado = resultado.sort_values(by=['index_amostra', 'similaridade']).drop_duplicates(subset = 'index_amostra', keep = 'last') #do MENOR PARA MAIOR

    #este "resultado" é a anotacao para todos os pontos.
    # resultado = resultado.drop_duplicates(subset = 'index_amostra', keep = 'last')

    #esta "anotacoes" já sao para ANOTACOES diferentes... retiro as duplicatas de smiles
    anotacoes = resultado.drop_duplicates(subset = 'smiles')

    anotacoes['erro'] = abs(((anotacoes['parent_mass'] - anotacoes['exact_mass']) / anotacoes['exact_mass']) *1000000)

    return resultado, anotacoes