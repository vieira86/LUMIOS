def comparation2(df1, df2): 
    import pandas as pd
    
    df2.reset_index(inplace=True)
        
    dereplication_H = [] # lista vazia que será preenchida com massas semelhantes à molécula
    for l, p in enumerate(df2['exact_mass'].values):
        r = df2['index'][l]

        a = (p-0.01)
        b = (p+0.01)
        #print(a, b)

        for x, j in enumerate (df1['parent_mass']):
            k = df1['index'][x]#buscar moléculas que estejam no nosso conjunto de dados

            if (j >= a and j <= b):
                #print(a, b)

                dereplication_H.append([j,k,l, r, p])

    resultado_H = pd.DataFrame(dereplication_H, columns = ['parent_mass_cacau', 'espectro_cacau', 'indice_base_atual', 'espectro_base', 'exact_mass' ])
    return resultado_H