def carregar_amostras(cacau_npy):  
    import pandas as pd
    import numpy as np
    # cacau_npy = np.load(r"C:\Users\PC\Desktop\cacau_inicio_npy\df_cacau_inicio.npy", allow_pickle=True)

    cacau_df = pd.DataFrame(cacau_npy)

    colunas = ['id', 'peaks', 'precursor_type', 'modo', 'precursor_mz', 'intensity', 'retention_time', 'parent_mass']
    for i, j in enumerate(colunas):
#         print(i, j)

        cacau_df.rename(columns = {i: j}, inplace=True)


    cacau_df.reset_index(inplace=True)
    return cacau_df