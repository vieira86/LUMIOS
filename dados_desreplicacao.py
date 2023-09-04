def dados_desreplicacao(df4, cacau_df, base):
    import pandas as pd    
    mask_amostra = df4['espectro_cacau'].values
    anotacoes_principais = cacau_df.loc[mask_amostra]
    anotacoes_principais['espectro_base'] = df4['espectro_base'].values
    anotacoes_principais.rename(columns = {'index': 'index_amostra', 'espectro_base': 'index_base'}, inplace=True)
    base.rename(columns ={'index': 'index_base'}, inplace=True)
    df = anotacoes_principais.merge(base, on=['index_base'], how='outer', suffixes=['', '_base'], indicator=True)
    both = df.loc[df['_merge']== 'both']
    df2=both
    df2['index_amostra'] = df2['index_amostra'].astype(int)
    df2['parent_mass'] = df2['parent_mass'].astype(float)
    #df2['index_base'] = int(float(df2['index_base']))
    df2['parent_mass_base'] = df2['parent_mass_base'].astype(float)
    df2.drop('_merge', axis=1, inplace=True)
    
    for col in df2.columns:
        df2[col] = df2[col].apply(lambda x: '{:.4f}'.format(x) if type(x) is float else x)
        
    return df2
