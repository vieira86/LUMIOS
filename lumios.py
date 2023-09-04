#!/usr/bin/env python
# coding: utf-8
# ######################
# Import libraries
######################
# 
import streamlit as st
import streamlit.components.v1 as stc
import numpy as np
# File Processing Pkgs
import pandas as pd
import numpy as np
import pandas as pd
import streamlit as st
import pickle
from PIL import Image
import matchms
import os
from matchms.importing import load_from_msp
from matchms.filtering import default_filters
from matchms.filtering import repair_inchi_inchikey_smiles
from matchms.filtering import derive_inchikey_from_inchi
from matchms.filtering import derive_smiles_from_inchi
from matchms.filtering import derive_inchi_from_smiles
from matchms.filtering import harmonize_undefined_inchi
from matchms.filtering import harmonize_undefined_inchikey
from matchms.filtering import harmonize_undefined_smiles
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import select_by_intensity
from matchms.filtering import select_by_mz
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import combinations
import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from PIL import Image 


image = Image.open('capa_LUMIOS.png')

# st.sidebar.image(image, use_column_width=True)
# st.write("""
#         #            LUMIOS
#         ###### Label Using Machine In Organic Samples
#         ***
#         """)
# st.markdown("<h1 style='text-align: center; color: black;'>LUMIOS</h1>", unsafe_allow_html=True)
# st.markdown("<h5 style='text-align: center; color: gray;'>Label Using Machine In Organic Samples</h5>", unsafe_allow_html=True)



# image = Image.open('lumios.jpg')

# st.image(image, use_column_width=True)


def main():
    # st.title("Molecular Desreplication - APP")
    menu = ["Home", "Data Processing", "Dereplicator", "Machine Learning / Deep Learning",  "Molecular Docking", "Storytelling with data", "About"]
    # image = Image.open('Lumios_aba3.png')
    # image2 = Image.open('pipeline_lumios.png')
    # st.sidebar.image(image, use_column_width=True)
    choice = st.sidebar.selectbox("Menu", menu)
    if choice == "Home":
        
        from PIL import Image 
        image = Image.open('capa_LUMIOS.png')
        image2 = Image.open('Lumios_aba2.png')
        image3 = Image.open('Lumios_logo.png')

        st.image(image, use_column_width=True)
        st.sidebar.image(image2, use_column_width=True)
        # st.sidebar.image(image3, use_column_width=True)
        # st.sidebar.markdown("<h1 style='text-align: center; color: red;'>LUMIOS is the acronym for Label Using Machine In Organic Samples, and it is a multitasking software created \
        # from Python language libraries (adding already consolidated libraries) that aims to help professionals and students of \
        #     area of ​​organic chemistry that perform computational approaches for rational exploration of matrices \
        #         complex nature of natural products, suggesting molecular labels, known as annotations</h1>")
        # st.sidebar.write("""
        # #            LUMIOS
        # ###### Label Using Machine In Organic Samples
        # ***
        # """)
        st.markdown("<h6 style='text-align: justify; color: black;'><b>LUMIOS</b> is the acronym for Label Using Machine In Organic Samples, \
            and it is a multitasking software created \
         from Python language libraries (adding already consolidated libraries) that aims to help professionals and students of \
             area of ​​organic chemistry that perform computational approaches for rational exploration of matrices \
                complex nature of natural products, suggesting molecular labels, known as annotations.</h6>", unsafe_allow_html=True)
    

    elif choice == ("Dereplicator"):
        from PIL import Image 
        image = Image.open('Lumios_logo.png')
        # st.sidebar.image(image, use_column_width=True)
        # image = Image.open('Desreplicator.png')
        st.sidebar.image(image, use_column_width=True)
        st.write("## LUMIOS - Desreplicator")
        

        data_file = st.file_uploader("Upload CSV",type=['npy'], accept_multiple_files=False)
	    
        if st.button("Process"):
            if data_file is not None:
                
                import numpy as np
                cacau_npy = np.load(data_file, allow_pickle=True)

                import carregamento_amostras
                cacau_df = carregamento_amostras.carregar_amostras(cacau_npy)

                float_intensity = []
                for i in cacau_df['intensity']:
                    a = float(i)
                    float_intensity.append(a)

                cacau_df['intensity'] = float_intensity

                import plotly.express as px
# df = px.data.gapminder()

                fig = px.scatter(cacau_df, x="retention_time", y="parent_mass",
                                size="intensity", color="id",
                                hover_name="id", log_x=False, size_max=60)
                
                st.plotly_chart(fig, use_container_width=True)
                
                # st.dataframe(cacau_df)

				# ##################
				# # Abrindo as bases
				# ##################
                base = np.load('base_completa_exact_mass_correto.npy', allow_pickle=True)

                import processar_base
                base = processar_base.processamento_base(base)

                st.write('Comparing spectra with', len(base), 'samples in database')
                
                # import time
                # with st.spinner('Processing...'):
                #     time.sleep(5)
                #     st.success('Done!')
                    
                import desreplicador
                base.reset_index(inplace=True)


                # st.write('Desreplicando seus dados...')
                # with st.spinner('Processando os dados...'):
                #     time.sleep(5)
                #     st.success('Done!')

                df4 = desreplicador.comparation2(cacau_df, base)

                # st.write('Dereplicating seus dados...')
                # with st.spinner('Processando os dados...'):
                #     time.sleep(5)
                #     st.success('Done!')
                # df4 = pd.read_csv('df4.csv')

                import time
                with st.spinner('Processing...'):
                    time.sleep(5)
                    st.success('Done!')

                st.write('Matches was found...')
                st.dataframe(df4)
            
                import dados_desreplicacao
                df2 = dados_desreplicacao.dados_desreplicacao(df4, cacau_df, base )
                
                # import time
                # with st.spinner('Processando os dados...'):
                #     time.sleep(5)
                #     st.success('Done!')

                st.dataframe(df2)

                import similaridade
                resultado, anotacoes = similaridade.similaridade(df2)
                # st.write('## possívies anotacoes para todos os pontos')
                # st.dataframe(resultado)

                st.write('## Annotations')
                st.dataframe(anotacoes)

                import pandas as pd
                lotus = pd.read_csv('Lotus_inchi.csv')

                import funcao_produto_natural
                anotacao_final = funcao_produto_natural.produtos_naturais(lotus, anotacoes)

                st.write('### Annotations - LOTUS filter')
                st.dataframe(anotacao_final)

                st.write('### Annotations - No filter')
                st.dataframe(anotacoes)

                ### utilizar o calculo de cosseno para corroborar com a similaridade:

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
                    return (np.dot(x, y) / (norm(x) * norm(y)))

                results = []
                anotacao_final.reset_index(inplace=True, drop=True)

                # st.markdown('resetou??')
                # st.dataframe(anotacao_final)

                for i in range(len(anotacao_final)):
                    a = anotacao_final['peaks'][i]
                    b = anotacao_final['peaks_base'][i]
                    
                    a, b = completando_array(a,b)
                    similarity = cosine_similarity(a,b)
                    if (similarity > 1):
                        similarity2 = cosine_similarity(b,a)
                        results.append(similarity2)
                    else:
                        results.append(similarity)

                anotacao_final['cosine'] = results


                results = []
                anotacoes.reset_index(inplace=True, drop=True)

                # st.markdown('resetou??')
                # st.dataframe(anotacao_final)

                for i in range(len(anotacoes)):
                    a = anotacoes['peaks'][i]
                    b = anotacoes['peaks_base'][i]
                    
                    a, b = completando_array(a,b)
                    similarity = cosine_similarity(a,b)
                    if (similarity > 1):
                        similarity2 = cosine_similarity(b,a)
                        results.append(similarity2)
                    else:
                        results.append(similarity)

                anotacoes['cosine'] = results

                st.markdown("## Results - Lotus Filter")
                # st.write("você tem", len(anotacao_final), "anotacoes moleculares (Lotus Filter)")

                # st.write('### Annotations')
                st.dataframe(anotacao_final)

                st.markdown("## Results - No Filter")
                # st.write("você tem", len(anotacoes), "anotacoes moleculares (Sem Lotus Filter)")

                # st.write('### Anotacao final')
                st.dataframe(anotacoes)


                @st.cache
                def convert_df(df):
                    # IMPORTANT: Cache the conversion to prevent computation on every rerun
                    return df.to_csv().encode('utf-8')

                csv = convert_df(anotacao_final)

                st.download_button(
                    label="Download Results - With Lotus Filter",
                    data=csv,
                    file_name='anotations_filter.csv',
                    mime='text/csv',
                )

                csv = convert_df(anotacoes)

                st.download_button(
                    label="Download Results - Without Lotus Filter",
                    data=csv,
                    file_name='anotations_nofilter.csv',
                    mime='text/csv',
                )


                st.write("### Anottations - With Lotus")
                import streamlit.components.v1 as components
                anotacao_final.rename(columns={'smiles': 'SMILES'}, inplace=True)
                raw_html = mols2grid.display(anotacao_final,  subset=["compound_name", 'img'])._repr_html_()
                # components.html(raw_html, width=900, height=900, scrolling=True)
                components.html(raw_html)

                dir = r'C:\Users\PC\Desktop\video_LUMIOS\anotacoes'
                import os
                if not os.path.exists(dir):
                    os.makedirs(dir)

                np.save(r"C:\Users\PC\Desktop\video_LUMIOS\anotacoes", anotacao_final)

                st.write("### Anottations - Without Lotus")
                import streamlit.components.v1 as components
                anotacoes.rename(columns={'smiles': 'SMILES'}, inplace=True)
                raw_html = mols2grid.display(anotacoes,  subset=["compound_name", 'img'])._repr_html_()
                components.html(raw_html, width=900, height=900, scrolling=True)
                # components.html(raw_html)

                np.save(r"C:\Users\PC\Desktop\video_LUMIOS\anotacoes", anotacoes)

                ## Verificar similaridade de anotacao específica:

               
            else:
                st.write('Insert your data')
                
                
    elif choice == ("Data Processing"):
        st.write("""
        ## LUMIOS - Data Processing""")
        from PIL import Image 
        image = Image.open('Lumios_logo.png')
        st.sidebar.image(image, use_column_width=True)
        title = st.text_input('Directory for your data in .MSP format')

        if len(title)==0:
            # image = Image.open('Data_processing.png')

            st.sidebar.write()

        else:
            if st.button("Process"):
                if title is not None:
                    st.write('Data in:', title)

                    import processamento_inicial_amostras
                    amostras = processamento_inicial_amostras.processamento_amostras(title)
                    # df = cacau.values

                    # cacau.to_csv('df_testando.csv', index=False)

                    # np.save("df_cacau_completo_com_plotly.npy", cacau)
                    import time

                    with st.spinner('Processing...'):
                        time.sleep(5)
                        st.success('Done!')

                        st.dataframe(amostras)

                    import numpy as np
                    import os
                    salvar_np = os.path.join(title, "processed_samples.npy")
                    np.save(salvar_np, amostras)
                     

                    @st.cache
                    def convert_df(df):
                        # IMPORTANT: Cache the conversion to prevent computation on every rerun
                        return df.to_csv().encode('utf-8')

                    csv = convert_df(amostras)

                    st.download_button(
                        label="Download data as CSV",
                        data=csv,
                        file_name='processed_samples.csv',
                        mime='text/csv',
                    )

    elif choice == ("Machine Learning / Deep Learning"):
        from PIL import Image 
        image = Image.open('Lumios_logo.png')
        st.sidebar.image(image, use_column_width=True)
        st.write("""
        ## LUMIOS - Machine Learning""")
        
        # image = Image.open('machine.png')

        # st.sidebar.image(image, use_column_width=True)

                
        st.write("##### Insert annotations")
        data_file = st.file_uploader("Upload CSV",type=['csv'], accept_multiple_files=False)
        if st.button("Start machine learning"):
            if data_file is not None:
                import pandas as pd
                from PIL import Image 
                df = pd.read_csv(data_file)
                # df.reset_index(inplace=True)
                num_vars = ['VSA_EState6', 'NumAromaticRings', 'fr_NH0', 'SlogP_VSA10', 'SMR_VSA3',
       'BCUT2D_MWHI', 'fr_aniline', 'fr_Ar_N', 'fr_Al_OH', 'fr_halogen',
       'FractionCSP3', 'VSA_EState3', 'BalabanJ', 'NumAromaticCarbocycles',
       'fr_allylic_oxid', 'SlogP_VSA1', 'fr_Al_OH_noTert', 'fr_benzene',
       'PEOE_VSA4', 'BCUT2D_CHGHI', 'NumAromaticHeterocycles', 'VSA_EState8',
       'EState_VSA10', 'RingCount', 'PEOE_VSA1', 'BCUT2D_MRHI', 'EState_VSA1',
       'SMR_VSA4', 'fr_Ar_OH', 'EState_VSA4']
                def calc_descriptors(df, num_vars):
                    from rdkit import Chem    # make sure to import it if you haven't done so
                    from rdkit.Chem import Descriptors 
                    from rdkit.ML.Descriptors import MoleculeDescriptors
                    df = df.copy()  # para não dar o SettingWithCopyWarning
                    #     descriptor_list = ["FractionCSP3", "MolLogP", "MolWt", "NumAromaticRings", "NumHAcceptors", 
                    #                        "NumHDonors", "NumRotatableBonds", "TPSA", "fr_azo", "LabuteASA"]  #  molecular surface area
                    #     descriptors_list = [x[0] for x in Descriptors._descList]
                    descriptors_list = num_vars
                    for desc in descriptors_list:
                    #         descriptors_list = features_top30.index
                        names = [desc]
                        calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
                        df[desc]  = df["smiles"].apply(lambda x: calc.CalcDescriptors(Chem.MolFromSmiles(x))[0])
                    return df

                anotacoes_descrip = calc_descriptors(df, num_vars)


                # from sklearn.impute import SimpleImputer
                # from sklearn.preprocessing import RobustScaler
                # from sklearn.preprocessing import OneHotEncoder
                # from sklearn.preprocessing import StandardScaler

                # from sklearn.pipeline import Pipeline
                # from sklearn.compose import ColumnTransformer

                # num_pipeline = Pipeline([
                #     ('imputer', SimpleImputer(strategy='mean')),
                #     ('robust_scaler', StandardScaler())
                # ])


                # # (name, transformer, columns)
                # preprocessed_pipeline = ColumnTransformer([
                #     ('numerical', num_pipeline, num_vars)
                # ])

                # anotation_proc = preprocessed_pipeline.fit_transform(anotacoes_descrip)
                anotacoes_descrip = anotacoes_descrip[num_vars]
                # import lightgbm as lgb
                import pickle
                # load_saved_model = load_model('Final_Model')
                loaded_model = pickle.load(open('modelo_lgb_novo_sem_processar', 'rb'))
                predictions_anotacoes = loaded_model.predict(anotacoes_descrip)

                anotacoes_descrip['predict'] = predictions_anotacoes

                # st.dataframe(anotacoes_descrip)

                # st.dataframe(df)
                # df.markdown('Vendo como é o df')
                # df.rename(columns = {'Unnamed: 0': 'index1'}, inplace=True)
                # st.markdown('Vendo como é o df')
                
                # st.dataframe(df)
                df.drop('Unnamed: 0', axis=1, inplace=True)
                
                # st.markdown('Vendo como é o df resetado')
                df.reset_index(inplace=True, drop=False)
                
                st.dataframe(df)
                anotacoes_descrip.reset_index(inplace=True)
                # anotacoes_descrip.rename(columns = {'index': 'index1'}, inplace=True)
                # st.markdown('Vendo como é o anotacoes_descrip')
                st.dataframe(anotacoes_descrip)

                anotacoes_finais = df.merge(anotacoes_descrip, how='inner', on='index')
                # st.markdown('juntos')
                st.dataframe(anotacoes_finais)
            

                # st.write('## Predição dos dados')

                # drugs = anotacoes_descrip.loc[anotacoes_descrip['predict'] == 'Drug']

                # @st.cache
                # def convert_df(df):
                #     # IMPORTANT: Cache the conversion to prevent computation on every rerun
                #     return df.to_csv(index=False).encode('utf-8')

                # csv = convert_df(anotacoes_finais)

                # st.download_button(
                # label="Download data as CSV",
                # data=csv,
                # file_name='amostras_machine_learning.csv',
                # mime='text/csv',
                # )
                anotacoes_finais_graph = anotacoes_finais.copy()
                import streamlit.components.v1 as components
                anotacoes_finais_graph.rename(columns={'smiles': 'SMILES'}, inplace=True)
                raw_html = mols2grid.display(anotacoes_finais_graph, subset=["predict", 'img'])._repr_html_()
                # components.html(raw_html, width=900, height=900, scrolling=True)
                components.html(raw_html, scrolling=True)
                

                # st.dataframe(df)

                st.markdown('Neural Networking')
                ##CRIANDO A PASTA E CONVERTENDO PARA AS ESTRUTURAS    
                import os
                # from PIL import Image
                import pandas as pd 
                import numpy as np
                import matplotlib.pyplot as plt
                import seaborn as sns
                from rdkit.Chem import Draw
                from itertools import combinations
                from rdkit import Chem
                from rdkit.Chem import Draw
                from rdkit.Chem import Descriptors
                from rdkit.ML.Descriptors import MoleculeDescriptors
                # rdkit.Chem.Draw 
                from rdkit.Chem.Draw import IPythonConsole

                dir = r'C:\Users\PC\Desktop\Projeto_final\pasta_temporaria'
                import os
                if not os.path.exists(dir):
                    os.makedirs(dir)

                # os.mkdir(r'C:\Users\PC\Desktop\Projeto_final\pasta_temporaria')
                for i in range(0, len(df)):
                    ms = [Chem.MolFromSmiles(df['smiles'][i])]
                    img = Draw.MolsToGridImage(ms, molsPerRow=1, maxMols=1, subImgSize=(100, 100), returnPNG=False)
                    #img
                    #identif = plants['Species'][i]
                    #identif.astype(str)
                    nome = ('teste_'+ repr(i)+'.png')
                    
                    path = r'C:\Users\PC\Desktop\Projeto_final\pasta_temporaria'
                    endereco = os.path.join(path, nome)

                    img.save(endereco)

                import warnings
                warnings.filterwarnings('ignore')
                import tensorflow as tf
                import os

                nomes_amostras_teste = os.listdir(r'C:\Users\PC\Desktop\Projeto_final\pasta_temporaria')

                import cv2
                imagens_teste = []
                for i in nomes_amostras_teste:
                    path = r'C:\Users\PC\Desktop\Projeto_final\pasta_temporaria'
                    nome = i
                #     path = r'C:\Users\PC\Desktop\Redes_Neurais_LUMIOS\imagens\drugs'
                    endereco = os.path.join(path, nome)
                    print(endereco)
                    img = cv2.imread(endereco)
                    imagens_teste.append(img)

            
                import numpy as np
                x_novo_teste = np.zeros((len(imagens_teste), 100, 100, 3)) # (n_amostras, (dimensao))   

                defeito_teste = []
                for i in range(len(imagens_teste)):
                    try:
                        x_novo_teste[i] = imagens_teste[i]
                    except:
                        print(i)
                        defeito_teste.append(i)  

                for i in range(x_novo_teste.shape[0]):
                    x_novo_teste[i] = x_novo_teste[i]/250  

                x_novo_teste = np.delete(x_novo_teste, defeito_teste, axis=0)    

                from tensorflow import keras
                model = keras.models.load_model(r'C:\Users\PC\Desktop\LUMIOS\Projeto_final\modelo_rede_neural_lumios_final') 

                y_test_proba_test = model.predict(x_novo_teste)
                
                y_test_pred = []
                for i in y_test_proba_test:
                    if i >= 0.6:
                        a = 'PN'
                    else:
                        a= 'Drug'
                    y_test_pred.append(a)

                import cv2
                imagens_end = []
                for i in nomes_amostras_teste:
                    path = r'C:\Users\PC\Desktop\Projeto_final\pasta_temporaria'
                    nome = i
                #     path = r'C:\Users\PC\Desktop\Redes_Neurais_LUMIOS\imagens\drugs'
                    endereco = os.path.join(path, nome)
                    print(endereco)
                    # img = cv2.imread(endereco)
                    imagens_end.append(endereco)

                lista = pd.DataFrame(imagens_end)
                ap=[]
                import re
                for i in lista[0]:
                    i = re.sub('[^0-9]', '', i)
                    ap.append(i)
                dados = pd.DataFrame(ap)
                dados[0] = dados[0].astype(int)
                dados['pred'] = y_test_pred

                dat = dados.sort_values(by=[0])
                pred = []
                for i in dat['pred']:
                    pred.append(i)
                                

                anotacoes_finais['neural_predict']  = pred

                st.markdown('Download')
                
                @st.cache
                def convert_df(df):
                    # IMPORTANT: Cache the conversion to prevent computation on every rerun
                    return df.to_csv(index=False).encode('utf-8')
                # anotacoes_finais.drop('smiles', axis=1, inplace=True)
                # anotacoes_finais.rename(columns={'SMILES': 'smiles'}, inplace=True)
                csv = convert_df(anotacoes_finais)

                st.download_button(
                label="Download data as CSV",
                data=csv,
                file_name='samples_machine_learning.csv',
                mime='text/csv',
                )
                # anotacoes_finais_graph = anotacoes_finais.copy()
                import streamlit.components.v1 as components
                anotacoes_finais.rename(columns={'smiles': 'SMILES'}, inplace=True)
                raw_html = mols2grid.display(anotacoes_finais, subset=["neural_predict", 'img'] )._repr_html_()
                components.html(raw_html, width=900, height=900, scrolling=True)
                # anotacoes_finais.rename(columns={'SMILES': 'smiles'}, inplace=True)

                import shutil
                shutil.rmtree(r'C:\Users\PC\Desktop\Projeto_final\\pasta_temporaria')


            else:
                st.write('Insert data')



    elif choice == ("Storytelling with data"):
        st.write("""
        ## LUMIOS - Storytelling""")
        from PIL import Image 
        image = Image.open('Lumios_logo.png')
        st.sidebar.image(image, use_column_width=True)

        # st.sidebar.markdown("**Story**")


        image = Image.open('Storytelling.png')
        st.sidebar.image(image, use_column_width=True)
                
        st.write("##### Insert annotations")
        data_file = st.file_uploader("Upload CSV",type=['csv'], accept_multiple_files=False)
        if st.button("Start storytelling"):
            if data_file is not None:
                import pandas as pd
                df = pd.read_csv(data_file)
                
                # Seleção das principais colunas:
    #             df = df[['smiles','compound_name','exact_mass', 'similaridade', 'erro',
    #    'cosine',
    #    'parent_mass_base', 'VSA_EState6', 'NumAromaticRings', 'fr_NH0', 'SlogP_VSA10', 'SMR_VSA3',
    #    'BCUT2D_MWHI', 'fr_aniline', 'fr_Ar_N', 'fr_Al_OH', 'fr_halogen',
    #    'FractionCSP3', 'VSA_EState3', 'BalabanJ', 'NumAromaticCarbocycles',
    #    'fr_allylic_oxid', 'SlogP_VSA1', 'fr_Al_OH_noTert', 'fr_benzene',
    #    'PEOE_VSA4', 'BCUT2D_CHGHI', 'NumAromaticHeterocycles', 'VSA_EState8',
    #    'EState_VSA10', 'RingCount', 'PEOE_VSA1', 'BCUT2D_MRHI', 'EState_VSA1',
    #    'SMR_VSA4', 'fr_Ar_OH', 'EState_VSA4','predict', '7P2G', '4DD8', '1NC6', '6VVU']]

                df = df[['smiles','compound_name','exact_mass', 'similaridade', 'erro',
       'cosine',
       'parent_mass_base', 'VSA_EState6', 'NumAromaticRings', 'fr_NH0', 'SlogP_VSA10', 'SMR_VSA3',
       'BCUT2D_MWHI', 'fr_aniline', 'fr_Ar_N', 'fr_Al_OH', 'fr_halogen',
       'FractionCSP3', 'VSA_EState3', 'BalabanJ', 'NumAromaticCarbocycles',
       'fr_allylic_oxid', 'SlogP_VSA1', 'fr_Al_OH_noTert', 'fr_benzene',
       'PEOE_VSA4', 'BCUT2D_CHGHI', 'NumAromaticHeterocycles', 'VSA_EState8',
       'EState_VSA10', 'RingCount', 'PEOE_VSA1', 'BCUT2D_MRHI', 'EState_VSA1',
       'SMR_VSA4', 'fr_Ar_OH', 'EState_VSA4','predict', '7P2G', '4DD8', '1NC6', '6VVU']]

                from sklearn.preprocessing import StandardScaler, MinMaxScaler
                dfa = df[['VSA_EState6', 'NumAromaticRings', 'fr_NH0', 'SlogP_VSA10', 'SMR_VSA3',
       'BCUT2D_MWHI', 'fr_aniline', 'fr_Ar_N', 'fr_Al_OH', 'fr_halogen',
       'FractionCSP3', 'VSA_EState3', 'BalabanJ', 'NumAromaticCarbocycles',
       'fr_allylic_oxid', 'SlogP_VSA1', 'fr_Al_OH_noTert', 'fr_benzene',
       'PEOE_VSA4', 'BCUT2D_CHGHI', 'NumAromaticHeterocycles', 'VSA_EState8',
       'EState_VSA10', 'RingCount', 'PEOE_VSA1', 'BCUT2D_MRHI', 'EState_VSA1',
       'SMR_VSA4', 'fr_Ar_OH', 'EState_VSA4']]


                import plotly.express as px
                from sklearn.decomposition import PCA

                # df = px.data.iris()
                # X = df[['sepal_length', 'sepal_width', 'petal_length', 'petal_width']]


                st.markdown("### Principal Component Analysis - PCA 3D")

                pca = PCA(n_components=3)
                components = pca.fit_transform(dfa)

                total_var = pca.explained_variance_ratio_.sum() * 100

                fig = px.scatter_3d(
                    components, x=0, y=1, z=2, color=df['compound_name'],
                    title=f'Total Explained Variance: {total_var:.2f}%',
                    labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
                )
                st.plotly_chart(fig)

                st.markdown("#### Principal Component Analysis - PCA 2D")

                import plotly.express as px
                from sklearn.decomposition import PCA

                # df = px.data.iris()
                # X = df[['sepal_length', 'sepal_width', 'petal_length', 'petal_width']]

                # pca = PCA(n_components=3)
                components = pca.fit_transform(dfa)

                total_var = pca.explained_variance_ratio_.sum() * 100

                fig = px.scatter(
                    components, x=0, y=1, color=df['compound_name'],
                    title=f'Total Explained Variance: {total_var:.2f}%',
                    labels={'0': 'PC 1', '1': 'PC 2'}
                )
                st.plotly_chart(fig)

                st.markdown("#### Hierarchical Clustering Analysis - HCA")

                import plotly.figure_factory as ff
                import numpy as np
                np.random.seed(1)

                # X = np.random.rand(15, 12) # 15 samples, with 12 dimensions each
                fig = ff.create_dendrogram(dfa)
                fig.update_layout(width=800, height=500)
                st.plotly_chart(fig)

                st.markdown('#### Correlations of descriptors')
                import matplotlib.pyplot as plt
                f = plt.figure(figsize=(25, 20))
                plt.matshow(df.corr(), fignum=f.number)
                plt.xticks(range(df.select_dtypes(['number']).shape[1]), df.select_dtypes(['number']).columns, fontsize=14, rotation=90)
                plt.yticks(range(df.select_dtypes(['number']).shape[1]), df.select_dtypes(['number']).columns, fontsize=14)
                cb = plt.colorbar()
                cb.ax.tick_params(labelsize=14)
                # plt.title('Correlation Matrix', fontsize=20)
                st.pyplot(plt)


                st.markdown("#### Anotation and similarity value")
                ### Anotations and similarity

                anotation_similarity = df.sort_values(by='similaridade', ascending=False)
                anotation_similarity = anotation_similarity.drop_duplicates(subset = 'smiles', keep = "first")

                # plotando o mesmo gráfico, mas na horizontal para um melhor entendimento
                import seaborn as sns
                import pandas as pd
                import seaborn as sns
                import matplotlib.pyplot as plt
                import numpy as np

                # estilos padrão para os plots/visualizações
                sns.set_theme(style="whitegrid")
                # plt.rcParams['figure.figsize'] = (20, 20)
                # plt.rcParams['axes.labelsize'] = 14
                # plt.rcParams['xtick.labelsize'] = 14
                # plt.rcParams['xtick.labelsize'] = 14
                # plt.rcParams['ytick.labelsize'] = 14

                fig = plt.figure(figsize = (15, 10))
                sns.barplot(data=anotation_similarity, x='similaridade', y='compound_name')
                plt.xlabel('Similarity',  fontsize=18)
                plt.ylabel('Compound Name' , fontsize=18)
                plt.tick_params(axis='both', which='major', labelsize=18)
                # plt.xticks([0, 0.2e9, 0.4e9, 0.6e9, 0.8e9, 1e9], ['0', '200M', '400M', '600M', '800M', '1B'])
                # plt.title('Similarity of anotations')
                st.pyplot(fig)

                st.markdown('#### Classification using Machine Learning')
                fig = plt.figure(figsize = (20, 10))
                sns.countplot(data=df, x='predict')
                # plt.xticks()
                plt.xlabel('Predict',  fontsize=18)
                plt.ylabel('Quantity of molecules in prediction' , fontsize=18)
                plt.tick_params(axis='both', which='major', labelsize=18)
                # plt.title('Prediction Distribution of anotations', size = 20)
                st.pyplot(fig)

                st.markdown('### Data distribution - exact mass of anotations')
                fig = plt.figure(figsize = (20, 10))
                # sns.set(font_scale=14)
             # kde ==> kernel density estimation
# É um jeito de estimar a função densidade de probabilidade de uma variável aleatória.

                sns.histplot(data=df, x='exact_mass', kde=True)
                plt.xlabel('Exact Mass - anotations', fontsize=18)
                plt.ylabel('Quantity', fontsize=18)
                plt.tick_params(axis='both', which='major', labelsize=18)
                # plt.title('Distribution of Exact Mass for anotations', fontsize=20)
                st.pyplot(fig)

                st.markdown('#### Data distribution (Violin Plot) - exact mass of anotations')
                fig = plt.figure(figsize = (20, 10))
                # sns.set(font_scale=14)
                sns.violinplot(data=df, x='exact_mass', kde=True)
                plt.xlabel('Exact Mass - anotations', fontsize=18)
                plt.ylabel('Quantity', fontsize=18)
                plt.tick_params(axis='both', which='major', labelsize=18)
                # plt.title('Distribution of Exact Mass for anotations', fontsize=20)
                st.pyplot(fig)


                st.markdown('#### Probability density function for exact mass and similarity')
                fig = plt.figure(figsize = (20, 10))
                ax = sns.kdeplot(data=df, x='exact_mass', hue='predict')
                plt.xlabel('Exact Mass - anotations', fontsize=18)
                plt.ylabel('Quantity', fontsize=18)
                plt.tick_params(axis='both', which='major', labelsize=18)
                plt.setp(ax.get_legend().get_texts(), fontsize='20') # for legend text
                plt.setp(ax.get_legend().get_title(), fontsize='20') # for legend title
                # plt.legend(fontsize=15)
                st.pyplot(fig)

                st.markdown('#### Comparition similarity or cosine score with exact mass of anotations')
                # fig = plt.figure(figsize = (15, 10))
                # fig, axs = plt.subplots(2, 1, figsize=(10, 6))
                sns.lmplot(data=df, x="exact_mass", y="similaridade", aspect=9, height=2)
                plt.xlabel('Exact Mass - anotations', fontsize=18)
                plt.ylabel('Similarity', fontsize=18)
                plt.tick_params(axis='both', which='major', labelsize=18)
                # sns.lmplot(data=df, x="exact_mass", y="cosine", aspect=9, height=2)
                st.pyplot(plt)

                sns.lmplot(data=df, x="exact_mass", y="cosine", aspect=9, height=2)
                plt.xlabel('Exact Mass - anotations', fontsize=18)
                plt.ylabel('Cosine', fontsize=18)
                plt.tick_params(axis='both', which='major', labelsize=18)
                # sns.lmplot(data=df, x="exact_mass", y="cosine", aspect=9, height=2)
                st.pyplot(plt)

                # st.markdown('#### Influency sp3 carbon in similarity')
                # fig = plt.figure(figsize = (15, 10))
                # # fig, axs = plt.subplots(2, 1, figsize=(10, 6))
                # sns.jointplot(data=df, x='FractionCSP3', y='similaridade', kind='hex', height=10)
                # plt.xlabel('Sp3 Carbon - anotations', fontsize=18)
                # plt.ylabel('Similarity', fontsize=18)
                # plt.tick_params(axis='both', which='major', labelsize=18)
                # # sns.lmplot(data=df, x="exact_mass", y="cosine", aspect=9, height=2)
                # st.pyplot(plt)

                st.markdown('#### Similarity in two classes (Drugs and Natural Products')
                fig, axs = plt.subplots(1, 2, figsize=(10, 6))
                sns.histplot(data=df.query('predict == "Drug"'), x="similaridade", kde=True, color="blue", ax=axs[0])
                sns.histplot(data=df.query('predict != "Drug"'), x="similaridade", kde=True, color="black", ax=axs[1])
                # plt.xlabel('Similarity comparatitions')
                st.pyplot(fig)

                st.markdown('#### Docking results - 7P2G - SARS-Cov-2')

                from plotly.subplots import make_subplots
                import plotly.graph_objects as go
                df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]

                import plotly.graph_objects as go
                df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]
                # df = df[['compound_name', '7P2G']]
                colors = ['lightslategray',] * len(df)
                
                novo = df.loc[df['7P2G']<=-6.1]
                indices = novo.index
                valores = indices.values
                for i in valores:
                    colors[i] = 'crimson'
                    
                fig = make_subplots(rows=1, cols=1, shared_yaxes=True, subplot_titles=("7P2G (SARS-CoV-2)"))

                fig = go.Figure(data=[go.Bar(
                    x=df['compound_name'],
                    y=df['7P2G'],
                    marker_color=colors # marker color can be a single color value or an iterable
                )])

                fig.add_hline(y=-6.1, line_width=3, line_dash = 'dash', line_color = 'black', annotation_text = "-6.2 baseline 7P2G")
                fig.update_layout(title_text='Protein 7P2G - Results Docking')
                st.plotly_chart(fig)

                st.markdown('#### Docking results - 4DD8 - Asthma')
                from plotly.subplots import make_subplots
                import plotly.graph_objects as go
                # df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]
                import plotly.graph_objects as go
                df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]
                colors = ['lightslategray',] * len(df)
                novo = df.loc[df['4DD8']<=-6.4]
                indices = novo.index
                valores = indices.values
                for i in valores:
                    colors[i] = 'sienna'
                    
                fig = make_subplots(rows=1, cols=1, shared_yaxes=True, subplot_titles=("4DD8 (Asthma)"))

                fig = go.Figure(data=[go.Bar(
                    x=df['compound_name'],
                    y=df['4DD8'],
                    marker_color=colors # marker color can be a single color value or an iterable
                )])

                fig.add_hline(y=-6.4, line_width=3, line_dash = 'dash', line_color = 'black', annotation_text = "-6.4 baseline 4DD8")
                fig.update_layout(title_text='Protein 4DD8 - Results Docking')
                st.plotly_chart(fig)

                st.markdown('#### Docking results - 1NC6 - Asthma')
                from plotly.subplots import make_subplots
                import plotly.graph_objects as go
                # df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]
                import plotly.graph_objects as go
                df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]
                colors = ['lightslategray',] * len(df)
                novo = df.loc[df['1NC6']<=-6.2]
                indices = novo.index
                valores = indices.values
                for i in valores:
                    colors[i] = 'mediumaquamarine'
                    
                fig = make_subplots(rows=1, cols=1, shared_yaxes=True, subplot_titles=("1NC6 (Asthma)"))

                fig = go.Figure(data=[go.Bar(
                    x=df['compound_name'],
                    y=df['1NC6'],
                    marker_color=colors # marker color can be a single color value or an iterable
                )])

                fig.add_hline(y=-6.2, line_width=3, line_dash = 'dash', line_color = 'black', annotation_text = "-6.2 baseline 1NC6")
                fig.update_layout(title_text='Protein 1NC6 - Results Docking')
                st.plotly_chart(fig)

                st.markdown('#### Docking results - 6VVU - Asthma')
                from plotly.subplots import make_subplots
                import plotly.graph_objects as go
                # df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]
                import plotly.graph_objects as go
                df = df[['compound_name', '7P2G', '4DD8', '1NC6', '6VVU']]
                colors = ['lightslategray',] * len(df)
                novo = df.loc[df['6VVU']<=-4.3]
                indices = novo.index
                valores = indices.values
                for i in valores:
                    colors[i] = 'steelblue'
                    
                fig = make_subplots(rows=1, cols=1, shared_yaxes=True, subplot_titles=("6VVU (Asthma)"))

                fig = go.Figure(data=[go.Bar(
                    x=df['compound_name'],
                    y=df['6VVU'],
                    marker_color=colors # marker color can be a single color value or an iterable
                )])

                fig.add_hline(y=-4.3, line_width=3, line_dash = 'dash', line_color = 'black', annotation_text = "-4.3 baseline 6VVU")
                fig.update_layout(title_text='Protein 6VVU - Results Docking')
                st.plotly_chart(fig)

            else:
                st.write('Insert your data')


    elif choice == ("Molecular Docking"):
        from PIL import Image 
        image = Image.open('Lumios_logo.png')
        st.sidebar.image(image, use_column_width=True)
        image = Image.open('Docking.png')
        st.sidebar.image(image, use_column_width=True)
        st.write("""
        ## LUMIOS - Molecular Docking""")
        # st.write("""
        #         # Desreplicação Molecular""")
                
        # st.write("###LUMIOS - Docking")

        # title = st.text_input('Diretório onde salvar as informacoes da docagem')

        data_file = st.file_uploader("Upload CSV",type=['csv'], accept_multiple_files=False)
        if st.button("Yes, Docking!"):
            if data_file is not None:
                import pandas as pd
                df = pd.read_csv(data_file)
                st.dataframe(df)
                # df.drop('Unnamed: 0.1', axis=1, inplace=True)
                # df.drop('Unnamed: 0', axis=1, inplace=True)
                df.rename(columns = {'index': 'index_base1'}, inplace=True)
                df.rename(columns = {'SMILES': 'smiles'}, inplace=True)
                st.dataframe(df)

                df.reset_index(inplace=True)
                st.dataframe(df)

                # df.rename(columns = {'index': 'index_base'}, inplace=True)
                # df.rename(columns = {'Unnamed: 0', 'index'}, inplace=True)
                # df.reset_index(inplace=True)

                # st.dataframe(df)

                ids = []
                for i in df['id']:
                    a = str(i)
                    ids.append(a)

                df['id'] = ids

                # st.dataframe(df)

                mols = []
                from rdkit import Chem
                from rdkit.Chem import AllChem

                for _, row in df.iterrows():
                    m = Chem.MolFromSmiles(row.smiles)
                    
                    m = Chem.AddHs(m)
                    
                    AllChem.EmbedMolecule(m, AllChem.ETKDG())
                    minimize_status = AllChem.UFFOptimizeMolecule(m, 2000)
                    
                    if not minimize_status == 0:
                        print(f"Failed to minimize_compound'{row['name']}")
                        
                    AllChem.ComputeGasteigerCharges(m)
                    
                    mols.append(m)            


                import pandas as pd
                from rdkit import Chem
                from rdkit.Chem import AllChem

                import numpy as np
                import os
                # salvar_np = os.path.join(title, "lig_dataset_novo.sdf")
                
                sdf_writer = Chem.SDWriter('lig_dataset_novo.sdf')

                lig_properties = df.columns.to_list()
                for i, mol in enumerate(mols):
                    data_ref = df.iloc[i]
                    mol.SetProp('index', '%s'%i)
                    mol.SetProp('_Name', data_ref['id'])
                    for p in lig_properties:
                        mol.SetProp(f"+{p}", str(data_ref[p]))    
                    sdf_writer.write(mol)
                sdf_writer.close()

                from openbabel import pybel
                import os
                # os.makedirs(out_dir_lig, exist_ok = True)
                import tempfile
                # tempfile.mkdtemp('teste_teste')
                # print(tempfile.mkdtemp('teste_teste'))
                dir = r'C:\Users\PC\Desktop\LUMIOS\Projeto_final\teste_teste'
                import os
                if not os.path.exists(dir):
                    os.makedirs(dir)
                # os.mkdir('./teste_teste')
                # if not os.path.exists('my_folder'):
                #     os.mkdir('./teste_teste')

                # import os
                # import shutil

                # path = './teste_teste'
                # if not os.path.exists(path):
                #     os.makedirs(path)
                # else:
                #     shutil.rmtree(path)           # Removes all the subdirectories!
                #     os.makedirs(path)

                for mol in pybel.readfile('sdf', 'lig_dataset_novo.sdf'):
                    mol.write('pdbqt', 'teste_teste/%s.pdbqt' %mol.data['index'], overwrite=True)

                # import os
                # import pandas as pd

                RECEPTOR = ['7P2G', '4DD8', '1NC6', '6VVU']
                # RECEPTOR = ['6VVU']


                import os

                df_docagem = [] 
                for i in  RECEPTOR:
                    WORK_DIR = r'C:\Users\PC\Desktop\LUMIOS\Projeto_final'
                    # LIG_DIR = os.path.join(WORK_DIR, 'teste_teste') #que será temporária
                    LIG_DIR = os.path.join(WORK_DIR, 'teste_teste')
                    RECEPTOR_DIR = os.path.join(WORK_DIR, 'receptors', i)
                    print(RECEPTOR_DIR)
                    # RECEPTOR_DIR = os.path.join(WORK_DIR, '7P2G') #pasta da proteina- FIXA

                    os.mkdir('./runs2/')
                    # OUT_DIR = os.path.join(WORK_DIR, "runs2", RECEPTOR)
                    OUT_DIR = os.path.join(WORK_DIR, "runs2")
                    
                    print(OUT_DIR)


                    ligands = []

                    for file in os.listdir(LIG_DIR):
                        if (os.path.isfile(os.path.join(LIG_DIR, file)) & file.endswith('.pdbqt')):
                            ligands.append(file)
                            
                    print(len(ligands))

                    from shutil import copy2

                    prepared_lig_dirs = []

                    try:

                        for lig in ligands:
                            lig_filename = os.path.splitext(lig)[0]
                            out_dir_lig = os.path.join(OUT_DIR, lig_filename)

                            os.makedirs(out_dir_lig, exist_ok = True)

                            copy2(os.path.join(RECEPTOR_DIR, f"{i}final.pdbqt"), out_dir_lig)
                            copy2(os.path.join(RECEPTOR_DIR, f"{i}config.txt"), out_dir_lig)

                            copy2(os.path.join(LIG_DIR, lig), out_dir_lig)

                            prepared_lig_dirs.append(out_dir_lig)

                    except OSError as error:
                        print(error)
                        print("Directory '% s' can not be removed" % OUT_DIR)

                    import shlex, subprocess
                    from datetime import datetime

                    # def run_docking (lig_out_dir, output_logs, i):

                    #     print(f"Running {lig_out_dir}")

                    #     output_logs = output_logs + f"\n[STARTRUN] {datetime.now()} OUTDIR {lig_out_dir}\n[STARTLOG]\n"

                    #     ligand = f"{os.path.basename(lig_out_dir)}.pdbqt"
                    #     vina = r'C:\Users\PC\Desktop\Vina\vina.exe'

                    #     #args = shlex.split(f"{vina} --config_7P2G.txt --ligand {ligand}")
                    #     pdbqt = i+'_final.pdbqt'
                    #     config = i+'config.txt'
                    #     args = [r'C:\Users\PC\Desktop\Curso_Docagem\Vina\vina.exe', '--receptor', pdbqt, '--config', 
                    #             config, '--ligand', ligand, '--log', 'results.txt']
                    #     print(args)
                    #     #args[0] = 'C:\Program Files (x86)\The Scripps Research Institute\Vina\vina.exe'
                    #     process = subprocess.Popen(
                    #     args,
                    #     stdout = subprocess.PIPE,
                    #     stderr = subprocess.PIPE,
                    #     cwd = os.path.join(lig_out_dir)
                    #     )
                    #     output, error = process.communicate()

                    #     if error:
                    #         print("Error: ", error)
                    #         output_logs = output_logs + output.decode("utf-8")

                    #     else:
                    #         print(output.decode("utf-8"))
                    #         output_logs = output_logs + output.decode("utf-8")
                    #     output_logs = output_logs + f"\n[ENDLOG]\n[ENDRUN] {datetime.now()}\n+++++++++++++\n"

                    #     print('\n+++++++++++++++\n')

                    #     return output_logs

                    output_logs = ""
                    print(i)
                    print(len(prepared_lig_dirs))
                    for j in prepared_lig_dirs:
                        print(j)
                        # output_logs = run_docking(j, output_logs, i)
                        output_logs = output_logs + f"\n[STARTRUN] {datetime.now()} OUTDIR {j}\n[STARTLOG]\n"

                        ligand = f"{os.path.basename(j)}.pdbqt"
                        vina = r'C:\Users\PC\Desktop\Vina\vina.exe'

                        #args = shlex.split(f"{vina} --config_7P2G.txt --ligand {ligand}")
                        # pdbqt = i+'final.pdbqt'
                        # config = i+'config.txt'
                        print(f"{i}final.pdbqt")
                        args = [r'C:\Users\PC\Desktop\Curso_Docagem\Vina\vina.exe', '--receptor', f"{i}final.pdbqt", '--config', 
                                f"{i}config.txt", '--ligand', ligand, '--log', 'results.txt']
                        print(args)
                        #args[0] = 'C:\Program Files (x86)\The Scripps Research Institute\Vina\vina.exe'
                        process = subprocess.Popen(
                        args,
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE,
                        cwd = os.path.join(j)
                        )
                        output, error = process.communicate()

                        if error:
                            print("Error: ", error)
                            output_logs = output_logs + output.decode("utf-8")

                        else:
                            print(output.decode("utf-8"))
                            output_logs = output_logs + output.decode("utf-8")
                        output_logs = output_logs + f"\n[ENDLOG]\n[ENDRUN] {datetime.now()}\n+++++++++++++\n"

                        print('\n+++++++++++++++\n')
                        
                        WORK_DIR = r'C:\Users\PC\Desktop\LUMIOS\Projeto_final\runs2'
                        ligands = []
                        for file in os.listdir(WORK_DIR):
                            ligands.append(file)
                            
                        print(ligands)

                        prepared_lig_dirs = []
                        for b in ligands:
                            #lig_filename = os.path.splitext(lig)[0]
                            out_dir_lig = os.path.join(WORK_DIR, b)
                            prepared_lig_dirs.append(out_dir_lig)
                        print(prepared_lig_dirs)

                        endereco = []
                        for q in prepared_lig_dirs:
                            a = q + '\\results.txt'
                            print(a)
                            endereco.append(a)
                            
                        print(endereco)

                    df_docagem1 = []
                    for t in endereco:
                        print(t)
                        data = pd.read_csv(t, header = None, on_bad_lines='skip')
                        print('data')
                        print(len(data))

                        if len(data) == 12:
                            df_docagem1.append(0)

                        else:

                            best_pose = pd.DataFrame(data[0][20].split('      ')).T
                            try:
                                a=best_pose.iat[0, 1].strip(" ")
                            except:
                                a=0


                            df_docagem1.append(a)
                    df_docagem.append(df_docagem1)

                    print(df_docagem1)
                    print(df_docagem)


                    import shutil
                    shutil.rmtree(r'C:\Users\PC\Desktop\LUMIOS\Projeto_final\runs2')

                    docking = pd.DataFrame(df_docagem).T

                    for i, j in enumerate(RECEPTOR):
                                # print(i, j)
                                #a = j
                        docking.rename(columns = {i: j}, inplace=True)

                    docking.reset_index(inplace=True)

                    docking['index'] = ligands
                    docking['index'] = docking['index'].astype(int)

                    anotacoes_com_docagem = df.merge(docking, how='inner', on='index')
                    st.write('Annotations with the best poses')
                    st.dataframe(anotacoes_com_docagem)


                @st.cache
                def convert_df(df):
                    # IMPORTANT: Cache the conversion to prevent computation on every rerun
                    return df.to_csv(index=False).encode('utf-8')

                csv = convert_df(anotacoes_com_docagem)

                st.download_button(
                label="Download data as CSV",
                data=csv,
                file_name='samples_machine_learning_docking.csv',
                mime='text/csv',
                )


                import shutil
                shutil.rmtree(r'C:\Users\PC\Desktop\LUMIOS\Projeto_final\\teste_teste')
                # shutil.rmtree(r'C:\Users\PC\Desktop\Projeto_final\runs2')
                # st.write('deu certo', i)

            else:
                st.write('Insert the annotations')

    else:
        st.subheader("About")
        # st.dataframe(amostras)
        from PIL import Image 
        image = Image.open('rafael.png')

        st.image(image, use_column_width=True)
        # st.info("Built with Streamlit")
        st.info("Rafael Vieira")
        



if __name__ == '__main__':
    main()
