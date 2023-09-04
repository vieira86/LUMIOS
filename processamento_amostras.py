def processamento_amostras(path_data):    
    import os
    import pandas as pd
#     path_data = title

    class_folders2 = sorted(os.listdir(path_data))
    class_folders2.remove('Lotus_inchi.csv')
    class_folders2.remove('results_vanila_02.json')

    from matchms.filtering import default_filters
    from matchms.filtering import repair_inchi_inchikey_smiles
    from matchms.filtering import derive_inchikey_from_inchi
    from matchms.filtering import derive_smiles_from_inchi
    from matchms.filtering import derive_inchi_from_smiles
    from matchms.filtering import harmonize_undefined_inchi
    from matchms.filtering import harmonize_undefined_inchikey
    from matchms.filtering import harmonize_undefined_smiles
    def metadata_processing(spectrum):
        spectrum = default_filters(spectrum)
        spectrum = repair_inchi_inchikey_smiles(spectrum)
        spectrum = derive_inchi_from_smiles(spectrum)
        spectrum = derive_smiles_from_inchi(spectrum)
        spectrum = derive_inchikey_from_inchi(spectrum)
        spectrum = harmonize_undefined_smiles(spectrum)
        spectrum = harmonize_undefined_inchi(spectrum)
        spectrum = harmonize_undefined_inchikey(spectrum)
        return spectrum


    from matchms.filtering import default_filters
    from matchms.filtering import normalize_intensities
    from matchms.filtering import select_by_intensity
    from matchms.filtering import select_by_mz
    def peak_processing(spectrum):
        spectrum = default_filters(spectrum)
        spectrum = normalize_intensities(spectrum)
        spectrum = select_by_intensity(spectrum, intensity_from=0.2)
        spectrum = select_by_mz(spectrum, mz_from=10, mz_to=1000)
        return spectrum



    from matchms.importing import load_from_msp
    import pandas as pd
    import numpy as np
    id_experiment = []
    peaks =[]
    precursortype = []
    modo = []
    precursor_mz = []
    intensity = []
    retention_time = []

    for i in class_folders2:
    #     id_experiment = []
    #     peaks =[]
    #     precursortype = []
    #     modo = []
    #     precursor_mz = []
        file_mgf = os.path.join(path_data, i)
        id_experiment.append(i)
        spectrums_MassBank = list(load_from_msp(file_mgf))
        spectrums_MassBank = [metadata_processing(s) for s in spectrums_MassBank]
        spectrums_MassBank = [peak_processing(s) for s in spectrums_MassBank]

        peaks.append(spectrums_MassBank[0].peaks.mz)
        if spectrums_MassBank[0].metadata.get('adduct') is None:
            precursortype.append(spectrums_MassBank[0].metadata.get('precursortype'))

        else:
            precursortype.append(spectrums_MassBank[0].metadata.get('adduct'))

        modo.append(spectrums_MassBank[0].metadata.get('ionmode'))
        precursor_mz.append(spectrums_MassBank[0].metadata.get('precursor_mz'))

        intensity.append(spectrums_MassBank[0].metadata.get('intensity'))
        retention_time.append(spectrums_MassBank[0].metadata.get('retention_time'))

    df = pd.DataFrame()
    df['id'] = id_experiment
    df['peaks'] = peaks
    df['precursor_type'] = precursortype
    df['modo'] = modo
    df['precursor_mz'] = precursor_mz
    df['intesity'] = intensity
    df['retention_time'] = retention_time

    cacau = df

    parent_mass =[]
    for i in range(len(df)):
        if (cacau['precursor_type'][i] == '[M+H]+'):
            a = cacau['precursor_mz'][i]-1.007825
            parent_mass.append(a)
        elif (cacau['precursor_type'][i] == '[M+Na]+'):
            b = cacau['precursor_mz'][i]-22.989770
            parent_mass.append(b)

        elif (cacau['precursor_type'][i] == '[M+K]+'):
            c = cacau['precursor_mz'][i]-38.963708
            parent_mass.append(c)

    cacau['parent_mass'] = parent_mass
    picos = []
    for i in cacau['peaks']:
        a = np.around(i,  decimals=2)
        picos.append(a)

    cacau['peaks'] = picos
    return cacau