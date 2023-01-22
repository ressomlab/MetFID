import pickle

import numpy as np
import pubchempy
from openbabel import pybel
from PyFingerprint.fingerprint import get_fingerprint


def calculate_adduct(pm, mode):
    """
    Given a precursor mass, ion mode, and an adduct(optional), returns all
    possible mass after considering the adducts. We only consider the
    followings:
    1. Positive ion mode
    M+H, M+NH4, M+Na, M+H-H2O, M+K, M+ACN+H, M+ACN+Na, M+2Na-H, M+2H, M+3H,
    M+H+Na, M+2H+Na, M+2Na, M+2Na+H, M+Li, M+CH3OH+H

    2. Negative ion mode
    M-H, M-H2O-H, M+Na-2H, M+Cl, M+K-2H, M+FA-H, M-2H, M-3H, M+CH3COO, M+F

    Adduct data was obtained from: https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/
    :param pm: precursor mass
    :param mode: ionization mode
    :return: List[all possible masses]
    """
    # adduct_table = {
    #     'M+H': pm - 1.007276, 'M+NH4': pm - 18.033823,
    #     'M+Na': pm - 22.989218, 'M+H-H2O': pm + 17.00687,
    #     'M+K': pm - 38.963158, 'M+ACN+H': pm - 42.033823,
    #     'M+ACN+Na': pm - 64.015765, 'M+2Na-H': pm - 44.971160,
    #     'M+2H': 2 * (pm - 1.007276), 'M+3H': 3 * (pm - 1.007276),
    #     'M-H': pm + 1.007276, 'M-H2O-H': pm + 19.01839,
    #     'M+Na-2H': pm - 20.974666, 'M+Cl': pm - 34.969402,
    #     'M+K-2H': pm - 36.948606, 'M+FA-H': pm - 44.998201,
    #     'M-2H': 2 * (pm + 1.007276), 'M-3H': 3 * (pm + 1.007276),
    #     'M+CH3COO': pm - 59.04078, 'M+F': pm - 18.99840
    # }
    pos = {'M+H': pm - 1.007276, 'M+NH4': pm - 18.033823, 'M+Na': pm - 22.989218}
    neg = {'M-H': pm + 1.007276, 'M+Cl': pm - 34.969402, 'M+FA-H': pm - 44.998201}

    if mode == 'positive':
        return [i for i in pos.values()]
    elif mode == 'negative':
        return [i for i in neg.values()]
    else:
        raise ValueError("The ion mode should be 'positive' or 'negative'.")


def data_process(msms_data_list):
    """
    Given a MS/MS data list which contains the first element in the list
    represents the precursor m/z. The remaining are m/z and intensity pairs.
    Returns a dict that contains precursor mass, ion mode, m/z and intensity
    pairs.
    :param msms_data_list: MS/MS data list
    :return: dict{[precursor masses], mode, [m/z], [intensity]}
    """
    msms_data_dict = {}
    mz = []
    intensity = []

    precursor_mode = msms_data_list[0].split(' ')
    msms_data_dict['precursor'] = calculate_adduct(float(precursor_mode[0]), precursor_mode[1])
    msms_data_dict['mode'] = precursor_mode[1]

    for i in msms_data_list[1:]:
        mz_intensity = i.split(' ')
        mz.append(float(mz_intensity[0]))
        intensity.append(float(mz_intensity[1]))

    msms_data_dict['m/z'] = mz
    msms_data_dict['intensity'] = intensity

    return msms_data_dict


def scaling(msms_data_dict):
    """
    Given a MSMS data dict, scaling the MSMS data if there exists some
    intensities greater than 100%.
    :param msms_data_dict: dict{precursor, rt, mode, [m/z], [intensity]}
    :return: dict{[precursor masses], mode, [m/z], [intensity]}
    """
    max_intensity = max(msms_data_dict['intensity'])
    rate = 100 / max_intensity
    if max_intensity > 100:
        for i in range(len(msms_data_dict['intensity'])):
            msms_data_dict['intensity'][i] *= rate

    return msms_data_dict


def filtering(msms_data_dict):
    """
    Given a MSMS data dict, raise a error when spectra that consisted of fewer
    than five peaks with relative intensity above 2%.

    Note: for package only, not for the testing.
    :param msms_data_dict: dict{precursor, rt, mode, [m/z], [intensity]}
    :return: dict{[precursor masses], mode, [m/z], [intensity]}
    """
    tem_list = []
    for i in msms_data_dict['intensity']:
        if i >= 2:
            tem_list.append(i)

    if len(tem_list) < 5:
        raise UserWarning('Too few peaks.')

    return msms_data_dict


# def denoise(msms_data_dict):
#     """
#     Given a MSMS data dict, remove the intensity pair that has m/z larger
#     than the precursor mass.
#     :param msms_data_dict: dict{precursor, rt, mode, [m/z], [intensity]}
#     :return: dict{[precursor masses], mode, [m/z], [intensity]}
#     """
#     mass = msms_data_dict['precursor']
#     for i in msms_data_dict['m/z']:
#         if i > mass:
#             index_to_remove = msms_data_dict['m/z'].index(i)
#             msms_data_dict['m/z'].pop(index_to_remove)
#             msms_data_dict['intensity'].pop(index_to_remove)
#
#     return msms_data_dict


def binning(msms_data_dict):
    """
    Given a MSMS data dict, binning the m/z range of each MS/MS spectrum into
    pre-specified bins, which indicate continuous integer m/z values, and
    calculate the accumulated intensities within each bin as feature values.
    :param msms_data_dict: dict{[precursor masses], mode, [m/z], [intensity]}
    :return: binned vector of length 40,088
    """
    first_digit = 5
    input_vec = [0] * 117331
    mz_list = [int(round((i - first_digit) * 100, 0)) for i in msms_data_dict['m/z']]

    for i in enumerate(mz_list):
        input_vec[i[1]] = msms_data_dict['intensity'][i[0]]

    with open('_files/spectra_to_add_indexes.p', 'rb') as index_file:
        index_list = pickle.load(index_file)

    spectra_vec = [input_vec[i] for i in index_list]
    # input_vec = one2two(input_vec)  # convert binned vector from 1d to 2d

    return spectra_vec


# TODO: for 2d CNN
def one2two(binned_vector):
    """
    Given a binned vector, convert it to matrix with shape (35, 34)
    :param binned_vector: binned vector with length 1174
    :return: a matrix with shape (35, 34)
    """
    data = np.array(binned_vector + [0] * 16)  # fill empty entry with 0's
    shape = (35, 34)
    return data.reshape(shape)


def get_binned(inchikey_list):
    """
    Given a text that contains inchikeys, returns a dict that has inchikeys as
    'key' and binned fingerprint as 'value'.
    :param inchikey_list: a text file that contains inchikeys
    :return: dict{inchikey: fingerprint}
    """
    inchikey_dict = {}

    for key in inchikey_list:
        key = key.rstrip()
        smiles = pubchempy.get_compounds(identifier=key, namespace='inchikey')[0].canonical_smiles
        mol = pybel.readstring('smi', smiles)
        fp_vec = fp_conversion(mol.calcfp('FP3').bits, mol.calcfp('FP4').bits, mol.calcfp('MACCS').bits)
        full_length_fp = fp_vec + obtain_fingerprint(smiles)
        inchikey_dict[key] = obtain_5618_fingerprint(full_length_fp)

    return inchikey_dict


def fp_conversion(fp3, fp4, macc):
    """
    Converts fp3, fp4, and MACCs to a new fingerprint.
    :param fp3: FP3 fingerprint (55 digits)
    :param fp4: FP4 fingerprint (307 digits)
    :param macc: MACCS (166 digits)
    :return: combined fingerprint with length of 528
    """
    fp3_temp = list(range(1, 56))
    fp4_temp = list(range(1, 308))
    macc_temp = list(range(1, 167))
    fp3_box = [0] * len(fp3_temp)
    fp4_box = [0] * len(fp4_temp)
    macc_box = [0] * len(macc_temp)

    for i in range(len(fp3_temp)):
        if fp3_temp[i] in fp3:
            fp3_box[i] = 1

    for i in range(len(fp4_temp)):
        if fp4_temp[i] in fp4:
            fp4_box[i] = 1

    for i in range(len(macc_temp)):
        if macc_temp[i] in macc:
            macc_box[i] = 1

    return fp3_box + fp4_box + macc_box


def obtain_fingerprint(smi):
    """
    Given a canonical SMILES, returns the extended part of the fingerprint. The
    extended part includes ECFP + PubChem + Klekota-roth, with length
    1024 + 881 + 4860 = 6765
    :param smi: canonical SMILES
    :return: extended fingerprint with length 6765
    """
    cdktypes = ['extended', 'pubchem', 'klekota-roth']

    output = {}
    for f in cdktypes:
        output[f] = get_fingerprint(smi, f)

    output_np = output.copy()
    for k, fp in output.items():
        output_np[k] = fp.to_numpy()

    fp = []
    for k, v in output_np.items():
        fp += v.tolist()

    return [int(i) for i in fp]


def obtain_5618_fingerprint(full_length_fp):
    """
    Given the full length fingerprint (length of 7293), trimmed it to length of
    5618.
    :param full_length_fp: full length fingerprint
    :return: 5618 fingerprint
    """
    new_fp = []
    with open('_files/fp_to_add_indexes.p', 'rb') as f:
        fp_indexes_list = pickle.load(f)

    for fp_i in fp_indexes_list:
        new_fp.append(full_length_fp[fp_i])

    if len(new_fp) == 5618:
        return new_fp
    else:
        print('fp length error')
        return
