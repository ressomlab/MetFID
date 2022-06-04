import pubchempy
from openbabel import pybel
import numpy as np


def data_process(msms_data_list):
    """
    Given a MS/MS data list which contains the first element in the list
    represents the precursor m/z, retention time(in minutes), and ion mode.
    The remaining are m/z and intensity pairs. Returns a dict that contains
    precursor mass, retention time, ion mode, m/z and intensity pairs.
    :param msms_data_list: MS/MS data list
    :return: dict{precursor, rt, mode, [m/z], [intensity]}
    """
    msms_data_dict = {}
    mz = []
    intensity = []

    precursor_mode = msms_data_list[0].split(' ')
    msms_data_dict['precursor'] = float(precursor_mode[0])
    mode = precursor_mode[1].rstrip()

    if mode == 'positive':
        msms_data_dict['mode'] = 'positive'
        msms_data_dict['precursor'] -= 1.007276
    elif mode == 'negative':
        msms_data_dict['mode'] = 'negative'
        msms_data_dict['precursor'] += 1.007276
    else:
        raise ValueError("The ion mode should be 'positive' or 'negative'.")

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
    :return: dict{precursor, rt, mode, [m/z], [intensity]}
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
    :return: dict{precursor, rt, mode, [m/z], [intensity]}
    """
    tem_list = []
    for i in msms_data_dict['intensity']:
        if i >= 2:
            tem_list.append(i)

    if len(tem_list) < 5:
        raise UserWarning('Too few peaks.')

    return msms_data_dict


def denoise(msms_data_dict):
    """
    Given a MSMS data dict, remove the intensity pair that has m/z larger
    than the precursor mass.
    :param msms_data_dict: dict{precursor, rt, mode, [m/z], [intensity]}
    :return: dict{precursor, rt, mode, [m/z], [intensity]}
    """
    mass = msms_data_dict['precursor']
    for i in msms_data_dict['m/z']:
        if i > mass:
            index_to_remove = msms_data_dict['m/z'].index(i)
            msms_data_dict['m/z'].pop(index_to_remove)
            msms_data_dict['intensity'].pop(index_to_remove)

    return msms_data_dict


def binning(msms_data_dict):
    """
    Given a MSMS data dict, binning the m/z range of each MS/MS spectrum into
    pre-specified bins, which indicate continuous integer m/z values, and
    calculate the accumulated intensities within each bin as feature values.
    :param msms_data_dict: dict{precursor, rt, [m/z], [intensity]}
    :return: binned vector of length 1174
    """
    first_digit = 5
    intensity_digit = 0
    input_vec = [0] * 1174  # 14 - 1178, 1174

    mz_list = [int(round(i, 0)) - first_digit for i in msms_data_dict['m/z']]

    for i in mz_list:
        input_vec[i] += msms_data_dict['intensity'][intensity_digit]
        intensity_digit += 1

    # input_vec = one2two(input_vec)  # convert binned vector from 1d to 2d

    return input_vec


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
        inchikey_dict[key] = fp_vec

    return inchikey_dict


def fp_conversion(fp3, fp4, macc):
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
