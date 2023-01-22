import re

import msms_data_process as mdp
import prediction
import retrieve_compund as rc


def msms_predict(msms_file_path):
    """
    Given a MS/MS data file which contains the first element in the file
    represents the precursor m/z, retention time(in minutes), and ion mode.
    The remaining are m/z and intensity pairs. Outputs the predictions in a
    text file.
    :param msms_file_path: file path for msms data
    :return: None
    """
    with open(msms_file_path) as msms_file:
        msms_data = msms_file.read()

    processed_msms = re.split('#', msms_data)[1:]  # split file by keyword '#'
    msms_counter = ['#' + i.split('\n')[0] for i in processed_msms]

    for i in processed_msms:
        msms_list = i.split('\n')[1:-1]
        msms_dict = mdp.data_process(msms_list)
        msms_dict_scale = mdp.scaling(msms_dict)
        # msms_dict_denoised = mdp.denoise(msms_dict_scale)
        binned_vec = mdp.binning(msms_dict_scale)
        predicted_fp = prediction.predict_fingerprint(binned_vec)
        compound_dict = rc.retrieve_compound(predicted_fp, '_files/massDB_SMILE_5618fp_CASMI16+22.csv', msms_dict['precursor'], 10)

        output_file = msms_file_path.split('.')[0] + '_prediction.txt'
        with open(output_file, 'a') as result:
            result.write(msms_counter.pop(0) + '\n')
            result.write(rc.visualize_compound_dict(compound_dict) + '\n')


def inchikey_predict(msms_file_path, inchikey_file_path):
    """
    Given a MS/MS data file which contains the first element in the file
    represents the precursor m/z, and ion mode.
    The remaining are m/z and intensity pairs. The other file is a .txt that
    includes some inchikeys. Ranks the inchikeys in order.
    :param msms_file_path: file path for msms data
    :param inchikey_file_path: file path for inchikey data
    :return:
    """
    with open(msms_file_path) as msms_file:
        msms_data = msms_file.read()

    with open(inchikey_file_path) as inchikey_file:
        inchikey_data = inchikey_file.read()

    processed_msms = re.split('#', msms_data)[1:]  # split file by keyword '#'
    processed_inchikey = re.split('#', inchikey_data)[1:]

    msms_counter = ['#' + i.split('\n')[0] for i in processed_msms]
    inchikey_counter = ['#' + i.split('\n')[0] for i in processed_inchikey]

    if msms_counter != inchikey_counter:
        raise ValueError('Counters for two files should be constant.')

    for i, j in zip(processed_msms, processed_inchikey):
        msms_list = i.split('\n')[1:-1]
        inchikey_list = j.split('\n')[1:-1]

        msms_dict = mdp.data_process(msms_list)
        inchikey_dict = mdp.get_binned(inchikey_list)

        msms_dict_scale = mdp.scaling(msms_dict)
        binned_vec = mdp.binning(msms_dict_scale)
        predicted_fp = prediction.predict_fingerprint(binned_vec)
        inchikey_score = rc.get_inchikey_score(inchikey_dict, predicted_fp)

        output_file = msms_file_path.split('.')[0] + '_prediction.txt'
        with open(output_file, 'a') as result:
            result.write(inchikey_counter.pop(0) + '\n')
            result.write(rc.visualize_compound_dict(inchikey_score, False) + '\n')
