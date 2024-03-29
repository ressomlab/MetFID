import pandas as pd
from tabulate import tabulate


def retrieve_compound(fingerprint, database, mass_list, ppm):
    """
    Given a search type, predicted fingerprint, a database, and mass(optional),
    returns a table that includes all the possible compounds.
    :param fingerprint: predicted fingerprint (len = 528)
    :param database: database.csv
    :param mass_list: list of precursor masses
    :param ppm: mass tolerance in ppm (default 10)
    :return: dict of possible compounds

    REQ: For now, the database is not provided by Ressom's lab, and it must be
    a .csv file.
    """
    compound_dict = {}

    candidate_dict = retrieve_candidate(database, mass_list, ppm)

    try:
        for candidate, fp in candidate_dict.items():
            tanimoto = calculate_tanimoto(fp, fingerprint)
            compound_dict[candidate] = tanimoto
    except AttributeError:
        print("Compound Not Found!")

    return compound_dict


def retrieve_candidate(database, mass_list, ppm):
    """
    Given a search type, a database, mass(optional), and inchikey(optional)
    returns a dict that contains all qualified candidates.
    :param database: database.csv
    :param mass_list: list of precursor mass
    :param ppm: mass tolerance in ppm
    :return: dict{(Name, Inchikey): candidate fingerprint}

    REQ: The search type must be one of 'mass' or 'formula'. For now, the
    database is not provided by Ressom's lab, and it must be a .csv file.
    """
    df = pd.read_csv(database)
    candidate_dict = {}
    for mass in mass_list:
        min_weight = mass * 1000000.0 / (1000000.0 + ppm)
        max_weight = mass * 1000000.0 / (1000000.0 - ppm)
        candidate_df = df.loc[df['Mass'].between(min_weight, max_weight, inclusive=True)]

        if len(candidate_df) == 0:
            continue

        candidate_list = candidate_df[['Name', 'Inchikey', 'Formula', '5618_fp']].values.tolist()

        for candidate in candidate_list:
            fp = candidate[3]
            try:
                fp = list(fp[1:len(fp) - 1].split(','))
            except:
                continue
            fp = [int(i) for i in fp]
            candidate_dict[candidate[0], candidate[1], candidate[2]] = fp

    return candidate_dict


def calculate_tanimoto(real_fp, predicted_fp):
    """
    Given a compound fingerprint and a predicted fingerprint, calculate the
    tanimoto score.
    :param real_fp: candidate's fingerprint
    :param predicted_fp: predicted fingerprint
    :return: tanimoto score
    """
    predicted_fingerprint = list(map(lambda x: int(x), predicted_fp))
    fp_length = len(predicted_fingerprint)

    if fp_length != len(real_fp):
        raise ValueError('Different size of vector!')

    identity_num = 0
    length = 0
    for i in range(fp_length):
        if predicted_fingerprint[i] == real_fp[i] == 1:
            identity_num += 1
            length += 1
        elif predicted_fingerprint[i] != real_fp[i]:
            length += 1

    return identity_num / length


def visualize_compound_dict(compound_dict, compound_name=True):
    """
    Given a compound dict, prints out the dict as a table in an order of
    descending tanimoto scores.
    :param compound_dict: compound dict
    :param compound_name: if we can obtain the compound name
    :return: None
    """
    sort_list = []
    sorted_dict = {k: v for k, v in sorted(compound_dict.items(), key=lambda item: item[1], reverse=True)}

    if compound_name:
        for k, v in sorted_dict.items():
            sort_list.append([str(k[0]), str(k[1]), str(k[2]), str(v)])

        return tabulate(sort_list, headers=['Compound Name', 'Inchikey', 'Formula', 'Score'])
    else:
        for k, v in sorted_dict.items():
            sort_list.append([str(k), str(v)])

        return tabulate(sort_list, headers=['Inchikey', 'Score'])


def get_inchikey_score(inchikey_dict, predicted_fp):
    """
    Given a dict that has inchikeys as 'key' and binned fingerprint as 'value',
    and the predicted fingerprint vector, returns another dict that has
    inchikeys as 'key' and tanimoto score as 'value'.
    :param inchikey_dict: dict{inchikey: fingerprint}
    :param predicted_fp: predicted fingerprint vector
    :return: dict{inchikey: score}
    """
    inchikey_score_dict = {}
    for key, fp_vec in inchikey_dict.items():
        inchikey_score_dict[key] = calculate_tanimoto(fp_vec, predicted_fp)

    return inchikey_score_dict
