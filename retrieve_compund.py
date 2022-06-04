import pandas as pd
import pubchempy
from tabulate import tabulate


def retrieve_compound(search_type, fingerprint, database, ppm=20, mass=None, inchikey=None):
    """
    Given a search type, predicted fingerprint, a database, and mass(optional),
    returns a table that includes all the possible compounds.
    :param search_type: search type
    :param fingerprint: predicted fingerprint (len = 528)
    :param database: database.csv
    :param ppm: mass tolerance in ppm (default 20)
    :param mass: precursor mass (optional, required when search_type is mass)
    :param inchikey: inchikey (optional, required when search_type is formula)
    :return: dict of possible compounds

    REQ: For now, the database is not provided by Ressom's lab, and it must be
    a .csv file.
    """
    candidate_dict = None
    compound_dict = {}

    if search_type == 'formula':
        candidate_dict = retrieve_candidate(search_type, database, inchikey=inchikey)
    elif search_type == 'mass':
        candidate_dict = retrieve_candidate(search_type, database, ppm, mass=mass)

    try:
        for candidate, fp in candidate_dict.items():
            tanimoto = calculate_tanimoto(fp, fingerprint)
            compound_dict[candidate] = tanimoto
    except AttributeError:
        print("Compound Not Found!")

    return compound_dict


def retrieve_candidate(search_type, database, ppm=None, mass=None, inchikey=None):
    """
    Given a search type, a database, mass(optional), and inchikey(optional)
    returns a dict that contains all qualified candidates.
    :param search_type: search type
    :param database: database.csv
    :param ppm: mass tolerance in ppm (optional, required when search_type is mass)
    :param mass: precursor mass (optional, required when search_type is mass)
    :param inchikey: inchikey (optional, required when search_type is formula)
    :return: dict{(Name, Inchikey): candidate fingerprint}

    REQ: The search type must be one of 'mass' or 'formula'. For now, the
    database is not provided by Ressom's lab, and it must be a .csv file.
    """
    df = pd.read_csv(database)
    candidate_df = None

    if search_type == 'formula':
        try:
            # omics_craft database
            formula = df.loc[df['Inchikey'] == inchikey]["Formula"].iloc[0]
            # print("get formula from database")
        except:
            # pubchempy
            formula = pubchempy.get_compounds(identifier=inchikey, namespace="inchikey")[0].molecular_formula
        candidate_df = df.loc[df["Formula"] == formula]
    elif search_type == 'mass':
        min_weight = mass * 1000000.0 / (1000000.0 + ppm)
        max_weight = mass * 1000000.0 / (1000000.0 - ppm)
        candidate_df = df.loc[df['Mass'].between(min_weight, max_weight, inclusive=True)]

    if len(candidate_df) == 0:
        return None

    candidate_list = candidate_df[['Name', 'Inchikey', 'fp_vec']].values.tolist()
    candidate_dict = {}
    for candidate in candidate_list:
        fp = candidate[2]
        try:
            fp = list(fp[1:len(fp) - 1].split(','))
        except:
            continue
        fp = [int(i) for i in fp]
        candidate_dict[candidate[0], candidate[1]] = fp

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
            sort_list.append([str(k[0]), str(k[1]), str(v)])

        return tabulate(sort_list, headers=['Compound Name', 'Inchikey', 'Score'])
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
