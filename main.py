import read_and_predict as rap
import warnings
warnings.filterwarnings(action="ignore")

if __name__ == '__main__':
    # for lab members (use the database to search)
    rap.msms_predict('_files/testing_compound.txt')

    # for other users
    # rap.inchikey_predict('_files/testing_compound.txt', '_files/inchikey_list.txt')
