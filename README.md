# MetFID
CNN-Based Compound Fingerprint Prediction for Metabolite Annotation

## Getting Started
Before running the program, you need to have the [PyFingerprint][7], [PubChemPy 1.0.4][1], [Open Babel 3.1.1][2], [tabulate 0.8.9][3] and [Tensorflow][4] installed.

And also, please download the latest CNN model -- [MetFID_CNN_40088_5618.h5][5] and put it in the root directory of MetFID.

## Prerequisites
This program works with Python 3.5+.

## Description

In the folder [`_files`][6], there are two required files. The first one is `testing_compound.txt`, which is a list of spectrum data; The second file is `inchikey_list.txt`, which contains lists of InChIKeys. 

Here is an example for `testing_compound.txt` and `inchikey_list.txt`.

* `testing_compound.txt`

The first row represents the precursor mass and ionization mode, followed by intensity pairs.
```
#1
389.1626 positive
58.0653 0.1393717711334875
70.0652 0.5051454118468733
72.0808 7.5707678873637185
84.0682 0.25425091544424155
129.0701 0.21914480352209237
130.0732 0.1374693989440672
145.0647 0.13117470906587816
165.0698 9.97754058994659
166.0777 29.57116152422922
183.0805 1.83382519459832
187.1077 1.3620057326501762
193.0761 1.1618155468334308
199.031 0.2417147860532003
201.0465 100
389.1626 17.431109736689535
#2
255.0299 negative
171.0449 1.405280299496928
183.0453 1.869843966679493
227.035 20.56615845658887
255.03 100
```
* `inchikey_list.txt`

This file contains the InChIKey list. The block followed by `#digit` are the lists of InChIKeys that user suspect as the true compounds.
```
#1
CECDPVOEINSAQG-UHFFFAOYSA-N
SGUAFYQXFOLMHL-UHFFFAOYSA-N
XSJLXZULZZXNJI-UHFFFAOYSA-N
ZKLPARSLTMPFCP-UHFFFAOYSA-N
HEWDOWUUTBCVJP-UHFFFAOYSA-N
#2
VILFVXYKHXVYAB-UHFFFAOYSA-N
XXOYUNQCBYNWNL-UHFFFAOYSA-N
JLZANJDGGAPGOF-UHFFFAOYSA-N
BBNQQADTFFCFGB-UHFFFAOYSA-N
RFHAOTPXVQNOHP-UHFFFAOYSA-N
```

Notice that, each list is separated by `#digit`. There will be an error pop up if the counter does not match in both files.


## Running the Program
* For Ressom Lab users  
Before running the program, please **uncomment** line 7 and **comment** line 10 in `main.py`, then run:

```
$ python3 main.py
```

* For public users  
After downloading the code, navigating to MetFID folder, then run:

```
$ python3 main.py
```

There will be an output file `testing_compound_prediction.txt` be created in the folder `_files`. Using the example above, we will get:

```
#1
Inchikey                        Score
---------------------------  --------
ZKLPARSLTMPFCP-UHFFFAOYSA-N  0.742547
SGUAFYQXFOLMHL-UHFFFAOYSA-N  0.381818
CECDPVOEINSAQG-UHFFFAOYSA-N  0.290323
XSJLXZULZZXNJI-UHFFFAOYSA-N  0.268844
HEWDOWUUTBCVJP-UHFFFAOYSA-N  0.201946
#2
Inchikey                        Score
---------------------------  --------
BBNQQADTFFCFGB-UHFFFAOYSA-N  0.586826
JLZANJDGGAPGOF-UHFFFAOYSA-N  0.373913
RFHAOTPXVQNOHP-UHFFFAOYSA-N  0.296629
VILFVXYKHXVYAB-UHFFFAOYSA-N  0.290419
XXOYUNQCBYNWNL-UHFFFAOYSA-N  0.167984
```

The first column represents the `InChIKeys`, and the second column represents the `Tanimoto similarity score`. Each table will be ranked in a descending order by score.


[1]:https://pubchempy.readthedocs.io/en/latest/guide/install.html
[2]:https://openbabel.org/wiki/Python
[3]:https://pypi.org/project/tabulate/
[4]:https://www.tensorflow.org/
[5]:https://drive.google.com/file/d/1KCvvjjRQk4Q9PUC7ajF20sTWxd6nBcnP/view?usp=sharing
[6]:https://github.com/aohongyu/MetFID/tree/main/_files
[7]:https://github.com/hcji/PyFingerprint
