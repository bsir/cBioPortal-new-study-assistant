## cBioPortal-new-study-assistant

### Required python packages
 * NumPy
 * SciPy
 * pandas
 * seaborn
 * python-levenshtein: [anaconda link](https://anaconda.org/conda-forge/python-levenshtein)

### Running the assistant
The default mode of the assistant selects a random study from the portal and searches other studies on the portal for matching attributes.
Default example:
```bash
python new_study_assistant.py
```

Example using acyc_mda_2015 raw data:
```bash
python new_study_assistant.py --study_to_drop='acyc_mda_2015' --new_study_path='./acyc_mda_2015/raw_data_clinical.txt' > similarity_output.txt
```

#### Options available
```bash
--new_study_path
```
 * Path to raw study data file that you want to analyze. Currently the code assumes that the file only contains attribute names in the header.  If the file contains a multi-line header the program will probably crash.


```bash
--study_to_drop
```
 * Excludes a study from the analysis, this is useful when analyzing a study already on cBioPortal

```bash
--specific_study
```
 * Use this option if you want to run the analysis on a study that is already in the cBioPortal

### Output
Running the script results in several files:
 * dendrogram.png - An image showing the dendrogram obtained for similarity detection based on attribute values.
 * n_attribute_distribution.png - An image showing the number of attributes contained in the test study relative to all studies on the cBioPortal.
 * n_common_attribute_distribution.png - An image showing the number of attributes in common with studies on the cBioPortal.
 * n_unique_attribute_distribution.png - An image showing the number of unique attributes in the test study relative to all studies on the cBioPortal.
