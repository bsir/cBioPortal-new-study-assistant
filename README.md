## cBioPortal-new-study-assistant

### Required python packages
 * NumPy
 * SciPy
 * pandas
 * seaborn
 * python-levenshtein: [anaconda link](https://anaconda.org/conda-forge/python-levenshtein)
 * Pylatex (required only for the '--output_pdf' option)

### Running the script
The default mode of the script selects a random study from the portal and searches other studies on the portal for matching attributes.  The current version of the script typically takes a few minutes to run and depends on internet access to download data from cBioPortal.

Default example:
```bash
python new_study_assistant.py
```

Example using acyc_mda_2015 raw data (this data is provided in the acyc_mda_2015 folder in this repository):
```bash
python new_study_assistant.py --study_to_drop='acyc_mda_2015' --new_study_path='./acyc_mda_2015/raw_data_clinical.txt' > similarity_output.txt
```

#### Options available
```bash
--new_study_path PATH
```
 * Path to raw study data file that you want to analyze. Currently the code assumes that the file only contains attribute names in the header.  If the file contains a multi-line header the program will probably crash.

```bash
--study_to_drop STUDY_ID
```
 * Excludes a study (specified by STUDY_ID) from the analysis, this is useful when analyzing a study already on cBioPortal

```bash
--specific_study STUDY_ID
```
 * Use this option if you want to run the analysis on a specific study (specified by STUDY_ID) that is already in the cBioPortal

```bash
--output_pdf
```
 * Activate this flag if you would like the report results printed to a pdf (requires Pylatex)

```bash
--datahub_path PATH
```
 * Specify path to a local version of the datahub, if this path is not specified the script will download data via the API instead


### Output
Similar attributes that are detected in the test study are printed to the screen.  The prefix "NEW_STUDY_" is added to each attribute in the test study to distinguish those attributes from those already on cBioPortal.  The text output from the script can be redirected to a file by using "> FILENAME" at the end of the python command.

Running the script also results in several image files:
 * dendrogram.png - An image showing the dendrogram obtained for similarity detection based on attribute values.
 * n_attribute_distribution.png - An image showing the number of attributes contained in the test study relative to all studies on the cBioPortal.
 * n_common_attribute_distribution.png - An image showing the number of attributes in common with studies on the cBioPortal.
 * n_unique_attribute_distribution.png - An image showing the number of unique attributes in the test study relative to all studies on the cBioPortal.
