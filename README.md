## cBioPortal-new-study-assistant

#### Required python packages
 * NumPy
 * SciPy
 * pandas
 * seaborn
 * python-levenshtein: [anaconda link](https://anaconda.org/conda-forge/python-levenshtein)

#### Running the assistant
The default mode of the assistant selects a random study from the portal and searches other studies on the portal for matching attributes.
```bash
python new_study_assistant.py
```

##### Options available
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

#### Output
Running the script results in several files:
 * dendrogram.png 
 * n_attribute_distribution.png
 * n_common_attribute_distribution.png
 * n_unique_attribute_distribution.png
