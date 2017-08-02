## cBioPortal-new-study-assistant

#### Required python packages
 * NumPy
 * SciPy
 * pandas
 * seaborn
 * python-levenshtein: [anaconda link](https://anaconda.org/conda-forge/python-levenshtein)

#### Running the assistant
The default mode of the assistant selects a random study from the portal and searches other studies on the portal for matching attributes.
```javascript
python new_study_assistant.py
```
##### Options available
```javascript
--study_to_drop
```

#### Output
Running the script results in several files:
 * dendrogram.png
 * n_attribute_distribution.png
 * n_common_attribute_distribution.png
 * n_unique_attribute_distribution.png
