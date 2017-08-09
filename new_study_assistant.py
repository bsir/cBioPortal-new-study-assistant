### load in packages
import argparse
import seaborn
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Levenshtein
import scipy
import urllib, json
from collections import Counter
import math
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

#define some functions
#-----------------------------------
#downloads main study information from the portal
def get_all_study_info():
    url = "http://www.cbioportal.org/api/studies"
    response = urllib.urlopen(url)
    json_data = json.loads(response.read())
    return json_data

#downloads all attribute data for a list of studies
def get_attribute_data(study_list):
    studies=[]
    study_attributes = []

    #loop over studies and download data from api
    for study in study_list:
        try:
            studyID = study['studyId']
        except:
            studyID = study
        study_api_url = 'http://cbioportal.org/api/studies/' + studyID + '/clinical-attributes'
        df = pd.read_json(study_api_url)

        #make sure the study is not empty
        if not df.empty:
            studies.append((studyID ,df))
            study_attributes.append((studyID ,df['clinicalAttributeId'].tolist()))
            
    return studies, study_attributes

def get_clinical_data_values(study_list, attribute_type='PATIENT'):
    studies_clinical_data=[]
    study_attributes = []

    #loop over studies and download data from api
    for study in study_list:
        try:
            studyID = study['studyId']
        except:
            studyID = study
        study_api_url = 'http://cbioportal.org/api/studies/' + studyID + '/clinical-data?clinicalDataType=' + attribute_type
        df = pd.read_json(study_api_url)
        #make sure the study is not empty
        if not df.empty:
            df=df.pivot(values='value', index='entityId', columns='clinicalAttributeId').fillna(value=np.nan)
            studies_clinical_data.append((studyID ,df))
            #study_attributes.append((studyID ,df['clinicalAttributeId'].tolist()))
    return studies_clinical_data

def get_study_names(study_attribute_list):
    names=[]
    for study in study_attribute_list:
        names.append(study[0])
    return names

def study_attributes_to_df(study_attribute_data):
    #transform data into boolean table of studies/attributes
    #entries in the table are 1 if a study contains an attribute and 0 otherwise
    study_attribute = []
    for i in study_attribute_data:
        for j in i[1]:
            study_attribute.append((i[0],j))

    study_attribute_pairs = pd.DataFrame.from_records(study_attribute, columns = ['study','attribute'])

    study_data_nolabel = pd.get_dummies(study_attribute_pairs['attribute'])

    study_data_combined = pd.concat([study_attribute_pairs['study'], study_data_nolabel], axis=1)
    study_data_combined = study_data_combined.groupby('study').sum()
        
    return study_data_combined

def drop_study(study_to_drop, combined_data):
    combined_data.drop(study_to_drop, axis=0, inplace=True)
    combined_data.drop([col for col, val in combined_data.sum().iteritems() if val == 0], axis=1, inplace=True)
    return combined_data

def find_attribute_name_matches(cBioPortal_attributes, new_attributes, cutoff=0.9):
    all_col_names = list(cBioPortal_attributes)
    lev_dist = np.zeros([len(all_col_names), len(new_attributes)])
    for i in range(len(all_col_names)):
        for j in range(len(new_attributes)):
            lev_dist[i,j]=Levenshtein.ratio(unicode(all_col_names[i].upper()), unicode(new_attributes[j].upper()))
            
    all_lev_distances = pd.DataFrame(data=lev_dist.T, index=new_attributes, columns=all_col_names)
    matches = pd.DataFrame(data=list(all_lev_distances[all_lev_distances > cutoff].stack().index))
            
    return matches

def output_results(exact_matches, possible_matches, non_matching_attributes):
    sample_or_patient = {True:"patient", False:"sample"}

    print "================================================" 
    print "attributes with exact matches:"
    print "================================================"

    for attribute in exact_matches:
        print attribute + " is present in " + str(int(study_data_combined.sum()[attribute])) + " other studies"
    
    ############################################################
    
    print ""
    print "================================================"
    print "attributes with possible matches:"
    print "================================================"

    no_matches = []
    possible_match_count = 0
    for attribute in non_matching_attributes:
        #print matching attributes
        try:
            filtered_matches = np.setdiff1d(possible_matches[possible_matches[0]==attribute][1].values, exact_matches)
        except:
            filtered_matches = np.empty(0)
        if (filtered_matches.size)>0:
            possible_match_count += 1
            #print attribute in new study along with datatype and patient/sample attribute
            print "================================================"
            print "new study attribute: " + attribute
            print '---------------------------------------'
            
            #print similar attributes found in existing studies along with the datatype
            print 'possible matches:'

            for matching_attribute in filtered_matches:
                print matching_attribute + ", which is present in " + str(int(study_data_combined.sum()[matching_attribute])) + " other studies"
                attribute_dtypes=[]
                attribute_ps_types=[]

        else:
            no_matches.append(attribute)

    if possible_match_count == 0:
        print "No similar matches detected."
    ############################################################3
    print ""
    print "================================================" 
    print "attributes with NO matches:"
    print "================================================"

    for attribute in no_matches:
        print attribute
        
def plot_attribute_distribution(cBioPortal_data, new_study_data):
    plt.rcParams['figure.figsize'] = (10.0, 8.0)
    plt.rcParams.update({'font.size':20})
    plt.figure()
    ax=seaborn.distplot(cBioPortal_data.sum(axis=1), kde=False)
    ax.set(ylabel='number of studies', xlabel='attributes in study')
    plt.axvline(len(new_study_data), color='k', linestyle='dashed', linewidth=2)
    plt.savefig('n_attribute_distribution.png')
    plt.close()

def plot_unique_and_common_attribute_distributions(cBioPortal_data, non_matching_attributes):
    #unique attributes plot
    plt.rcParams['figure.figsize'] = (10.0, 8.0)
    plt.rcParams.update({'font.size':20})
    plt.figure()
    unique_attributes_cBio_studies = cBioPortal_data.T[(cBioPortal_data.sum(axis=0)==1).values].sum(axis=0)
    ax=seaborn.distplot(unique_attributes_cBio_studies, kde=False)
    ax.set(ylabel='number of studies', xlabel='unique attributes in study')
    plt.axvline(non_matching_attributes.size, color='k', linestyle='dashed', linewidth=2)
    plt.savefig('n_unique_attribute_distribution.png')
    plt.close()
    
    #common attributes plot
    plt.rcParams['figure.figsize'] = (10.0, 8.0)
    plt.rcParams.update({'font.size':20})
    plt.figure()
    common_attributes = cBioPortal_data.sum(axis=1) - unique_attributes_cBio_studies
    ax=seaborn.distplot(common_attributes, kde=False)
    ax.set(ylabel='number of studies', xlabel='common attributes in study')
    plt.axvline(exact_matching_attributes.size, color='k', linestyle='dashed', linewidth=2)
    plt.savefig('n_common_attribute_distribution.png')
    plt.close()
    
def get_new_study_attributes(test_study_data):
    if "PATIENT_ID" in test_study_data:
        del test_study_data["PATIENT_ID"]
    if "SAMPLE_ID" in test_study_data:   
        del test_study_data["SAMPLE_ID"]
    test_study_attribute_names = map(unicode,map(str.upper,list(test_study_data)))
    return test_study_attribute_names

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    
    try:
        import unicodedata
        unicodedata.numeric(s)
    except (TypeError, ValueError):
        pass
    
    return False

def process_clinical_data(all_study_clinical_data, study_to_drop=''):
    attribute_data = []
    for study in all_study_clinical_data:
        if study[0] != study_to_drop:
            study_data_attributes = list(study[1])

            #get attribute data to filter based on datatype
            if api_flag:
                study_api_url = 'http://cbioportal.org/api/studies/' + study[0] + '/clinical-attributes'
                df = pd.read_json(study_api_url)
                for attribute in study_data_attributes:
                    if not attribute.endswith('_ID'):
                        if df['datatype'][df['clinicalAttributeId']==attribute].values[0] == 'STRING':
                            data = set(study[1][attribute])
                            for d in data:
                                if not is_number(d):
                                    attribute_data.append((attribute.upper(), d.upper()))
            
            else:
                for attribute in study_data_attributes:
                    if (not unicode(attribute).endswith("ID")) and ("FILE_NAME" not in attribute)  and ("DATE" not in attribute):
                        data = set(study[1][attribute])
                        for d in data:
                            if not is_number(d):
                                attribute_data.append((attribute.upper(), d.upper()))

    attribute_data_pairs = pd.DataFrame.from_records(attribute_data, columns = ['attribute','data'])
    attribute_data_nolabel = pd.get_dummies(attribute_data_pairs['data'])
    attribute_data_combined = pd.concat([attribute_data_pairs['attribute'], attribute_data_nolabel], axis=1)
    attribute_data_combined = attribute_data_combined.groupby('attribute').sum()
    attribute_data_combined[attribute_data_combined>0]=1

    #drop values which only occur in single attribute
    attribute_data_combined.drop([col for col, val in attribute_data_combined.sum().iteritems() if val > 10], axis=1, inplace=True)
    attribute_data_combined.drop([row for row, val in attribute_data_combined.sum(axis=1).iteritems() if val < 1], axis=0, inplace=True)
    
    return attribute_data_combined

def process_new_study_data(new_study_data):
    ns_attribute_data = []
    for attribute in list(new_study_data):
        if (attribute != u'SAMPLE_ID') and (attribute != u'PATIENT_ID'):
            data = set(new_study_data[attribute])
            for d in data:
                if not is_number(d):
                    ns_attribute_data.append((attribute.upper(), d.upper()))
    
    ns_attribute_data_pairs = pd.DataFrame.from_records(ns_attribute_data, columns = ['attribute','data'])

    ns_attribute_data_nolabel = pd.get_dummies(ns_attribute_data_pairs['data'])

    ns_attribute_data_combined = pd.concat([ns_attribute_data_pairs['attribute'], ns_attribute_data_nolabel], axis=1)
    ns_attribute_data_combined = ns_attribute_data_combined.groupby('attribute').sum()

    ns_attribute_data_combined[ns_attribute_data_combined>0]=1
    ns_attribute_data_combined = ns_attribute_data_combined.T.add_prefix('NEW_STUDY_').T
    
    return ns_attribute_data_combined

def get_clusters(attribute_value_data):
    data_link = scipy.cluster.hierarchy.linkage(attribute_value_data, method='complete', metric='cosine') # computing the linkage
    thresh = 0.7*max(data_link[:,2])
    clusters = fcluster(data_link, thresh, 'distance')
    clustered_attributes=[]

    for label in np.unique(clusters):
        clustered_attributes.append(list(np.asarray(attribute_value_data.index.values)[clusters==label]))
        
    plt.rcParams['figure.figsize'] = (180.0, 25.0)
    plt.rcParams.update({'font.size':40})
    plt.figure()
    den=dendrogram(data_link, labels = combined_attribute_values.index, leaf_font_size = 15, above_threshold_color='#AAAAAA')
    plt.savefig('dendrogram.png', bbox_inches='tight')
    plt.close()
        
    return clustered_attributes

def output_cluster_matches(all_clusters):
    print ""
    print "================================================"
    print "possible attribute matches based on data values"
    print "================================================"


    for cluster in all_clusters:
        if len(cluster)>1:
            for attribute in cluster:
                if 'NEW_STUDY_' in attribute:
                    matches = np.setdiff1d(cluster, attribute)
                    match_count = 0
                    for match in matches:
                        if 'NEW_STUDY_' not in match:
                            match_count += 1
                            if match_count == 1:
                                print "================================================"
                                print "possible matches for: " + attribute
                                print "--------------------"
                            print "possible match: " + match

def find_data_files(basedir):
    import os
    import os.path

    file_list = []
    for dirpath, dirnames, filenames in os.walk(basedir):
        for filename in [f for f in filenames if ((("clinical" in f) and ("data" in f)) and (f.endswith(".txt")))]:
            file_list.append(os.path.join(dirpath, filename))
    return file_list

def load_data_from_files(file_list, basedir):
    study_data = []
    attribute_data = []
    for clinical_file in file_list:
        #load in individual data from files and append to a list
        df = pd.read_table(clinical_file, skiprows=0)
        if (list(df)[0]!='SAMPLE_ID' and list(df)[0]!='PATIENT_ID'):
            rows_to_skip=df[(df[df.columns[0]]=="SAMPLE_ID") | (df[df.columns[0]]=="PATIENT_ID") | (df[df.columns[0]]=="OTHER_PATIENT_ID") | (df[df.columns[0]]=="OTHER_SAMPLE_ID")].index[0]+1
            df = pd.read_table(clinical_file, skiprows=rows_to_skip)
        df.columns = map(str.upper, df.columns)
        col_names = list(df)
        tcga_provisional_flag = False
        study_name = clinical_file.split(basedir+'/')[1].split('/')[0]
        if clinical_file.split(basedir+'/')[1].split('/')[1] == 'tcga':
            study_name = study_name + '_tcga'
        if study_name == '':
            study_name = clinical_file.split(basedir+'/')[1].split('/')[1]
            if clinical_file.split(basedir+'/')[1].split('/')[2] == 'tcga':
                study_name = study_name + '_tcga'
        study_data.append((study_name, col_names))
        attribute_data.append((study_name, df))
    return study_data, attribute_data

def process_levenshtein_matches(exact_matches, possible_matches, non_matching_attributes):
    match_pairs = []
    no_matches = []
    for attribute in non_matching_attributes:
        try:
            filtered_matches = np.setdiff1d(possible_matches[possible_matches[0]==attribute][1].values, exact_matches)
        except:
            filtered_matches = np.empty(0)
        if (filtered_matches.size)>0:
            for matching_attribute in filtered_matches:
                match_pairs.append((attribute, matching_attribute))
        else:
            no_matches.append(attribute)
    return match_pairs, no_matches

def process_cluster_matches(all_clusters):
    match_pairs = []
    for cluster in all_clusters:
        if len(cluster)>1:
            for attribute in cluster:
                if 'NEW_STUDY_' in attribute:
                    matches = np.setdiff1d(cluster, attribute)
                    for match in matches:
                        if 'NEW_STUDY_' not in match:
                            match_pairs.append((attribute[10:], match))
    return match_pairs

def get_unique_attributes(unique_attributes, value_matches):
    no_matches = []
    for attribute in unique_attributes:
        if sum(value_matches['New study attribute']==attribute)==0:
            no_matches.append(attribute)
    return no_matches

def generate_latex_report(match_table,):
    from pylatex import NoEscape, Document, Section, Subsection, Tabular, Math, TikZ, Axis, Plot, Figure, Matrix, Alignat
    from pylatex.utils import italic
    import os

    file_directory = ''

    geometry_options = {"tmargin": "1in", "lmargin": "1in", "bmargin": "1in", "rmargin": "1in"}
    doc = Document(geometry_options=geometry_options)

    with doc.create(Section('cBioPortal new study report')):
        doc.append('This is an example report')

        with doc.create(Subsection('Attribute matches')):
            with doc.create(Tabular('|c|c|')) as table:
                table.add_hline()
                table.add_row(list(match_table))
                table.add_hline()
                table.add_hline()

                for row in match_table.index:
                    table.add_row(list(match_table.loc[row,:]))
                    table.add_hline()

        doc.append(NoEscape(r'\newpage'))

        with doc.create(Subsection('Number of attribute distributions')):
            with doc.create(Figure(position='h!')) as n_attributes:
                n_attributes.add_image(file_directory + 'n_attribute_distribution.png', width=NoEscape(r'0.99\textwidth'))
                n_attributes.add_caption('Number of attributes')

        #with doc.create(Subsection('Number of unique attributes')):
            with doc.create(Figure(position='h!')) as n_unique_attributes:
                n_unique_attributes.add_image(file_directory + 'n_unique_attribute_distribution.png', width=NoEscape(r'0.99\textwidth'))
                n_unique_attributes.add_caption('Number of unique attributes')

        #with doc.create(Subsection('Number of common attributes')):
            with doc.create(Figure(position='h!')) as n_common_attributes:
                n_common_attributes.add_image(file_directory + 'n_common_attribute_distribution.png', width=NoEscape(r'0.99\textwidth'))
                n_common_attributes.add_caption('Number of common attributes')

    doc.generate_pdf('report', clean_tex=True)


#values that should be read in
# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--new_study_path", help="path to new study data")
parser.add_argument("--study_to_drop", help="if study being tested is already on cBioPortal it may be helpful to drop that study from the analysis")
parser.add_argument("--specific_study", help="name of specific study on the portal to test")
parser.add_argument("--datahub_path", help="directory where the datahub is located")
parser.add_argument("--output_pdf", action='store_true', help="flag to output pdf")
args = parser.parse_args()

new_study_path = args.new_study_path
study_to_drop = args.study_to_drop
specific_study = args.specific_study
datahub_path = args.datahub_path
output_pdf = args.output_pdf
random_study=False
if (new_study_path is None) and (specific_study is None):
    random_study = True

api_flag=True
if datahub_path is not None:
    api_flag = False
    print "Using cBioPortal data from local datahub: " + datahub_path
else:
    print "cBioPortal data will be obtained via API"

latex_flag=False
if output_pdf:
    latex_flag=True
    print "PDF file will be output"


#main part of the script below
#main function
if api_flag:
    all_cBioPortalStudies = get_all_study_info()
    studies, study_attributes = get_attribute_data(all_cBioPortalStudies)
else:
    list_of_files = find_data_files(datahub_path)
    study_attributes, attribute_value_data = load_data_from_files(list_of_files, datahub_path)
study_names = get_study_names(study_attributes)
study_data_combined = study_attributes_to_df(study_attributes)


#choose a study to serve as a "new" study
if random_study or specific_study is not None:
    if random_study:
        test_study = random.choice(study_names)
    else:
        test_study = specific_study
    print "Study being analyzed: " + test_study
    study_to_drop = test_study
    study_data_combined = drop_study(test_study, study_data_combined)
    test_study_data, test_study_clin_attributes = get_attribute_data([test_study])
    test_study_attribute_data = test_study_data[0][1]
    test_study_attribute_names = test_study_attribute_data['clinicalAttributeId'].values
    test_study_patient_attribute_values = get_clinical_data_values([test_study], attribute_type="PATIENT")
    test_study_sample_attribute_values = get_clinical_data_values([test_study], attribute_type="SAMPLE")
    processed_test_study_patient_attribute_values = process_new_study_data(test_study_patient_attribute_values[0][1])
    processed_test_study_sample_attribute_values = process_new_study_data(test_study_sample_attribute_values[0][1])
    new_study_clinical_attribute_values = pd.concat([processed_test_study_patient_attribute_values,
                                                     processed_test_study_sample_attribute_values],
                                                     axis=0).fillna(value=0)

#otherwise read in data from a new study
else:
    new_study_data = pd.read_table(new_study_path)
    test_study_attribute_names = get_new_study_attributes(new_study_data)
    if study_to_drop is not None:
        study_data_combined = drop_study(study_to_drop, study_data_combined)
    new_study_clinical_attribute_values = process_new_study_data(new_study_data)
    
#check matching attributes (via name only)
exact_matching_attributes = np.intersect1d(test_study_attribute_names, list(study_data_combined))
non_matching_attributes = np.setdiff1d(test_study_attribute_names, list(study_data_combined))
possible_matches = find_attribute_name_matches(study_data_combined, non_matching_attributes)

#check matching attributes based on data values
if api_flag:
    patient_data = get_clinical_data_values(all_cBioPortalStudies, attribute_type="PATIENT")
    patient_data_processed = process_clinical_data(patient_data, study_to_drop)
    sample_data = get_clinical_data_values(all_cBioPortalStudies, attribute_type="SAMPLE")
    sample_data_processed = process_clinical_data(sample_data, study_to_drop)
    combined_PS_attribute_data = pd.concat([patient_data_processed, sample_data_processed], axis=0).fillna(value=0)
else:
    combined_PS_attribute_data = process_clinical_data(attribute_value_data, study_to_drop)
combined_attribute_values = pd.concat([combined_PS_attribute_data, new_study_clinical_attribute_values], axis=0).fillna(value=0)

#get clusters and plot dendrogram
clusters = get_clusters(combined_attribute_values)

#process matches
processed_name_matches, no_name_matches = process_levenshtein_matches(exact_matching_attributes, 
                                                                      possible_matches,
                                                                      non_matching_attributes)
processed_value_matches = process_cluster_matches(clusters)

#combine matches into one dataframe
exact_matches_df = pd.DataFrame(np.column_stack([exact_matching_attributes, exact_matching_attributes]), 
                                columns=['New study attribute', 'Possible matches'])
name_matches_df = pd.DataFrame(processed_name_matches, columns=['New study attribute', 'Possible matches'])
value_matches_df = pd.DataFrame(processed_value_matches, columns=['New study attribute', 'Possible matches'])
no_name_or_value_matches = get_unique_attributes(no_name_matches, value_matches_df)
no_matches_df = pd.DataFrame(np.asarray(([no_name_or_value_matches, ['No matches found']*len(no_name_or_value_matches)])).T,
                             columns=['New study attribute', 'Possible matches'])
all_matches_df = (pd.concat([exact_matches_df, value_matches_df, no_matches_df])
                  .drop_duplicates()
                  .sort_values('New study attribute',axis=0)
                  .reset_index(drop=True))

#print all_matches_df
if latex_flag:
    generate_latex_report(all_matches_df)
else:
    print all_matches_df


#make and save figures
plot_attribute_distribution(study_data_combined, test_study_attribute_names)
plot_unique_and_common_attribute_distributions(study_data_combined, non_matching_attributes)
