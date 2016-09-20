# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 08:55:41 2016

Reads cohort from Cambridge and extracts NSPN variables

@author: maria j. rosa
"""

# Read subject's information from main NSPN database
import csv
import numpy as np

# Get all baseline scans
f = open('/Users/maria/Documents/NSPN/code/nspn_subjects_baseline.csv', 'rU')
csv_f = csv.reader(f)
all_baseline_scans = []
all_baseline_centre = []
all_baseline_ids = []
for row in csv_f:
    all_baseline_scans.append(row[2])
    all_baseline_centre.append(row[1])
    all_baseline_ids.append(row[0])

# Get gender
f = open('/Users/maria/Documents/NSPN/code/nspn_subjects_baseline_gender_info.csv', 'rU')
csv_f = csv.reader(f)
all_gender = []
for row in csv_f:
    all_gender.append(row[1])

# Get cohort
f = open('/Users/maria/Documents/NSPN/code/nspn_subjects_baseline_cohort_info.csv', 'rU')
csv_f = csv.reader(f)
all_cohort = []
for row in csv_f:
    all_cohort.append(row[1])

# Get cohort
f = open('/Users/maria/Documents/NSPN/code/nspn_subjects_baseline_age_info.csv', 'rU')
csv_f = csv.reader(f)
all_age = []
for row in csv_f:
    all_age.append(row[1])

# Get scans
f = open('/Users/maria/Documents/NSPN/code/NSPN_scans_all.csv', 'rU')
csv_f = csv.reader(f)
ids_mri = []
ids_centre = []
ids_gender = []
ids_cohort = []
ids_age = []
depressed_id = []
k = 0
for row in csv_f:
    indx = all_baseline_scans.index(row[0])
    nspn_id = all_baseline_ids[indx]
    centre_id = all_baseline_centre[indx]
    gender_id = all_gender[indx]
    cohort_id = all_cohort[indx]
    k = k + 1
    if cohort_id == 'Depression':
        depressed_id.append(k)
        print k
        print cohort_id
    age_id = all_age[indx]
    ids_mri.append(nspn_id)
    ids_centre.append(centre_id)
    ids_gender.append(gender_id)
    ids_cohort.append(cohort_id)
    ids_age.append(age_id)


# Read HPQ1
f = open('/Users/maria/Documents/NSPN/docs/NSPN_HPQ1_processed.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_hpq = []
data_hpq_ids = []
for row in csv_f:
    data_hpq.append(row)
    data_hpq_ids.append(row[0])

# Read IUA
f = open('/Users/maria/Documents/NSPN/docs/NSPN_IUA1_processed_110816v5.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_iua = []
data_iua_ids = []
for row in csv_f:
    data_iua.append(row)
    data_iua_ids.append(row[0])

# Get weight and height
f = open('/Users/maria/Documents/NSPN/docs/NSPN_weight_height_processed.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_wh = []
data_wh_ids = []
for row in csv_f:
    data_wh.append(row)
    data_wh_ids.append(row[0])

# Read baseline subjects
f = open('/Users/maria/Documents/NSPN/docs/UCHANGE_CompleteCohort_BaselineOnly_20150722.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_baseline = []
baseline_vars = []
data_wh_baseline = []
baseline_wh = []
age_cambridge = []
for row in csv_f:
    indx_hpq = data_hpq_ids.index(row[0])
    indx_iua = data_iua_ids.index(row[0])
    indx_wh = data_wh_ids.index(row[0])
    tmp = data_hpq[indx_hpq] + data_iua[indx_iua][1::]
    tmp_wh = data_wh[indx_wh]
    baseline_vars.append(tmp[0])
    data_baseline.append(tmp)
    baseline_wh.append(tmp_wh[0])
    data_wh_baseline.append(tmp_wh)
    age_cambridge.append(row[4])

# Save variables data - Cambridhe cohort
ids_good = []
mri_good = []
gender_good = []
cohort_good = []
age_good = []
with open("/Users/maria/Documents/NSPN/docs/NSPN_vars_baseline_test.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_mri)[0]):
        if (baseline_vars.count(ids_mri[i])) > 0:
            indx = baseline_vars.index(ids_mri[i])
            writer.writerow(data_baseline[indx])
            ids_good.append(i+1)
            mri_good.append(ids_centre[i])
            gender_good.append(ids_gender[i])
            cohort_good.append(ids_cohort[i])
            age_good.append(ids_age[i])

# Save age Cambridge
with open("/Users/maria/Documents/NSPN/docs/NSPN_age_baseline_cambridge.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_mri)[0]):
        if (baseline_vars.count(ids_mri[i])) > 0:
            indx = baseline_vars.index(ids_mri[i])
            writer.writerow([i,age_cambridge[indx]])


# Save weight and height data - Cambridhe cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_weight_height_baseline.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_mri)[0]):
        if (baseline_wh.count(ids_mri[i])) > 0:
            indx = baseline_wh.index(ids_mri[i])
            writer.writerow(data_wh_baseline[indx])

# Save IDs - Cambridhe cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_MRIids_baseline.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_good)[0]):
        writer.writerow([ids_good[i],ids_good[i]])

# Save centre information - Cambridhe cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_gender_bin_baseline.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_good)[0]):
        if gender_good[i] == 'Female':
            writer.writerow([i,0])
        else:
            writer.writerow([i,1])

# Save gender information - Cambridhe cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_cohort_baseline.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_good)[0]):
        writer.writerow([i,cohort_good[i]])

# Save gender information - Cambridhe cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_age_baseline.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_good)[0]):
        writer.writerow([i,age_good[i]])

# Save depressed IDs
with open("/Users/maria/Documents/NSPN/docs/NSPN_IDs_depressed.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,len(depressed_id)):
        writer.writerow([i,depressed_id[i]])
