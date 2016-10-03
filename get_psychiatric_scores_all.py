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
f = open('/Users/maria/Documents/NSPN/docs/nspn_subjects_baseline.csv', 'rU')
csv_f = csv.reader(f)
all_baseline_scans = []
all_baseline_ids = []
for row in csv_f:
    all_baseline_scans.append(row[2])
    all_baseline_ids.append(row[0])

# Get cohort
f = open('/Users/maria/Documents/NSPN/docs/nspn_subjects_baseline_cohort_info.csv', 'rU')
csv_f = csv.reader(f)
all_cohort = []
for row in csv_f:
    all_cohort.append(row[1])

# Get psychiatric scores
f = open('/Users/maria/Documents/NSPN/docs/dat4JDMMa.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
all_psy = []
all_psy_ids = []
for row in csv_f:
    all_psy.append(row[1::])
    all_psy_ids.append(row[0])

# Get scans
f = open('/Users/maria/Documents/NSPN/docs/NSPN_scans_all.csv', 'rU')
csv_f = csv.reader(f)
ids_mri = []
depressed_id = []
nspn_ids_depressed = []
data_psy_dep = []
k = 0
for row in csv_f:
    indx = all_baseline_scans.index(row[0])
    nspn_id = all_baseline_ids[indx]
    cohort_id = all_cohort[indx]
    k = k + 1
    if cohort_id == 'Depression':
        depressed_id.append(k)
        nspn_ids_depressed.append(nspn_id)
        if (all_psy_ids.count(nspn_id) >=1):
            indx_psy = all_psy_ids.index(nspn_id)
            tmp = all_psy[indx_psy]
            data_psy_dep.append(tmp)
        else:
            data_psy_dep.append(['NA','NA','NA','NA'])

        print nspn_id
    ids_mri.append(nspn_id)

# Read questionnaire data
f = open('/Users/maria/Documents/NSPN/docs/NSPN_HPQ1_processed.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_hpq = []
data_hpq_ids = []
for row in csv_f:
    data_hpq.append(row)
    data_hpq_ids.append(row[0])

# Read baseline subjects
f = open('/Users/maria/Documents/NSPN/docs/UCHANGE_CompleteCohort_BaselineOnly_20150722.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_psy = []
psy_ids = []
baseline_vars = []
for row in csv_f:
    indx_hpq = data_hpq_ids.index(row[0])
    tmp = data_hpq[indx_hpq]
    baseline_vars.append(tmp[0])
    if (all_psy_ids.count(row[0]) >=1):
        indx_psy = all_psy_ids.index(row[0])
        tmp = all_psy[indx_psy]
        data_psy.append(tmp)
        psy_ids.append(all_psy_ids[indx_psy])
    else:
        print row[0]
        data_psy.append(['NA','NA','NA','NA'])
        psy_ids.append(row[0])

# Save variables data - baseline cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_psy_baseline.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_mri)[0]):
        if (psy_ids.count(ids_mri[i])) > 0:
            indx = psy_ids.index(ids_mri[i])
            writer.writerow(data_psy[indx])

# Save variables data - depressed cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_psy_depressed.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_mri)[0]):
        if (nspn_ids_depressed.count(ids_mri[i])) > 0:
            indx = nspn_ids_depressed.index(ids_mri[i])
            writer.writerow(data_psy_dep[indx])
