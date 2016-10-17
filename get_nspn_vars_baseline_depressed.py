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
all_baseline_centre = []
all_baseline_ids = []
for row in csv_f:
    all_baseline_scans.append(row[2])
    all_baseline_centre.append(row[1])
    all_baseline_ids.append(row[0])

# Get gender
f = open('/Users/maria/Documents/NSPN/docs/nspn_subjects_baseline_gender_info.csv', 'rU')
csv_f = csv.reader(f)
all_gender = []
for row in csv_f:
    all_gender.append(row[1])

# Get cohort
f = open('/Users/maria/Documents/NSPN/docs/nspn_subjects_baseline_cohort_info.csv', 'rU')
csv_f = csv.reader(f)
all_cohort = []
for row in csv_f:
    all_cohort.append(row[1])

# Get cohort
f = open('/Users/maria/Documents/NSPN/docs/nspn_subjects_baseline_age_info.csv', 'rU')
csv_f = csv.reader(f)
all_age = []
for row in csv_f:
    all_age.append(row[1])

# Get scans
f = open('/Users/maria/Documents/NSPN/docs/NSPN_scans_all.csv', 'rU')
csv_f = csv.reader(f)
ids_mri = []
ids_centre = []
ids_gender = []
ids_cohort = []
ids_age = []
depressed_id = []
nspn_ids_depressed = []
k = 0
for row in csv_f:
    indx = all_baseline_scans.index(row[0])
    nspn_id = all_baseline_ids[indx]
    centre_id = all_baseline_centre[indx]
    gender_id = all_gender[indx]
    cohort_id = all_cohort[indx]
    age_id = all_age[indx]
    k = k + 1
    if cohort_id == 'Depression':
        depressed_id.append(k)
        nspn_ids_depressed.append(nspn_id)
        print k
        print cohort_id
    ids_mri.append(nspn_id)
    ids_centre.append(centre_id)
    ids_gender.append(gender_id)
    ids_cohort.append(cohort_id)
    ids_age.append(age_id)

# Read HPQ1
f = open('/Users/maria/Documents/NSPN/docs/NSPN_HPQ1_processed_110816.csv', 'rU')
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

# Read SCID
f = open('/Users/maria/Documents/NSPN/docs/NSPN_SCID_depressed.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_scid = []
data_scid_ids = []
for row in csv_f:
    data_scid.append(row)
    data_scid_ids.append(row[0])

# Get weight and height
f = open('/Users/maria/Documents/NSPN/docs/NSPN_weight_height_waist_processed.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_wh = []
data_wh_ids = []
for row in csv_f:
    data_wh.append(row)
    data_wh_ids.append(row[0])

# Read depressed subjects
data_depressed = []
data_wh_depressed = []
scid_data = []
for row in nspn_ids_depressed:
    indx_hpq = data_hpq_ids.index(row)
    indx_iua = data_iua_ids.index(row)
    indx_wh = data_wh_ids.index(row)
    indx_scid = data_scid_ids.index(row)
    tmp = data_hpq[indx_hpq] + data_iua[indx_iua][1::]
    tmp_wh = data_wh[indx_wh]
    tmp_scid = data_scid[indx_scid]
    scid_data.append(tmp_scid)
    data_depressed.append(tmp)
    data_wh_depressed.append(tmp_wh)

# Save variables data - depressed cohort
# ids_good = []
# mri_good = []
# gender_good = []
# cohort_good = []
# age_good = []
# with open("/Users/maria/Documents/NSPN/docs/NSPN_vars_depressed.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0,np.shape(ids_mri)[0]):
#         if (nspn_ids_depressed.count(ids_mri[i])) > 0:
#             indx = nspn_ids_depressed.index(ids_mri[i])
#             writer.writerow(data_depressed[indx])
#             ids_good.append(i+1)
#             mri_good.append(ids_centre[i])
#             gender_good.append(ids_gender[i])
#             cohort_good.append(ids_cohort[i])
#             age_good.append(ids_age[i])

# Save variables data - depressed cohort
with open("/Users/maria/Documents/NSPN/docs/NSPN_SCID_baseline_depressed.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0,np.shape(ids_mri)[0]):
        if (nspn_ids_depressed.count(ids_mri[i])) > 0:
            indx = nspn_ids_depressed.index(ids_mri[i])
            writer.writerow(scid_data[indx])

# Save weight and height data - depressed cohort
# with open("/Users/maria/Documents/NSPN/docs/NSPN_weight_height_waist_depressed.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0,np.shape(ids_mri)[0]):
#         if (nspn_ids_depressed.count(ids_mri[i])) > 0:
#             indx = nspn_ids_depressed.index(ids_mri[i])
#             writer.writerow(data_wh_depressed[indx])

# # Save IDs - depressed cohort
# with open("/Users/maria/Documents/NSPN/docs/NSPN_MRIids_depressed.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0,np.shape(ids_good)[0]):
#         writer.writerow([ids_good[i],ids_good[i]])
#
# # Save centre information - depressed cohort
# with open("/Users/maria/Documents/NSPN/docs/NSPN_gender_bin_depressed.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0,np.shape(ids_good)[0]):
#         if gender_good[i] == 'Female':
#             writer.writerow([i,0])
#         else:
#             writer.writerow([i,1])
#
# # Save gender information - depressed cohort
# with open("/Users/maria/Documents/NSPN/docs/NSPN_cohort_depressed.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0,np.shape(ids_good)[0]):
#         writer.writerow([i,cohort_good[i]])
#
# # Save gender information - depressed cohort
# with open("/Users/maria/Documents/NSPN/docs/NSPN_MRIcentre_depressed.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0,np.shape(ids_good)[0]):
#         writer.writerow([i,mri_good[i]])
#
# # Save gender information - depressed cohort
# with open("/Users/maria/Documents/NSPN/docs/NSPN_age_depressed.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0,np.shape(ids_good)[0]):
#         writer.writerow([i,age_good[i]])
