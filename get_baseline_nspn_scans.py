# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 10:53:49 2016

@author: maria
"""

# Read subject's information from main NSPN database
import csv
import numpy as np

# Read main NSPN database
f = open('NSPN-datasetDownload-2016-02-09_MR.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data = []
data_cohort = []
data_gender = []
data_age = []
for row in csv_f:
    for l in range(0,len(row)):
        if (row[l] == 'iua_baseline_arm_1'):
            tmp = np.asarray(row)
            if (tmp[0] == tmp[l-1]):
                if (tmp[l+2] != ''):
                    if (tmp[l+2][0:2] == 'CB'):
                        tmp[l+2] = tmp[l+2][3::]
                    elif (tmp[l+2][0:2] == 'MQ'):
                        tmp[l+2] = tmp[l+2][2::]
                    data.append(tmp[(l-1):(l+9)])
                    data_gender.append(tmp[7])
                    data_cohort.append(tmp[2])
                    data_age.append(tmp[28])
            else:
                print "NSPN IDs don't match!!!"

# Check if information is read correctly
for i in range(0, len(data)):
    if (i == 0):
        c = len(data[i])
    if (i != 0):
	newc = len(data[i])
	if (newc != c):
	    print "Colum error!!!!!!"
	    print i


# data = np.asarray(data)
# with open("nspn_subjects_baseline_all_info.csv", "wb") as nf:
#     writer = csv.writer(nf)
#     for i in range(0, data.shape[0]):
#         writer.writerow(data[i])

data_cohort = np.asarray(data_cohort)
with open("nspn_subjects_baseline_cohort_info.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0, data_cohort.shape[0]):
        writer.writerow([i,data_cohort[i]])

data_gender = np.asarray(data_gender)
with open("nspn_subjects_baseline_gender_info.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0, data_gender.shape[0]):
        writer.writerow([i,data_gender[i]])

data_age = np.asarray(data_age)
with open("nspn_subjects_baseline_age_info.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0, data_age.shape[0]):
        writer.writerow([i,data_age[i]])
