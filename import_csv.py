# Read subject's information from main NSPN database
import csv
import numpy as np

# Read main NSPN database
f = open('NSPN-datasetDownload-2016-02-09_MR.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data = []
for row in csv_f:
    data.append(row)

# Check if information is read correctly
for i in range(0, len(data)):
    if (i == 0):
        c = len(data[i])
    if (i != 0):
	newc = len(data[i])
	if (newc != c):
	    print "Colum error!!!!!!"
	    print i

# Get fMRI subject IDs
g = open('column_1.csv','rU')
csv_g = csv.reader(g)
data_col1 = []
for r in csv_g:
    data_col1.append(r[0])

# Get IDs from main database
col_1 = []
for i in range(0,len(data)):
    col_1.append(data[i][0])

# Find fMRI subjects IDs on main database
fmri_subs = []
for j in range(0, len(data_col1)):
    idx = col_1.index(data_col1[j])
    fmri_subs.append(data[idx])

# Create new database with only fMRI subjects
data_fmri = np.asarray(fmri_subs)
data_fmri = data_fmri[::,[0,2,4,7,8,10,11,28,29,41,42,43,44,45,46,48,49,50,52,53,54,78,79,80,81,83,84,85,87,88,89,113,114,115,116,118,119,120,122,123,124]]
            
with open("nspn_fmri_subjects.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0, data_fmri.shape[0]):
        writer.writerow(data_fmri[i])
        
# Create new database with only baseline fMRI subjects        
f = open('nspn_fmri_subjects.csv', 'rU')
data_fmri_baseline = []
csv_f = csv.reader(f)
for row in csv_f:
    for l in range(0,len(row)):
        if (row[l] == 'iua_baseline_arm_1'):
            tmp = np.asarray(row)
            if (tmp[0] == tmp[l-1]):
                data_fmri_baseline.append(tmp[(l-1):(l+9)])
            else:
                print "NSPN IDs don't match!!!"
                
# Create new database with only baseline subjects from Gabriel data 
data_fmri_baseline = np.asarray(data_fmri_baseline)
with open("nspn_fmri_subjects_baseline.csv", "wb") as nf:
    writer = csv.writer(nf)
    for i in range(0, data_fmri_baseline.shape[0]):
        writer.writerow(data_fmri_baseline[i])
        

# Get fMRI subject IDs
# Read main NSPN database
f = open('GabrielIDs.csv', 'rU')
csv_f = csv.reader(f)
header = next(csv_f)
data_gabriel_scans = []
data_gabriel_ids = []
for row in csv_f:
    if (row[1] == 'iua_baseline_arm_1'):
        data_gabriel_scans.append(row[3])
        data_gabriel_ids.append(row[0])
    

## Get fMRI subject IDs from Gabriel
#g = open('IDs_T1s.csv','rU')
#csv_g = csv.reader(g)
#data_col1 = []
#for r in csv_g:
#    data_col1.append(r[0])
#
## Get IDs from main database
#col_1 = []
#for i in range(0,len(data)):
#    col_1.append(data[i][0])
#
## Find fMRI subjects IDs on main database
#fmri_subs = []
#for j in range(0, len(data_col1)):
#    idx = col_1.index(data_col1[j])
#    fmri_subs.append(data[idx])
#
## Create new database with only fMRI subjects
#data_fmri = np.asarray(fmri_subs)
#data_fmri = data_fmri[::,[0,2,4,7,8,10,11,28,29,41,42,43,44,45,46,48,49,50,52,53,54,78,79,80,81,83,84,85,87,88,89,113,114,115,116,118,119,120,122,123,124]]
            
#with open("nspn_fmri_subjects_gabriel.csv", "wb") as nf:
#    writer = csv.writer(nf)
#    for i in range(0, data_fmri.shape[0]):
#        writer.writerow(data_fmri[i])
#        
## Create new database with only baseline fMRI subjects        
#f = open('nspn_fmri_subjects_gabriel.csv', 'rU')
#data_fmri_baseline_gab = []
#csv_f = csv.reader(f)
#for row in csv_f:
#    for l in range(0,len(row)):
#        if (row[l] == 'iua_baseline_arm_1'):
#            tmp = np.asarray(row)
#            if (tmp[0] == tmp[l-1]):
#                data_fmri_baseline_gab.append(tmp[(l-1):(l+9)])
#            else:
#                print "NSPN IDs don't match!!!"
#                
## Create new database with only baseline subjects from Gabriel data 
#data_fmri_baseline_gab = np.asarray(data_fmri_baseline_gab)
#with open("nspn_fmri_subjects_baseline_gabriel.csv", "wb") as nf:
#    writer = csv.writer(nf)
#    for i in range(0, data_fmri_baseline_gab.shape[0]):
#        writer.writerow(data_fmri_baseline_gab[i])
#        
# Check fMRI scanner IDs with Gabriel list
data_gabriel_scans = np.asarray(data_gabriel_scans)
g = 0
t = 0
for r in data_gabriel_ids:
    if (r in list(data_fmri_baseline[::,0])):
        idx = list(data_fmri_baseline[::,0]).index(r)
        if (data_fmri_baseline[idx][3] == data_gabriel_scans[g]):
            print "fMRI scan ID does match T1 ID!!!"
            t = t + 1
    else:
        print "Subject not in fMRI list!!"
        print str(r)
    g = g+1



        

