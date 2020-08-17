import pandas as pd
import os.path
import os
import allel

#Specify working directory
working_dir = "/home/nick_rose/nick/Downloads/vcf/"

#Create list of file names within direcory
for root, dirs, files in os.walk(working_dir):
    file_list = []
    for filename in files:
        if filename.endswith('.vcf'):
            file_list.append(os.path.join(root, filename))

#Create master file
df_master = pd.DataFrame()

#Extract data from VCF files and Append to Master
for f in file_list:
    #Create header of loci
    df = allel.vcf_to_dataframe(f)
    df["LOCI"] = df["CHROM"] + "_" + df["POS"].astype(str)
    df = df[['LOCI']]
    df = df.T
    new_header = df.iloc[0]
    #Define Callset Data
    callset = allel.read_vcf(f, fields = ['calldata/REPCN'])
    callset = callset['calldata/REPCN']
    #Create list of Sample names (Temporary naming system, only works for samples 7n long)
    sam = f[len(working_dir):(len(working_dir) + 7)]
    df2 = pd.DataFrame(callset, columns=[sam])
    #Insert REPCN data
    df2 = df2.T
    df2.columns = (new_header)
    #Append file to master
    df_master = pd.concat([df_master,df2])

#Create CSV file
df_master = df_master.reset_index()
df_master.to_csv('/home/nick_rose/nick/Downloads/output/df_master.csv', index = False)
print('Done')
