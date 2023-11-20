"""
    Project :  variants intersection df for VHL
    Description: merge raw SGE data and the UKB allele table.
    Name : Chloé Terwagne
    date : 19th Jan 2023
    Python version 3.10
"""
# IMPORT ---------------------------------------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import datetime

pd.set_option('display.width', 3000)
pd.set_option('display.max_columns', 3000)
pd.set_option("display.max_rows", None)

# CONSTANT -------------------------------------------------------------------------------------------------------------
DATE = str(datetime.date.today()).replace('-', '_')
GENE = "vhl"
PATH_ASSAY = "VHL_SGE_ST1.1_no1c_20230911.csv"

vhl_sge_df = pd.read_csv(PATH_ASSAY, low_memory=False)
# vhl_sge_df.sge_region != "exon 1c"
vhl_sge_df.rename({'Chrom': 'chr', 'Ref': 'ref', 'conseq': 'consequence'},
                  axis=1, inplace=True)
vhl_sge_df["variant_id"] = vhl_sge_df['chr'].astype(str) + '_' + vhl_sge_df['hg38_pos'].astype(str) + '_' + vhl_sge_df['ref'] + '_' + vhl_sge_df['alt']
vhl_sge_df = vhl_sge_df[["variant_id", "chr", "hg38_pos", "ref", "alt",'nAA', 'sge_region','protPos', 'cHGVS', 'pHGVS',
     "function_score_final", "q_value", "tier_class", "rna_score", "rna_score_day_20", "delta_rna",
     "consequence",'clinvar_simple',
     "VARITY_R", "REVEL", "CADD.phred", 'max_spliceAI', 'Cancer_type_single']]

vhl_sge_df = vhl_sge_df.replace({"consequence": {"SYNONYMOUS": "Synonymous", "NON_SYNONYMOUS": "Non synonymous",
                                 "STOP_GAINED": "Stop gained", "SPLICE_SITE": "Splice site", "INTRONIC": "Intronic",
                                 "CANONICAL_SPLICE": "Canonical splice", "STOP_LOST": "Stop lost", "3PRIME_UTR": "3' UTR"}})

vhl_sge_df["clinvar_simple"] = vhl_sge_df["clinvar_simple"].replace('absent', 'Absent ')
vhl_sge_df["Cancer_type_single"] = vhl_sge_df["Cancer_type_single"].replace(np.nan, 'Absent')

# add position columns for overview graph
g_pos = (max(vhl_sge_df.function_score_final) - min(vhl_sge_df.function_score_final))/3
print(g_pos)
c_pos = g_pos*2
print(c_pos)
vhl_sge_df['alt_pos'] = [x for x in vhl_sge_df['alt'].map({"A": max(vhl_sge_df.function_score_final), "C": min(vhl_sge_df.function_score_final)+c_pos, "G": min(vhl_sge_df.function_score_final)+g_pos, "T": min(vhl_sge_df.function_score_final)})]
vhl_sge_df['ref_pos'] = [x for x in vhl_sge_df['ref'].map({"A": max(vhl_sge_df.function_score_final), "C": min(vhl_sge_df.function_score_final)+c_pos, "G": min(vhl_sge_df.function_score_final)+g_pos, "T": min(vhl_sge_df.function_score_final)})]


print("\n MISSENSE MEAN")

# Filter rows where 'consequence' is "Non synonymous"
vhl_missense_only = vhl_sge_df[vhl_sge_df['consequence'] == "Non synonymous"].copy()
vhl_missense_only_rnagreaterthanmin2 = vhl_sge_df[(vhl_sge_df['consequence'] == "Non synonymous")& (vhl_sge_df['rna_score'] >= -2)].copy()
print('vhl_missense_only: ', vhl_missense_only.shape)
print("vhl_missense_only_rnagreaterthanmin2:",  vhl_missense_only_rnagreaterthanmin2.shape)

# Calculate the mean based on 'protPos' and store it in a new column in the original DataFrame
vhl_sge_df['average_fs_missense_at_aa'] = vhl_missense_only.groupby('protPos')['function_score_final'].transform('mean')
vhl_sge_df['average_fs_missense_at_aa_rna'] = vhl_missense_only_rnagreaterthanmin2.groupby('protPos')['function_score_final'].transform('mean')

#vhl_sge_df.loc[vhl_sge_df['consequence'] == "Non synonymous", 'average_fs_missense_at_aa'] = vhl_missense_only.groupby('protPos')['function_score_final'].transform('mean')
vhl_missense_only['average_fs_missense_at_aa'] = vhl_missense_only.groupby('protPos')['function_score_final'].transform('mean')
vhl_missense_only_rnagreaterthanmin2['average_fs_missense_rna'] = vhl_missense_only_rnagreaterthanmin2.groupby('protPos')['function_score_final'].transform('mean')

# Get txt mapping file missense only
vhl_missense_only = vhl_missense_only[['protPos', 'average_fs_missense_at_aa']]
vhl_missense_only.drop_duplicates(inplace=True)
vhl_missense_only['protPos'] = [int(x) for x in list(vhl_missense_only['protPos'])]
vhl_missense_only.to_csv("Average_FS_at_AA_missense_only.txt", sep='\t', header=None, index=None)

# Get txt mapping file missense only with RNA greater than -2
vhl_missense_only_rnagreaterthanmin2 = vhl_missense_only_rnagreaterthanmin2[['protPos', 'average_fs_missense_rna']]
vhl_missense_only_rnagreaterthanmin2.drop_duplicates(inplace=True)
vhl_missense_only_rnagreaterthanmin2['protPos'] = [int(x) for x in list(vhl_missense_only_rnagreaterthanmin2['protPos'])]
vhl_missense_only_rnagreaterthanmin2.to_csv("Average_FS_at_AA_missense_only_rna_not_below_min2.txt", sep='\t', header=None, index=None)


vhl_sge_df['average_fs_at_aa'] = vhl_sge_df.groupby('protPos')['function_score_final'].transform('mean')



print(vhl_sge_df.head())
print(vhl_sge_df.shape)
print(len(vhl_sge_df["variant_id"]))
print(len(vhl_sge_df["variant_id"].unique()))
vhl_sge_df.to_csv("vhl_preprocess_df.csv")
print(vhl_sge_df.tier_class.unique())
print(vhl_sge_df.consequence.unique())
print(vhl_sge_df.clinvar_simple.unique())

print(vhl_sge_df.sge_region.unique())

print("min max average_fs_missense_at_aa_rna: ")
print(min(vhl_sge_df.average_fs_missense_at_aa_rna), max(vhl_sge_df.average_fs_missense_at_aa_rna))


print("min max function score: ")
print(min(vhl_sge_df.function_score_final), max(vhl_sge_df.function_score_final))
print(min(vhl_sge_df.delta_rna), max(vhl_sge_df.delta_rna))

exon_lim_dict = {}
for exon in vhl_sge_df.sge_region.unique():
    print(exon, end=" ")
    df_temp = vhl_sge_df[vhl_sge_df.sge_region == exon]
    print(min(df_temp.hg38_pos), max(df_temp.hg38_pos))
    exon_lim_dict[exon] = [min(df_temp.hg38_pos), max(df_temp.hg38_pos)]
# Create index
prev = 0
i = 0
list_indx=[]
for pos in vhl_sge_df["hg38_pos"]:
    if pos!=prev:
        i+=1
    list_indx.append(i)
    prev=pos

vhl_sge_df["index"]=list_indx
# create gap between exon 1p 2 et 3 for better vizualisation -> still not at scale
# exon1p idx +10, exon2 idx +30, exon3a and 3b idx +30
vhl_sge_df.loc[vhl_sge_df.sge_region == 'exon 1p', 'index'] = vhl_sge_df[vhl_sge_df.sge_region == 'exon 1p']['index'] + 30
vhl_sge_df.loc[vhl_sge_df.sge_region == 'exon 2', 'index'] = vhl_sge_df[vhl_sge_df.sge_region == 'exon 2']['index'] + 60
vhl_sge_df.loc[vhl_sge_df.sge_region == 'exon 3a', 'index'] = vhl_sge_df[vhl_sge_df.sge_region == 'exon 3a']['index'] + 90
vhl_sge_df.loc[vhl_sge_df.sge_region == 'exon 3b', 'index'] = vhl_sge_df[vhl_sge_df.sge_region == 'exon 3b']['index'] + 90

print(vhl_sge_df.head(50))
print(exon_lim_dict)
print("TGCA")
print(min(vhl_sge_df.function_score_final), ",",  min(vhl_sge_df.function_score_final)+g_pos, ",",min(vhl_sge_df.function_score_final)+c_pos,",", max(vhl_sge_df.function_score_final))

vhl_sge_df.to_csv("vhl_preprocess_df.csv")