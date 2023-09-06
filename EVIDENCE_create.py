#!/usr/bin/env python3

# GET CIViC INFO for each variant

import pandas as pd
import gzip
import sys

# READ THE FILES
f1 = sys.argv[1] # CIViC evidence reference
f2 = sys.argv[2] # VARIANTS.tsv
f3 = sys.argv[3] # sample name

cols = ["molecular_profile","disease","therapies","therapy_interaction_type","evidence_type","evidence_direction","significance","rating","citation_id","evidence_statement"] # get the columns of interest
CIVICdf = pd.read_csv(f1, sep="\t", usecols=cols)
VARdf = pd.read_csv(f2, sep="\t", usecols=["molecular_profile"])

# FILTER OUT VARIANTS NOT PRESENT IN CIVIC
filtered_var = VARdf[~VARdf["molecular_profile"].str.contains("\.\s\.", na=False)]
var_list = filtered_var["molecular_profile"].tolist()

# FETCH DATA IN CIVIC REFERENCE FILE FOR THE FILTERED VARS
link_string = "https://pubmed.ncbi.nlm.nih.gov/"
stars_rating = {5.0:"*****",
                4.0:"****",
                3.0:"***",
                2.0:"**",
                1.0:"*"}

outputdf = pd.DataFrame(columns=cols) # empty dataframe to store the final output

for var in var_list:
    RESULTdf = CIVICdf[CIVICdf["molecular_profile"] == var].copy() # fetch data for that variant
    RESULTdf["citation_id"] = link_string + RESULTdf["citation_id"].astype(str) # create pubmed links
    RESULTdf["rating"] = RESULTdf["rating"].replace(stars_rating) # replace score rating for asterisks
    RESULTdf.sort_values(by="rating", ascending=False, inplace=True) # sort
    outputdf = pd.concat([outputdf, RESULTdf], ignore_index=True) # append to the output

# write the output dataframe to a file
output = f"{f3}.ONCOZOOM_EVIDENCE.tsv"
outputdf.to_csv(output, sep="\t", index=False) # write output
    
#####################
