#!/usr/bin/env python3

# GET CMC TIER, SOMATIC FREQ AND COSMIC URL FOR CODING SOMATIC VARIANTS

import pandas as pd
import gzip
import sys

# READ THE FILES
f1 = sys.argv[1] # COSMIC CMC ref file
f2 = sys.argv[2] # COSMIC FREQ ref file
f3 = sys.argv[3] # VARIANTS generated file

cols = ["LEGACY_MUTATION_ID","MUTATION_SIGNIFICANCE_TIER"] # get the columns of interest
CMCdf = pd.read_csv(f1, sep="\t", usecols=cols)

FREQdf = pd.read_csv(f2, sep="\t", usecols=["LEGACY_MUTATION_ID","FREQ"])

VARdf = pd.read_csv(f3, sep="\t")

# create the url column 
# follow the COSMIC recommendations for creating a link
VARdf["MUTATION_URL"] = "https://cancer.sanger.ac.uk/cosmic/search?genome=37&q=" + VARdf["LEGACY_MUTATION_ID"].astype(str) 

# MERGE DATAFRAMES
MERGEDdfA = VARdf.merge(FREQdf)
MERGEDdfB = MERGEDdfA.merge(CMCdf)
MERGEDdfB["%_cfDNA_or_Amplification"] = MERGEDdfB["%_cfDNA_or_Amplification"].apply(lambda x: round(x*100,1))

# SAVE THE MERGE TO A NEW TSV FILE
MERGEDdfB.to_csv(f"VARIANTS.tsv", sep="\t", index=False)
