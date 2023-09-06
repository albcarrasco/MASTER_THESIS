# GET 'RELEVANT FOR THERAPY SELECTION' INFO

import pandas as pd
from pandas import notna
import gzip
import sys


# READ THE FILES
f1 = sys.argv[1] # VARIANTS file
f2 = sys.argv[2] # EVIDENCE file
f3 = sys.argv[3] # sample name

colsVAR = ["Gene","molecular_profile","AA_mutation","%_cfDNA_or_Amplification","MUTATION_SIGNIFICANCE_TIER"]
colsEVI = ["molecular_profile","therapies","therapy_interaction_type","evidence_direction","significance"] # get the columns of interest

VARdf = pd.read_csv(f1, sep="\t", usecols=colsVAR)
EVIdf = pd.read_csv(f2, sep="\t", usecols=colsEVI)

# KEEP VARIANTS IN A HIGH CMC TIER AND/OR PRESENT IN CIViC
valid_tiers = ["1","2"]
RELVARdf = VARdf[VARdf["MUTATION_SIGNIFICANCE_TIER"].isin(valid_tiers) | (VARdf["molecular_profile"] != ". .")]

# OUTPUT FILE
output = f"{f3}.ONCOZOOM_RELEVANT_SELECTION.tsv"
colsOUT = ["Gene","Alteration","%_cfDNA_or_Amplification","Resistance","Sensitivity/Response"]

# PROCESS VARIANTS FOR THERAPY SELECTION
def modify_aa(aa): # change format of the AA mutation to get a molecular profile of variants absenst in CIViC
    aa = aa.replace("p.","")
    return aa

def process_RELVAR(row): # extract field values from the rows in the RELEVANT VARIANTS dataframe
    gene = row["Gene"].iloc[0]
    alteration = modify_aa(row["AA_mutation"].iloc[0]) 
    amplification = row["%_cfDNA_or_Amplification"].iloc[0]
    return gene, alteration, amplification
    
def process_EVI(row): # process the rows in the EVIDENCE dataframe
    if row["therapy_interaction_type"] == "Combination":
        therapies = row["therapies"].replace(",", "+")
    #elif row["therapy_interaction_type"] == "Substitutes":
    #     therapies = row["therapies"].replace(",", "/")
    else:
        therapies = row["therapies"]
    return therapies

# check if there are relevant variants
if RELVARdf.empty:
    empty_msg = "No variants relevant for therapy selection were found.\nPlease check the full results for more detail."
    # write output
    with open(output, "w") as file:
        file.write("\t".join(colsOUT) + "\n")
        file.write(empty_msg) 
else:
    # create an empty dataframe
    OUTPUTdf = pd.DataFrame(columns=colsOUT)
    
    # RETRIEVE THERAPY INFO FOR EACH VARIANT
    # process VARIANTS.tsv
    RELVARdf = RELVARdf.copy()
    
    # create a new ID column called 'variant', this column will be identical to molecular_profile for those variants present in CIViC
    RELVARdf.loc[:,"variant"] = RELVARdf['Gene'] + " " + RELVARdf['AA_mutation'].apply(modify_aa) 
    var_list = (RELVARdf["variant"]).tolist() 
    
    for var in var_list: # for each variant
        vardf = RELVARdf[RELVARdf["variant"] == var].copy() # get only the row for that variant
        gene, alteration, amplification = process_RELVAR(vardf) # get info for that variant from the VARIANTS df
    
        # process EVIDENCE.tsv
        EVIsubset = EVIdf[EVIdf["molecular_profile"] == var].copy() # subset EVIDENCE.tsv for that specific variant
        # filter Resitance and Sensitivity therapies
        resistance_therapies = EVIsubset[(EVIsubset["significance"] == "Resistance") & (EVIsubset["evidence_direction"] != "Does not support")].copy()
        sensitivity_therapies = EVIsubset[(EVIsubset["significance"] == "Sensitivity/Response") & (EVIsubset["evidence_direction"] != "Does not support")].copy()
        # change format
        resistance_therapies.loc[:, "therapies"] = resistance_therapies.apply(process_EVI, axis=1)
        sensitivity_therapies.loc[:, "therapies"] = sensitivity_therapies.apply(process_EVI, axis=1)
        # remove nas, get unique values, sort (kinda messy lines)
        resistance_list = sorted(list(set((",".join(list(filter(notna,list(resistance_therapies["therapies"].tolist()))))).split(sep=","))))
        sensitivity_list = sorted(list(set((",".join(list(filter(notna,list(sensitivity_therapies["therapies"].tolist()))))).split(sep=","))))
        
        # replace empty lists (variants not found in CIViC) for 'None'
        for ls in resistance_list, sensitivity_list:
            if len(ls) == 1 and ls[0] == "":
                ls[0] = "None"
    
        # ADD info to OUTPUT
        OUTPUTdf.loc[len(OUTPUTdf)] = [gene, alteration, amplification, ", ".join(resistance_list), ", ".join(sensitivity_list)]
    # write
    OUTPUTdf.to_csv(output, sep="\t", index=False)
    
#######
