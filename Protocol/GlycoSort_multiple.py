# Import packages
import pandas as pd 
import sys
import re 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import scipy
import math
import getopt

from statsmodels.sandbox.stats.multicomp import multipletests
from functools import partial, reduce

pd.options.mode.chained_assignment = None
def getArgs(argv):
    bionic = ""
    maxQuant = ""
    experiment = "Experiment"
    argHelp = "{0} -b <bionic file> -m <maxquant file> -o <output file name> -p <plot decision>".format(argv[0])
	

    if len(argv) == 1:
        print(argHelp)
        sys.exit(2)

    try:
        opts, args = getopt.getopt(argv[1:], "h:b:m:o:", ["help", "bionic=", "maxquant=", "output="])
        #print(opts)
    except:
        print(argHelp)
        sys.exit(2)
    
    for opt, arg in opts:
        #print(arg)
        if opt in ["-h", "--help"]:
            print(argHelp)
            sys.exit(2)
        elif opt in ["-b", "--bionic"]:
            bionic = arg
        elif opt in ["-m", "--maxquant"]:
            maxQuant = arg
        elif opt in ["-o", "--output"]:
            experiment = arg


    return bionic, maxQuant, experiment

if __name__ == "__main__":
    getArgs(sys.argv)

bionic, maxQuant, experiment = getArgs(sys.argv)

#create folder for plots 
parent = os.getcwd()
directory = experiment
gly = "Glyco"
sit = "Site"
path = os.path.join(parent, directory) 
pathGly = os.path.join(path, gly)
pathSite = os.path.join(path, sit)
print(path)

try:
	os.mkdir(path)
	os.mkdir(pathGly)
	os.mkdir(pathSite)
except (FileExistsError):
    pass
except Exception as e:
    raise e


bioDf = pd.read_excel(bionic)
maxQuantDf = pd.read_excel(maxQuant)

#Lists containing columns that will be used to filter the dataset and create subsets containing only the desired information
desired = ["ProteinName", "Position", "Sequence", "Glycopeptide", "VariableFinalModNameList", "CompositionList", "PEP", "Charge", 
"Cleavage", "ObservedMH", "CalcMH", "ScanTimeList", "TotalIonCurrent", "fullIonCurrent", "icPercentageSite", "icPercentageGly"]
desired2 = ["ProteinName", "Position", "Glycopeptide", "TotalIonCurrent","icPercentageGly", "fullIonCurrent"]
desired3 = ["ProteinName", "Position", "Glycopeptide","TotalIonCurrent", "icPercentageSite", "fullIonCurrent"]

#Function for Total Ion Current 
def getTIC(bionic, maxQuant):
    """
    This function will compare the Scan Number from
    the bionic and maxQuant files, then will retrieve
    the Total Ion Current from the maxQuant file.
    Paramenters:
    bionic = the bionic file
    maxQuant = the maxQuant file
    """
    for i in bionic.index:
        #print(maxQuant['Scan_number'])
        scanNumber = bionic.loc[i,"ScanNumber"]
        if maxQuant.loc[maxQuant['Scan_number'] == scanNumber, 'Total_ion_current'].empty:
            continue
        else:
            bionic.loc[i, "TotalIonCurrent"] = maxQuant.loc[maxQuant['Scan_number'] == scanNumber, 'Total_ion_current'].iloc[0]
    
    bionic['TotalIonCurrent'] = bionic.groupby(['Glycopeptide'])['TotalIonCurrent'].transform("sum")
    newDf = bionic.drop_duplicates(subset=['Glycopeptide'])
    return newDf

d1 = {}
d2 = {}
d3 = {}
for x in range(1,len(set(list(bioDf['type'])))+1):
    d1["dataSeparatedTRE{0}".format(x)] = "treatment {0}".format(x)
    d2["dataGlycoMergedTRE{0}".format(x)] = "treatment {0}".format(x)
    d3["dataSiteTRE{0}".format(x)] = "treatment {0}".format(x)

nestDici = {}
typeCheck = bioDf["type"][1]
tempd1 = {}
tempd2 = {}
tempd3 = {}
for i in range(len(bioDf['file'])):

    fileName = bioDf['name'][i]
    #print(fileName)
    fileType = bioDf['type'][i]
    #print(fileType)
    if fileType != typeCheck:
        tempd1 = {}
        tempd2 = {}
        tempd3 = {}
    
    typeCheck = fileType
    bDf = pd.read_csv(bioDf['file'][i], sep = ",")
    bDf.dropna(how = "all", inplace = True)
    #print(bDf)
    mqDf = pd.read_csv(maxQuantDf['file'][i], sep = "\t")
    mqDf.columns = mqDf.columns.str.replace(' ', '_')
    mqDf.dropna(how = "all", inplace = True)
    bDf['CompositionList'].replace('', np.nan, inplace=True)
    bDf.dropna(subset=['CompositionList'], inplace=True)
    bDf=bDf.loc[bDf['PEP'] <= 0.00099]
    bDf[~bDf.ProteinName.str.contains("Common contaminant")]
    #Creating glycopeptide column by merging the Sequence and CompositionList columns
    bDf['Sequence'] = bDf['Sequence'].str.replace("[\(\[].*?[\)\]]", "", regex = True)
    bDf['Glycopeptide'] = bDf['Sequence'] + "." + bDf['CompositionList']
    bDf['VariableFinalModNameList'] = bDf['VariableFinalModNameList'].str.replace(r"\;","_", regex = True)
    bDf['Position'] = bDf['Position'].astype(int) + bDf['VariableFinalModNameList'].str.split(pat = "N", expand = True)[1].str.split(pat = "(", expand = True)[0].astype(int) - 1
        
    #Separating the Scan Number from the ScanNumberList column
    bDf['ScanNumber'] = pd.to_numeric(bDf['ScanNumberList'].str.split("scan=", expand = True)[1])

    #Getting the merged total ion current
    final = getTIC(bDf, mqDf)
    #Calculating the percentages of each site per protein and each glycopeptide per site per protein
    ticCol = final.columns[-1]
    final['fullIonCurrent'] = final.groupby(["ProteinName"])[ticCol].transform("sum")
    final['icPercentageSite'] = (final.groupby(["ProteinName", "Position"])[ticCol].transform("sum")/final.groupby(["ProteinName"])[ticCol].transform("sum"))*100
    final['icPercentageGly'] = (final.groupby(["ProteinName", "Position","Glycopeptide"])[ticCol].transform("sum")/final.groupby(["ProteinName", "Position"])[ticCol].transform("sum"))*100

    #Selecting specific columns, creating sub dataframes and renaming columns based on sample IDs
    finalSelected = final[desired]
    finalSelected2 = final[desired2]
    finalSelected3 = final[desired3]
    name = str("TotalIonCurrent_" + fileName)
    fullIon = str("fullIonCurrent_" + fileName)
    icPerSite = str("icPercentageSite_" + fileName)
    icPerGly = str("icPercentageGly_" + fileName)
    finalRenamed = finalSelected.rename(columns = {"TotalIonCurrent": name, "fullIonCurrent": fullIon, "icPercentageSite": icPerSite, "icPercentageGly": icPerGly})
    finalRenamed2 = finalSelected2.rename(columns = {"TotalIonCurrent": name, "icPercentageGly": icPerGly, "fullIonCurrent": fullIon})
    finalRenamed3 = finalSelected3.rename(columns = {"TotalIonCurrent": name, "icPercentageSite": icPerSite, "fullIonCurrent": fullIon})
    
    tmp1 = list(d1.keys())[list(d1.values()).index(fileType)]
    tmp2 = list(d2.keys())[list(d2.values()).index(fileType)]
    tmp3 = list(d3.keys())[list(d3.values()).index(fileType)]
    
    tempd1.update({fileName + "_" + fileType: finalRenamed})
    tempd2.update({fileName + "_" + fileType: finalRenamed2})
    tempd3.update({fileName + "_" + fileType: finalRenamed3.drop_duplicates(subset=['Position'])})
    
    nestDici.update({tmp1: tempd1,
               tmp2: tempd2,
               tmp3: tempd3})

tempGlyco = partial(pd.merge, on=['ProteinName','Glycopeptide', 'Position'], how='outer')
tempSite = partial(pd.merge, on=['ProteinName', 'Position'], how='outer')

#Formating and saving files
count = 1
filDi = {}
z = 1
for i in nestDici.keys():
    if count > 3:
        z = z + 1
        count = 1
    if "Glyco" in i:
        filDi.update({"dfGlycopepTRE{0}".format(z): reduce(tempGlyco, nestDici[i].values())})
        filDi["dfGlycopepTRE{0}".format(z)].to_csv(os.path.join(path,"glycopepCombinedTRE{0}.csv".format(z)), index = False)
    if "Site" in i: 
        filDi.update({"dfProtSiteTRE{0}".format(z): reduce(tempSite, nestDici[i].values())})
        filDi["dfProtSiteTRE{0}".format(z)].to_csv(os.path.join(path,"proteinSiteCombinedTRE{0}.csv".format(z)), index = False)
    else:
        for a in nestDici["dataSeparatedTRE{0}".format(z)]:
            nestDici["dataSeparatedTRE{0}".format(z)][a].to_csv(os.path.join(path,str(a) + "_processed.csv"), index = False)
    count = count + 1
