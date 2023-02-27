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
    plotDecision = "relevant"
    threshNaN = 0.7
    argHelp = "{0} -b <bionic file> -m <maxquant file> -o <output file name> -p <plot decision>".format(argv[0])
	

    if len(argv) == 1:
        print(argHelp)
        sys.exit(2)

    try:
        opts, args = getopt.getopt(argv[1:], "h:b:m:o:p:t:", ["help", "bionic=", "maxquant=", "output=", "plot_decision=", "thresh_NaN="])
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
        elif opt in ["-p", "--plot"]:
            plotDecision = arg
        elif opt in ["-t", "--thresh"]:
            threshNaN = arg

    return bionic, maxQuant, experiment, plotDecision, threshNaN

if __name__ == "__main__":
    getArgs(sys.argv)

bionic, maxQuant, experiment, plotDecision, threshNaN = getArgs(sys.argv)

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
        #print(maxQuant['Scan'])
        scanNumber = bionic.loc[i,"ScanNumber"]
        bionic.loc[i, "TotalIonCurrent"] = maxQuant.loc[maxQuant['Scan_number'] == scanNumber, 'Total_ion_current'].iloc[0]
    
    bionic['TotalIonCurrent'] = bionic.groupby(['Glycopeptide'])['TotalIonCurrent'].transform("sum")
    newDf = bionic.drop_duplicates(subset=['Glycopeptide'])
    return newDf

dataSeparatedCTR = {}
dataGlycoMergedCTR = {}
dataSiteCTR = {}
dataSeparatedTRE = {}
dataGlycoMergedTRE = {}
dataSiteTRE = {}
for i in range(len(bioDf['file'])):

    fileName = bioDf['name'][i]
    fileType = bioDf['type'][i]

    bDf = pd.read_csv(bioDf['file'][i], sep = ",")
    bDf.dropna(how = "all", inplace = True)
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
    if fileType == "control":
        dataSeparatedCTR[fileName + "_" + fileType] = finalRenamed
        dataGlycoMergedCTR[fileName + "_" + fileType] = finalRenamed2
        dataSiteCTR[fileName + "_" + fileType] = finalRenamed3.drop_duplicates(subset=['Position'])
    else:
        dataSeparatedTRE[fileName + "_" + fileType] = finalRenamed
        dataGlycoMergedTRE[fileName + "_" + fileType] = finalRenamed2
        dataSiteTRE[fileName + "_" + fileType] = finalRenamed3.drop_duplicates(subset=['Position'])

#Merging dataframes by group
tempGlyco = partial(pd.merge, on=['ProteinName','Glycopeptide', 'Position'], how='outer')
tempSite = partial(pd.merge, on=['ProteinName', 'Position'], how='outer')

dfGlycopepCTR = reduce(tempGlyco, dataGlycoMergedCTR.values())
dfGlycopepTRE = reduce(tempGlyco, dataGlycoMergedTRE.values())
dfProtSiteCTR = reduce(tempSite, dataSiteCTR.values())
dfProtSiteTRE= reduce(tempSite, dataSiteTRE.values())

#Saving files
dfGlycopepCTR.to_csv(os.path.join(path, "glycopepCombinedCTR.csv"), index = False)
dfProtSiteCTR.to_csv(os.path.join(path,"proteinSiteCombinedCTR.csv"), index = False)
dfGlycopepTRE.to_csv(os.path.join(path,"glycopepCombinedTRE.csv"), index = False)
dfProtSiteTRE.to_csv(os.path.join(path,"proteinSiteCombinedTRE.csv"), index = False)

for i in dataSeparatedCTR:
	dataSeparatedCTR[i].to_csv(os.path.join(path, str(i) + "_processed.csv"), index = False)

for i in dataSeparatedTRE:
	dataSeparatedTRE[i].to_csv(os.path.join(path, str(i) + "_processed.csv"), index = False)



#Statistics using the TIC value
#Remove the percentages
def arrangeData(df):
    """
    This function keeps the
    rows that has atleast 70% of values (not NaNs).
    """
    
    df_stat = df[df.columns.drop(list(df.filter(regex='icPercentageGly')))].copy()
    df_stat2 = df_stat[df_stat.columns.drop(list(df_stat.filter(regex='fullIonCurrent')))].copy()
    df_stat2.dropna(axis=0, thresh=(int(len(df_stat2.filter(regex=r'TotalIonCurrent').columns)*float(threshNaN))), subset = list(df_stat2.filter(regex=r'TotalIon').columns), inplace = True)
    
    return df_stat2

dfGlycopepCTR_stat = arrangeData(dfGlycopepCTR)
dfGlycopepTRE_stat = arrangeData(dfGlycopepTRE)

dfProtSiteCTR_stat = arrangeData(dfProtSiteCTR)
dfProtSiteTRE_stat = arrangeData(dfProtSiteTRE)

dfProtSiteCTR_stat = dfProtSiteCTR[dfProtSiteCTR.columns.drop(list(dfProtSiteCTR.filter(regex='Glycopeptide')))].copy()
dfProtSiteTRE_stat = dfProtSiteTRE[dfProtSiteTRE.columns.drop(list(dfProtSiteTRE.filter(regex='Glycopeptide')))].copy()

dfGlycopepCTR_stat["group"] = "CONTROL"
dfGlycopepTRE_stat["group"] = "TREATMENT"
dfProtSiteCTR_stat["group"] = "CONTROL"
dfProtSiteTRE_stat["group"] = "TREATMENT"

mergedGlyco_stat = pd.merge(dfGlycopepCTR_stat, dfGlycopepTRE_stat, on = ["ProteinName", "Position", "Glycopeptide","group"], how = "outer")
mergedGlycoPlot_stat = dict(list(mergedGlyco_stat.groupby(["ProteinName", "Position"])))

mergedSite_stat = pd.merge(dfProtSiteCTR_stat, dfProtSiteTRE_stat, on = ["ProteinName", "Position","group"], how = "outer")
mergedSitePlot_stat = dict(list(mergedSite_stat.groupby(["ProteinName"])))

def getInfo(mergedData, source):
    """
    This function will get the IDs, fold-change
    and pvalues for proteins/glycopeptides found 
    only in control/treatment or both.
    Parameters:
    mergedData = the merged dataset used for plots.
    source = file type, can be "glyco" or "site"
    """
    statDf = pd.DataFrame()
    ctUniques = pd.DataFrame()
    trUniques = pd.DataFrame()
    for i, (k, v) in enumerate(mergedData.items()):
        if source == "glyco":
            meltdf = pd.melt(v, id_vars = ["group", "ProteinName", "Glycopeptide", "Position"])
            filterList = meltdf['Glycopeptide'].unique().tolist()
        elif source == "site":
            meltdf = pd.melt(v, id_vars = ["group", "ProteinName", "Position"])
            filterList = meltdf['Position'].unique().tolist()
        
        for a in filterList:
            if source == "glyco":
                control = meltdf[meltdf['Glycopeptide'] == a].where(meltdf.group == "CONTROL").dropna()['value'].astype(float)
                treatment = meltdf[meltdf['Glycopeptide'] == a].where(meltdf.group == "TREATMENT").dropna()['value'].astype(float)
                site = k[1]
                prot = k[0]
                protID = str(k[0]).split("|")[1]
                gene = str(k[0]).split("|")[2].split("_")[0]
                glycopep = a
                
            elif source == "site":
                control = meltdf[meltdf['Position'] == a].where(meltdf.group == "CONTROL").dropna()['value'].astype(float)
                treatment = meltdf[meltdf['Position'] == a].where(meltdf.group == "TREATMENT").dropna()['value'].astype(float)
                site = a
                prot = str(k)
                protID = str(k).split("|")[1]
                gene = str(k).split("|")[2].split("_")[0]
                
                
            if control.empty and not treatment.empty:
                if source == "glyco":
                    tempTr = pd.DataFrame({'Site' : [site], 'Protein': [prot], 'Protein-ID': [protID], 
                                           'Gene': [gene], 'Glycopeptide': [glycopep]})
                else:
                    tempTr = pd.DataFrame({'Site' : [site], 'Protein': [prot], 'Protein-ID': [protID], 
                                           'Gene': [gene]})
            
                trUniques = pd.concat([trUniques,tempTr])

            elif treatment.empty and not control.empty:
                if source == "glyco":
                    tempCt = pd.DataFrame({'Site' : [site], 'Protein': [prot], 'Protein-ID': [protID], 
                                           'Gene': [gene], 'Glycopeptide': [glycopep]})
                else:
                    tempCt = pd.DataFrame({'Site' : [site], 'Protein': [prot], 'Protein-ID': [protID], 
                                           'Gene': [gene]})

                ctUniques = pd.concat([ctUniques,tempCt])

            elif not control.empty and not treatment.empty:

                fc = np.log2(treatment.mean()/control.mean())
                pval = scipy.stats.ttest_ind(np.log2(control),np.log2(treatment))[1]
                
                if source == "glyco":
                    temp = {'Site' : [site], 'Protein': [prot], 'Protein-ID': [protID], 
                            'Gene': [gene],'Glycopeptide': [glycopep], 'pval': [pval], 'log2-FC': [fc]}
                elif source == "site":
                    temp = {'Site' : [site], 'Protein': [prot], 'Protein-ID': [protID], 
                            'Gene': [gene], 'pval': [pval], 'log2-FC': [fc]}
                    
                tempDf = pd.DataFrame(data = temp)
                statDf = pd.concat([statDf,tempDf])
        
                
    return trUniques, ctUniques, statDf

trUniquesGly, ctUniquesGly, statDfGly = getInfo(mergedGlycoPlot_stat, "glyco")
trUniquesSite, ctUniquesSite, statDfSite = getInfo(mergedSitePlot_stat, "site")

trUniquesGly.to_csv(os.path.join(path,"UniquesGlycopeptidesTreatment.csv"), index = False)
ctUniquesGly.to_csv(os.path.join(path,"UniquesGlycopeptidesControl.csv"), index = False)
statDfGly.to_csv(os.path.join(path,"GlycopepStatistics.csv"), index = False)
trUniquesSite.to_csv(os.path.join(path,"UniquesSitesTreatment.csv"), index = False)
ctUniquesSite.to_csv(os.path.join(path,"UniquesSitesControl.csv"), index = False)
statDfSite.to_csv(os.path.join(path,"SiteStatistics.csv"), index = False)

#Get regulated (without p-val correction)
regGlycos = statDfGly['Glycopeptide'].loc[statDfGly['pval'] < 0.05].tolist()
regGlycos.append(trUniquesGly['Glycopeptide'].tolist())
regGlycos.append(ctUniquesGly['Glycopeptide'].tolist())
regSites = statDfSite['Site'].loc[statDfSite['pval'] < 0.05].tolist()
regSites.append(trUniquesSite['Site'].tolist())
regSites.append(ctUniquesSite['Site'].tolist())

#Make fold change of uniques to be 6.64 or -6.64 and set p-value to simbolic value of zero
trUniquesGly['pval'] = float(0)
trUniquesGly['log2-FC'] = float(6.64)
ctUniquesGly['pval'] = float(0)
ctUniquesGly['log2-FC'] = float(-6.64)
trUniquesSite['pval'] = float(0)
trUniquesSite['log2-FC'] = float(6.64)
ctUniquesSite['pval'] = float(0)
ctUniquesSite['log2-FC'] = float(-6.64)
temp_df = pd.merge(trUniquesGly, statDfGly[statDfGly['pval'] < 0.05], on = ["Site", "Protein", "Protein-ID", "Gene", "Glycopeptide", "pval", "log2-FC"], how = "outer")
completeRelevantGly = pd.merge(temp_df, ctUniquesGly, on = ["Site", "Protein", "Protein-ID", "Gene", "Glycopeptide",  "pval", "log2-FC"], how = "outer")
temp_df = pd.merge(trUniquesSite, statDfSite[statDfSite['pval'] < 0.05], on = ["Site", "Protein", "Protein-ID", "Gene",  "pval", "log2-FC"], how = "outer")
completeRelevantSite = pd.merge(temp_df, ctUniquesSite, on = ["Site", "Protein", "Protein-ID", "Gene",  "pval", "log2-FC"], how = "outer")

fuc = 0
fuc_log = []
highm = 0
highm_log = []
pauci = 0
pauci_log = []
sial = 0
sial_log = []
hycomp = 0
hycomp_log = []
for i, (k, v) in enumerate(dict(list(completeRelevantGly.groupby(["Glycopeptide"]))).items()):
    glyp = k.split(".")[3]
    manno = glyp.split('Hex(')
    hexnac = glyp.split('HexNAc(')[1].split(')')[0]
    if "Fuc" in glyp:
        fuc = fuc + 1
        fuc_log.append(v.iloc[0]['log2-FC'])
    elif "Fuc" not in glyp and "NeuAc" in glyp or "Fuc" not in glyp and "NeuGc" in glyp:
        sial = sial + 1
        sial_log.append(v.iloc[0]['log2-FC'])
    elif "Fuc" not in glyp and "NeuAc" not in glyp and "NeuGc" not in glyp and len(manno) > 1:
        hexnac = int(glyp.split('HexNAc(')[1].split(')')[0])
        manno_count = int(manno[1].split(')')[0]) 
        #print(manno_count)
        if manno_count > 3 and hexnac == 2:
            #print(str(manno_count) + "  " + str(hexnac))
            highm = highm + 1
            highm_log.append(v.iloc[0]['log2-FC'])
        if manno_count <= 3 and hexnac == 2:
            pauci = pauci + 1
            pauci_log.append(v.iloc[0]['log2-FC'])
    else: 
        hycomp = hycomp + 1
        hycomp_log.append(v.iloc[0]['log2-FC'])


upreg_fuc = pd.Series(fuc_log).gt(0).sum()
downreg_fuc = pd.Series(fuc_log).lt(0).sum()
upreg_highm = pd.Series(highm_log).gt(0).sum()
downreg_highm = pd.Series(highm_log).lt(0).sum()
upreg_sial = pd.Series(sial_log).gt(0).sum()
downreg_sial = pd.Series(sial_log).lt(0).sum()
upreg_pauci = pd.Series(pauci_log).gt(0).sum()
downreg_pauci = pd.Series(pauci_log).lt(0).sum()
upreg_hycomp = pd.Series(hycomp_log).gt(0).sum()
downreg_hycomp = pd.Series(hycomp_log).lt(0).sum()

#Pie chart
glycos = ['Fuc', 'High. Mann', 'Pauci','Sial', "Hybrid/Complex"]
data = [fuc, highm, pauci, sial, hycomp]
devi = (0.05, 0.15, 0.22, 0.1, 0.15)
 
colors = ( "darkgreen", "cadetblue", "violet",
          "salmon", "orangered")

wp = { 'linewidth' : 1, 'edgecolor' : "black" }
def func(pct, allvalues):
    absolute = int(pct / 100.*np.sum(allvalues))
    return "{:.1f}%\n({:d} )".format(pct, absolute)
 
fig, ax = plt.subplots(figsize =(10, 7))
wedges, texts, autotexts = ax.pie(data,
                                  autopct = lambda pct: func(pct, data),
                                  explode = devi,
                                  labels = glycos,
                                  shadow = False,
                                  colors = colors,
                                  startangle = 90,
                                  wedgeprops = wp,
                                  textprops = dict(color ="black"))
 
ax.legend(wedges, glycos,
          title ="Gly",
          loc =2,
          bbox_to_anchor =(1, 0, 0.5, 1),  borderaxespad=0.)
 
plt.setp(autotexts, size = 8, weight ="bold")
ax.set_title("Glycans distribution")
plt.savefig(os.path.join(pathGly, "glycansPieChart"), dpi = 300)
plt.close('all')

#Barplot
ups = [upreg_fuc, upreg_highm, upreg_pauci, upreg_sial, upreg_hycomp]
downs = [downreg_fuc, downreg_highm, downreg_pauci, downreg_sial, downreg_hycomp]
width = 0.35

df_bar = pd.DataFrame(data={'UP': ups, 'DOWN': downs}, index=glycos)

ax = df_bar.plot(kind='barh', color = ["firebrick", "navy"], ylabel='Glycans', title='Count by Glycans (Treatment/Control)')
ax.set(xlabel='Count')

total = sum(df_bar.sum())
for c, col in zip(ax.containers, df_bar.columns):
    ax.bar_label(c, fontsize = 6,label_type='edge', labels=[f'{val}\n{val / total * 100.0:.1f} %' for val in df_bar[col]])

ax.legend(title='Group', bbox_to_anchor=(1, 1.02), loc='upper left')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()

plt.savefig(os.path.join(pathGly, "glycansBarPlot"), dpi = 600, bbox_inches="tight")
plt.close('all')

#Plotting the Glycopeptides percentage per site per protein among all samples
def statCalc(old, source):
	"""
	This function will calculate the standard deviation
	and mean values among the replicates of each group.
	Parameter:
	old = the original dataset.
	type = if it is glyco or site.
	source = file type, can be "glyco" or "site"
	"""
	new = old.copy()
	#new = new[new.columns.drop(list(new.filter(regex='TotalIonCurrent')))]
	new.dropna(axis=0, thresh=(int(len(new.filter(regex=r'icPerc').columns)*float(threshNaN))), subset = list(new.filter(regex=r'icPercentage').columns), inplace = True)
	new = new.fillna(0)
	if source == "glyco":
		new["std"] = new.filter(regex=r'icPercentage').std(axis = 1)
		new["mean"] = new.filter(regex=r'icPercentage').mean(axis = 1)
	elif source == "site":
		new["std"] = new.filter(regex=r'icPercentage').std(axis = 1)
		new["mean"] = new.filter(regex=r'icPercentage').mean(axis = 1)

	return new

dfGlycopepCTRstat = statCalc(dfGlycopepCTR, "glyco")
dfGlycopepTREstat = statCalc(dfGlycopepTRE, "glyco")
dfProtSiteCTRstat = statCalc(dfProtSiteCTR, "site")
dfProtSiteTREstat = statCalc(dfProtSiteTRE, "site")
dfProtSiteCTRstat.drop(list(dfProtSiteCTRstat.filter(regex='Glycope')), axis = 1, inplace = True)
dfProtSiteCTRstat.drop(list(dfProtSiteCTRstat.filter(regex='fullIon')), axis = 1, inplace = True)
dfProtSiteCTRstat.drop(list(dfProtSiteCTRstat.filter(regex='TotalIon')), axis = 1, inplace = True)
dfProtSiteTREstat.drop(list(dfProtSiteTREstat.filter(regex='Glycope')), axis = 1, inplace = True)
dfProtSiteTREstat.drop(list(dfProtSiteTREstat.filter(regex='fullIon')), axis = 1, inplace = True)
dfProtSiteTREstat.drop(list(dfProtSiteTREstat.filter(regex='TotalIon')), axis = 1, inplace = True)

#Merge groups for plotting 
dfGlycopepCTRstat["group"] = "CONTROL"
dfGlycopepTREstat["group"] = "TREATMENT"
dfProtSiteCTRstat["group"] = "CONTROL"
dfProtSiteTREstat["group"] = "TREATMENT"

mergedGlyco = pd.merge(dfGlycopepCTRstat, dfGlycopepTREstat, on = ["ProteinName", "Position", "Glycopeptide", "group"], how = "outer")
mergedGlyco = mergedGlyco.fillna(0)
mergedGlyco = mergedGlyco[mergedGlyco.columns.drop(list(mergedGlyco.filter(regex='TotalIonCurrent')))]
mergedGlyco = mergedGlyco[mergedGlyco.columns.drop(list(mergedGlyco.filter(regex='fullIonCurrent')))]
mergedGlyco.drop(['std_x', 'std_y', 'mean_x', 'mean_y'], axis = 1, inplace = True )
mergedGlycoPlot = dict(list(mergedGlyco.groupby(["ProteinName", "Position"])))
mergedGlycoSave = pd.merge(dfGlycopepCTRstat, dfGlycopepTREstat, on = ["ProteinName", "Position", "Glycopeptide"], how = "outer")
mergedGlycoSave.to_csv(os.path.join(path,"glyControlxTreatment.csv"), index = False)

for i, (k, v) in enumerate(mergedGlycoPlot.items()):
	if plotDecision == "complete":
		df = v
	if plotDecision == "relevant":
		df = v[v['Glycopeptide'].isin(regGlycos)]
	if len(df) > 0:
	    saveName = str(k[0]).split("|")[1] +  "-Position " + str(int(k[1])) + ".tiff"
	    savePath = (os.path.join(pathGly, saveName))
	    plotdf = pd.melt(df, id_vars = ["group", "ProteinName", "Glycopeptide", "Position"])
	    plt.figure(figsize = (5,5))
	    ax = sns.boxplot(x='Glycopeptide', y='value', hue='group',hue_order = ['CONTROL', 'TREATMENT'], data=plotdf)
	    ax.set_xticklabels(ax.get_xticklabels(),rotation = 45, ha = 'right')
	    ax.set(title=k)
	    ax.set_ylabel("Percentage", fontsize = 15)
	    ax.set_xlabel("Glycopeptide", fontsize = 15)
	    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	    plt.savefig(savePath, dpi = 300, bbox_inches="tight")
	    plt.close("all")

mergedSite = pd.merge(dfProtSiteCTRstat, dfProtSiteTREstat, on = ["ProteinName", "Position",  "group"], how = "outer")
mergedSite = mergedSite.fillna(0)
mergedSite.drop(['std_x', 'std_y', 'mean_x', 'mean_y'], axis = 1, inplace = True )
mergedSitePlot = dict(list(mergedSite.groupby(["ProteinName"])))
mergedSiteSave = pd.merge(dfProtSiteCTRstat, dfProtSiteTREstat, on = ["ProteinName", "Position"], how = "outer")
mergedSiteSave.to_csv(os.path.join(path,"siteControlxTreatment.csv"), index = False)

for i, (k, v) in enumerate(mergedSitePlot.items()):
	if plotDecision == "complete":
		df = v
	if plotDecision == "relevant":
		df = v[v['Position'].isin(regSites)]
	if len(df) > 0:
	    saveNameSite = str(k).split("|")[1] + ".tiff"
	    savePathSite = (os.path.join(pathSite, saveNameSite))
	    plotdf = pd.melt(df, id_vars = ["group", "ProteinName", "Position"])
	    plt.figure(figsize = (5,5))
	    ax = sns.boxplot(x='Position', y='value', hue='group',hue_order = ['CONTROL', 'TREATMENT'], data=plotdf)
	    ax.set_xticklabels(ax.get_xticklabels(),rotation = 45, ha = 'right')
	    ax.set(title=k)
	    ax.set_ylabel("Percentage", fontsize = 15)
	    ax.set_xlabel("Site", fontsize = 15)
	    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	    plt.savefig(savePathSite, dpi = 300, bbox_inches="tight")
	    plt.close("all")


fuc = 0
fuc_log = []
highm = 0
highm_log = []
pauci = 0
pauci_log = []
sial = 0
sial_log = []
hycomp = 0
hycomp_log = []
for i, (k, v) in enumerate(dict(list(dfGlycopepCTR.groupby(["Glycopeptide"]))).items()):
    glyp = k.split(".")[3]
    manno = glyp.split('Hex(')
    hexnac = glyp.split('HexNAc(')[1].split(')')[0]
    if "Fuc" in glyp:
        fuc = fuc + 1
        
    elif "Fuc" not in glyp and "NeuAc" in glyp or "Fuc" not in glyp and "NeuGc" in glyp:
        sial = sial + 1
        
    elif "Fuc" not in glyp and "NeuAc" not in glyp and "NeuGc" not in glyp and len(manno) > 1:
        hexnac = int(glyp.split('HexNAc(')[1].split(')')[0])
        
        
        if manno_count > 3 and hexnac == 2:
            
            highm = highm + 1
            
        if manno_count <= 3 and hexnac == 2:
            pauci = pauci + 1
            
    else: 
        hycomp = hycomp + 1
        

#Pie chart
glycos = ['Fuc', 'High. Mann', 'Pauci','Sial', "Hybrid/Complex"]
data = [fuc, highm, pauci, sial, hycomp]
devi = (0.05, 0.15, 0.22, 0.1, 0.15)
 
colors = ( "darkgreen", "cadetblue", "violet",
          "salmon", "orangered")

wp = { 'linewidth' : 1, 'edgecolor' : "black" }
def func(pct, allvalues):
    absolute = int(pct / 100.*np.sum(allvalues))
    return "{:.1f}%\n({:d} )".format(pct, absolute)
 
fig, ax = plt.subplots(figsize =(10, 7))
wedges, texts, autotexts = ax.pie(data,
                                  autopct = lambda pct: func(pct, data),
                                  explode = devi,
                                  labels = glycos,
                                  shadow = False,
                                  colors = colors,
                                  startangle = 90,
                                  wedgeprops = wp,
                                  textprops = dict(color ="black"))
 
ax.legend(wedges, glycos,
          title ="Gly",
          loc =2,
          bbox_to_anchor =(1, 0, 0.5, 1),  borderaxespad=0.)
 
plt.setp(autotexts, size = 8, weight ="bold")
ax.set_title("Glycans distribution (Total Control)")
plt.savefig(os.path.join(pathGly, "glycansPieChartTotalControl"), dpi = 300)
plt.close('all')


fuc = 0
fuc_log = []
highm = 0
highm_log = []
pauci = 0
pauci_log = []
sial = 0
sial_log = []
hycomp = 0
hycomp_log = []
for i, (k, v) in enumerate(dict(list(dfGlycopepTRE.groupby(["Glycopeptide"]))).items()):
    glyp = k.split(".")[3]
    manno = glyp.split('Hex(')
    hexnac = glyp.split('HexNAc(')[1].split(')')[0]
    if "Fuc" in glyp:
        fuc = fuc + 1
        
    elif "Fuc" not in glyp and "NeuAc" in glyp or "Fuc" not in glyp and "NeuGc" in glyp:
        sial = sial + 1
        
    elif "Fuc" not in glyp and "NeuAc" not in glyp and "NeuGc" not in glyp and len(manno) > 1:
        hexnac = int(glyp.split('HexNAc(')[1].split(')')[0])
        
        
        if manno_count > 3 and hexnac == 2:
            
            highm = highm + 1
            
        if manno_count <= 3 and hexnac == 2:
            pauci = pauci + 1
            
    else: 
        hycomp = hycomp + 1
        

#Pie chart
glycos = ['Fuc', 'High. Mann', 'Pauci','Sial', "Hybrid/Complex"]
data = [fuc, highm, pauci, sial, hycomp]
devi = (0.05, 0.15, 0.22, 0.1, 0.15)
 
colors = ( "darkgreen", "cadetblue", "violet",
          "salmon", "orangered")

wp = { 'linewidth' : 1, 'edgecolor' : "black" }
def func(pct, allvalues):
    absolute = int(pct / 100.*np.sum(allvalues))
    return "{:.1f}%\n({:d} )".format(pct, absolute)
 
fig, ax = plt.subplots(figsize =(10, 7))
wedges, texts, autotexts = ax.pie(data,
                                  autopct = lambda pct: func(pct, data),
                                  explode = devi,
                                  labels = glycos,
                                  shadow = False,
                                  colors = colors,
                                  startangle = 90,
                                  wedgeprops = wp,
                                  textprops = dict(color ="black"))
 
ax.legend(wedges, glycos,
          title ="Gly",
          loc =2,
          bbox_to_anchor =(1, 0, 0.5, 1),  borderaxespad=0.)
 
plt.setp(autotexts, size = 8, weight ="bold")
ax.set_title("Glycans distribution (Total Treatment)")
plt.savefig(os.path.join(pathGly, "glycansPieChartTotalTreatment"), dpi = 300)
plt.close('all')
