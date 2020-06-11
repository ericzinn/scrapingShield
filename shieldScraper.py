# For Data Management (csv files, dataframes, directory handling)
import pandas as pd
import numpy as np
import csv
import argparse
from pathlib import Path

# For interfacing with the website
import requests
from bs4 import BeautifulSoup
import urllib.request

# For making pretty pictures
import seaborn as sns
import matplotlib.pyplot as plt

# Begin Function Definitions
def combineTablesByCondition(listOfDataframes, condition):
    """Given a list of dataframes and a condition, combine those
       dataframes into a single, aggregated dataframe"""
    
    # Initialize a new, empty dataframe
    aggregateDf = pd.DataFrame()

    for tmpDf in listOfDataframes:
        aggregateDf[tmpDf.columns[0]] = tmpDf[tmpDf[tmpDf.columns[0]] == condition].squeeze()

    return aggregateDf.T.drop(aggregateDf.T.columns[0], axis = 1).fillna(0)

def heatMap(dataFrame, conditionName):
    """Render a heatmap from a dataframe, save it with the provided
       condition name"""

    # Set figure size
    # Going to allocate 0.25" per line to ensure readability
    plt.figure(figsize = (5, 0.25*len(dataFrame)))
    ax = sns.heatmap(dataFrame, cmap = "bwr")
    plt.savefig("heatmaps/"+conditionName+".png")

def clusterMap(dataFrame, conditionName, logScale = False, saveFig = True):
    """Render clustermap from a dataframe, save it (by default) with the
       provided condition name.  Includes options to pre-process dataframe
       into log-transformed data"""
    if logScale == True:
        dataFrame = np.log10(dataFrame + 1e-1)
    # Set figure size
    # Going to allocate 0.25" per line to ensure readability
    g = sns.clustermap(dataFrame, method = 'ward', col_cluster = False, cmap = "bwr", z_score = None,
                       cbar_pos=(1.0, .2, .03, .4), dendrogram_ratio = (0.2,0.01), figsize = (5, 0.33*len(dataFrame)))
    #g.cax.set_visible(False)
    #g.ax_row_dendrogram.set_visible(False)
    if saveFig == True:
        if logScale == False:
            plt.savefig("heatmaps/"+conditionName+"(Clustered).png", bbox_inches='tight',pad_inches = 0)
        else:
            plt.savefig("heatmaps/"+conditionName+"(Clustered, log10).png", bbox_inches='tight',pad_inches = 0)
# End Function Definitions

# Set up argument parser
parser = argparse.ArgumentParser()

# Required Argument(s)
requiredNamed = parser.add_argument_group("Required Named Arguments")
parser.add_argument("-i", "--input-file", default = "",
                    help = "Path to gene list (csv format)")

# Store argument variables                
args = vars(parser.parse_args())

# Allocate arguments into local vars
pathToFile = args['input_file']

# Read genes from the provided list (Utf-8 encoded CSV file)
geneFile = open(pathToFile, 'r', encoding='utf-8-sig')
reader = csv.reader(geneFile)
geneList = [row[0] for row in reader]

# Make new directories to hold our output data
Path("charts").mkdir(parents=True, exist_ok=True)
Path("heatmaps").mkdir(parents=True, exist_ok=True)
Path("csvFiles/individualGenes/").mkdir(parents=True, exist_ok=True)

# Define the 'base' URL which we're going to use to request pages corresponding
# to each of the genes we're interested in
baseUrl = "https://shield.hms.harvard.edu/viewgene.html?gene="

# Scrape the charts
print("Scraping charts...")
for gene in geneList:
    gene = gene.capitalize()
    page = requests.get(baseUrl+gene)
    soup = BeautifulSoup(page.content, 'html.parser')
    chart = soup.find(id = "FACS_chart")
    # Un-comment the following line if you want to get realtime progress
    #print("Downloading chart for gene: "+gene)
    #print(chart['href'])
    try:
        urllib.request.urlretrieve(chart['href'], "charts/"+gene+".png")
    except:
        print("Error: No RNASeq data found for "+gene)

# Scrape the raw data
print("Scraping FACS Sorted RNASeq Data...")

# Initialize an empty list to store all of our tables
dfList = []

for gene in geneList:
    gene = gene.capitalize()
    # Uncomment the following line to get real-time progress updates
    #print("Downloading table for gene: "+gene)
    page = requests.get(baseUrl+gene)
    soup = BeautifulSoup(str(page.content), 'lxml')
    try:
        table = soup.find(id = "FACS_data_table")
    except:
        table = None
    if table:
        df = pd.read_html(str(table))[0].rename(columns = {"Unnamed: 0":gene})
        df.to_csv("csvFiles/individualGenes/"+gene+".csv")
        dfList.append(df)
    else:
        print("Error: No RNASeq data found for "+gene)

# Begin massaging Data by getting a list of conditions from the first
# table that we scraped
listOfConditions = dfList[0][dfList[0].columns[0]].unique()

print("Massaging data...")
# Combine Tables to create new dataframes for each experimental condition
UtricleGFPPosDf = combineTablesByCondition(dfList, listOfConditions[0])
CochleaGFPPosDf = combineTablesByCondition(dfList, listOfConditions[1])
UtricleGFPNegDf = combineTablesByCondition(dfList, listOfConditions[2])
CochleaGFPNegDf = combineTablesByCondition(dfList, listOfConditions[3])

print("Outputting aggregated data to csv files...")
UtricleGFPPosDf.to_csv("csvFiles/"+listOfConditions[0]+".csv")
CochleaGFPPosDf.to_csv("csvFiles/"+listOfConditions[1]+".csv")
UtricleGFPNegDf.to_csv("csvFiles/"+listOfConditions[2]+".csv")
CochleaGFPNegDf.to_csv("csvFiles/"+listOfConditions[3]+".csv")
# Render basic heatmaps
print("Rendering heatmaps...")
heatMap(UtricleGFPPosDf, listOfConditions[0])
heatMap(CochleaGFPPosDf, listOfConditions[1])
heatMap(UtricleGFPNegDf, listOfConditions[2])
heatMap(CochleaGFPNegDf, listOfConditions[3])

# Render clustered heatmaps
print("Rendering clustered heatmaps...")
clusterMap(UtricleGFPPosDf, listOfConditions[0])
clusterMap(CochleaGFPPosDf, listOfConditions[1])
clusterMap(UtricleGFPNegDf, listOfConditions[2])
clusterMap(CochleaGFPNegDf, listOfConditions[3])

# Render log-transformed clustered heatmaps
print("Rendering log-transformed clustered heatmaps...")

clusterMap(UtricleGFPPosDf, listOfConditions[0], logScale=True)
clusterMap(CochleaGFPPosDf, listOfConditions[1], logScale=True)
clusterMap(UtricleGFPNegDf, listOfConditions[2], logScale=True)
clusterMap(CochleaGFPNegDf, listOfConditions[3], logScale=True)

print("Finished!")