# usage: python3 slimModifier.py slimFle csvFile
# python3 slimModifier.py neutral.slim simulation_freqs.csv 

import sys
import pandas as pd
import os


os.makedirs("slim_scripts")

#Read files from command line
slimFileName = sys.argv[1]
csvFile = sys.argv[2]
#Open CSV extract columns
df = pd.read_csv(csvFile)
number = df.iloc[:,0]
name = df.iloc[:,1]

#Read .slim file which is to be modified
slimFile = open(slimFileName,"r+").read().split("\n")
pos1 = 0
pos2 = 0
for index, data in enumerate(slimFile):
    if "FIRST TARGET" in data:
        pos1 = index + 1
    if "SECOND TARGET" in data:
        pos2 = index + 1

#Create a slim file for each entry in the .csv file, with 'name + .slim' 
for i in range(len(number)):
    first = slimFile[pos1].index(",")
    slimFile[pos1] = slimFile[pos1][:27] + str(number[i]) + ");"

    second = slimFile[pos2].index(",")
    slimFile[pos2] = slimFile[pos2][:27] + str(number[i]) + ");"
    
    #Outputs will be in the result directory
    with open("slim_scripts/" + str(name[i]) + ".slim", "w") as file:
        for j in slimFile:
            file.write(str(j) + "\n")
