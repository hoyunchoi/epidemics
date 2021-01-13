import numpy as np
import pandas as pd
import glob

rootDirectory = "../data/epidemics/"

def readCSV(t_fileName):
    data = pd.read_csv(t_fileName, sep=',', header=None)
    data = data.values.transpose()
    if (len(data) == 1):
        return data[0]
    else:
        return tuple([row for row in data])

#* Naming Conventions
def SIRFileName(t_networkSize, t_meanDegree, t_SI_II, t_I_R, t_randomEngineSeed=None):
    fileName = "N{:.1e},M{:d},SIII{:.2f},IR{:.2f}".format(t_networkSize, t_meanDegree, t_SI_II, t_I_R)
    if (t_randomEngineSeed == None):
        return fileName + ".txt"
    else:
        return fileName + "-{:d}.txt".format(t_randomEngineSeed)

def SEIRFileName(t_networkSize, t_meanDegree, t_SI_EI, t_E_I, t_I_R, t_randomEngineSeed=None):
    fileName = "N{:.1e},M{:d},SIEI{:.2f},EI{:.2f},IR{:.2f}".format(t_networkSize, t_meanDegree, t_SI_EI, t_E_I, t_I_R)
    if (t_randomEngineSeed == None):
        return fileName + ".txt"
    else:
        return fileName + "-{:d}.txt".format(t_randomEngineSeed)

def GSEIRFileName(t_networkSize, t_meanDegree, t_SE_EE, t_SI_EI, t_E_I, t_I_R, t_randomEngineSeed=None):
    fileName = "N{:.1e},M{:d},SEEE{:.2f},SIEI{:.2f},EI{:.2f},IR{:.2f}".format(t_networkSize, t_meanDegree, t_SE_EE, t_SI_EI, t_E_I, t_I_R)
    if (t_randomEngineSeed == None):
        return fileName + ".txt"
    else:
        return fileName + "-{:d}.txt".format(t_randomEngineSeed)


#* Read Data
def readSIR(t_networkType, t_spreadType, t_networkSize, t_meanDegree, t_SI_II, t_I_R, t_randomEngineSeed=None):
    FullFileName = rootDirectory + "SIR/" + t_networkType + "/" + t_spreadType + "/" + SIRFileName(t_networkSize, t_meanDegree, t_SI_II, t_I_R, t_randomEngineSeed)
    return readCSV(FullFileName)

def readSEIR(t_networkType, t_spreadType, t_networkSize, t_meanDegree, t_SI_EI, t_E_I, t_I_R, t_randomEngineSeed=None):
    FullFileName = rootDirectory + "SEIR/" + t_networkType + "/" + t_spreadType + "/" + SEIRFileName(t_networkSize, t_meanDegree, t_SI_EI, t_E_I, t_I_R, t_randomEngineSeed)
    return readCSV(FullFileName)

def readGSEIR(t_networkType, t_spreadType, t_networkSize, t_meanDegree, t_SE_EE, t_SI_EI, t_E_I, t_I_R, t_randomEngineSeed=None):
    FullFileName = rootDirectory + "GSEIR/" + t_networkType + "/" + t_spreadType + "/" + GSEIRFileName(t_networkSize, t_meanDegree, t_SE_EE, t_SI_EI, t_E_I, t_I_R, t_randomEngineSeed)
    return readCSV(FullFileName)


if __name__ == "__main__":
    print("This is a moduel readData.py")