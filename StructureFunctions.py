import FileFunctions as FF, conf
import os, shutil
import subprocess as Sp
from collections import defaultdict
import numpy as np, math
# Base Pairs from dot bracket Secondary structure
def ListBasePairsFromStruct(Struct):  # return dic={structure:[liste de pairs de base ],....}
    lista = []
    stack = []
    for i in range(len(Struct)):  # sequence length
        if Struct[i] == '(':  # Opening base pair
            stack.append(i)
        elif Struct[i] == ')':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))
    return lista

# Parse an RNAsubopt file to extract Base pairs
def GetBasePairsFromStructFile(faPath):   #return dic={structure:[liste de pairs de base ],....}
    #print faPath
    DicStruct={}
    lines = FF.Parsefile(faPath)
    #print lines
    SeqLen=len(lines[1])-1
    #print SeqLen,"seq length"
    for j in range(len(lines)):
        DicStruct[j]=ListBasePairsFromStruct(lines[j].strip().split(' ')[0])
    return len(lines),DicStruct

#Calculate distance between each 2 structures and generate a matrix 2 by 2 , stored in the file SVMLFile
def DistanceTwoBPlist(Struct1,Struct2):
    return len(set(Struct1).symmetric_difference(set(Struct2)) )
    #return len(set(Struct1).intersection(set(Struct2)) )

def dd():
    return 0
def aa():
    return defaultdict(dd)
def DistanceStruct(StructFile, SVMlFile, numberofsruct, MFESnbrstruct, constrainte):
    Redondantestructure = defaultdict(aa)
    MatDist = defaultdict(aa)
    Redondantestructure1 = []
    DicStruct = {}
    Dicnumberofsruct = {}

    for i in range(len(constrainte) - 1):
        Dicnumberofsruct[constrainte[i]] = numberofsruct
    Dicnumberofsruct[constrainte[len(constrainte) - 1]] = MFESnbrstruct

    nb, DicStruct = GetBasePairsFromStructFile(StructFile)

    for i in range(0, nb):
        for j in range(i + 1, nb):
            MatDist[i][j] = DistanceTwoBPlist(DicStruct[i], DicStruct[
                j])
            if MatDist[i][j] == 0:
                if j not in Redondantestructure1:
                    if j > numberofsruct * (len(constrainte) - 1):
                        Dicnumberofsruct[constrainte[len(constrainte) - 1]] -= 1
                    else:
                        Dicnumberofsruct[constrainte[int(j / numberofsruct)]] -= 1
                    Redondantestructure1.append(j)

            MatDist[j][i] = MatDist[i][j]


    for elem in Redondantestructure1:
        if elem < numberofsruct * (len(constrainte) - 1):
            ConditionNumber = int((elem) / numberofsruct)
        else:
            ConditionNumber = len(constrainte) - 1
        StructureNumber = elem - ConditionNumber * numberofsruct
        Redondantestructure[constrainte[ConditionNumber]][StructureNumber] = 1

    # strore the distance matrix in the file SVMLFile
    o = open(os.path.join("output",SVMlFile), "w")
    for i in range(len(MatDist)):
        o.write("%i\t" % (i + 1))
        for j in range(len(MatDist)):
            if (i != j):
                o.write("%i:%.4f\t" % (j + 1, MatDist[i][j]))
        o.write("\n")
    o.close()

    if Redondantestructure!=0:
        print "Warning! redundant structures"
    FF.PickleVariable(MatDist, os.path.join(conf.PickledData,"dissmatrix.pkl"))
    FF.PickleVariable(Redondantestructure1, os.path.join(conf.PickledData,"Redondantestructures.pkl"))
    FF.PickleVariable(Redondantestructure, os.path.join(conf.PickledData, "Redondantestructures_Id.pkl"))
    FF.PickleVariable(Dicnumberofsruct,os.path.join(conf.PickledData,"Dicnumberofsruct.pkl"))
    return 0


def FromStructFiletoRNAEvalInput(StructFile, InputRNAeval, rna):
    lines = FF.Parsefile(StructFile)
    o = open(InputRNAeval,
             "w")  # geneate intermediate file with sequence+strcuture , seq+strcture .... as the input format  to use RNAeval
    # print "sdfqspkojr",len(lines)
    for i in range(1, len(lines)):
        o.write("%s%s\t" % (rna, lines[i]))
    o.close()

#StructFile contains the RNA sequence in the first line and list of correponding structures by line
def ENERGY_VALUES_STRUCTURES(StructFile,rna):
        Energy=[]
        #generate the rnaeval input file
        FromStructFiletoRNAEvalInput(StructFile,"InputRNAeval",rna)
     # launch the RNaeval command
        os.system('RNAeval <' + "InputRNAeval" + '>' + "energyvalues")
        #Sp.call("RNAeval " , stdin="InputRNAeval", stdout="energyvalues", shell=True)
        # Parse the RNAevaloutput to extract energy values
        lines=FF.Parsefile("energyvalues")
        for i in xrange(1,len(lines),2):
            # i is the stucture number and 'lines[i].split(" ")[1][1:-2]' is  the  corresponding  energy value
        #print 'holla',(lines[i].split(" ")[1][1:-2])
            Energy.append(lines[i].split(" ")[1][1:-2])
        return Energy

##Boltzmman energy   according to the formula B=exp^\frac{-e}{RT}
def BoltzmannEnergy(Energy):
    #print Energy,"fdfd"
    T=37+273.15
    R=0.0019872370936902486
    return np.exp(-float(Energy)/float(100.*R*T))

def Boltzmann_Calc(constraintes, StructfileRepos, numberofsruct, MFESnbrstruct, rna, Redondantestructure):
    Energy = defaultdict(aa)
    Boltzman = defaultdict(aa)
    ConditionalBoltzman = defaultdict(aa)
    ZBolzman = defaultdict(aa)

    for Condition in constraintes:
        FileStructure = StructfileRepos + '/' + Condition
        #print FileStructure,"ffftft"
        Energy[Condition] = ENERGY_VALUES_STRUCTURES(FileStructure,rna)  # list of energy values for the structures present in the Condition
    #print Energy,"llllllllllllllllllllll","done"
    for Condition in constraintes:
        #print MFESnbrstruct
        if Condition == "MFES":

            Boltzman[Condition] = [BoltzmannEnergy(Energy[Condition][i]) for i in range(2)]
            #print "heehrh",Boltzman[Condition]
            #for i in range(MFESnbrstruct):
            #   print i, Boltzman[Condition][i]
            #print Boltzman[Condition],"mfe"
        else:
            listawithoutRedonddnace = []
            for i in range(numberofsruct):
                Boltzman[Condition][i] = BoltzmannEnergy(Energy[Condition][i])
                if Redondantestructure[Condition][i]==0:  # if the structure is not redundant
                    listawithoutRedonddnace.append(BoltzmannEnergy(Energy[Condition][i]))

        #print Boltzman, "eeeeeeeeeeeeeeeeeeeeeeeeee"
        ZBolzman[Condition] = sum(listawithoutRedonddnace)  # Partition function

    #FF.PickleVariable(Boltzman, os.path.join(conf.PickledData, "Boltzman.pkl"))
    listall = []
    for Condition in constraintes[:-1]:  # to not count MFES
        lista = []
        for i in range(numberofsruct):
            if Redondantestructure[Condition][i]==0:
                lista.append(BoltzmannEnergy(Energy[Condition][i]) / ZBolzman[Condition])
            else:
                lista.append(0)  # to solve the problem of the number of structure variation
        listall += lista
        ConditionalBoltzman[Condition] = lista

    FF.PickleVariable(ConditionalBoltzman, os.path.join(conf.PickledData, "ConditionalBoltzman.pkl"))
    FF.PickleVariable(ZBolzman, os.path.join(conf.PickledData, "ZBolzman.pkl"))

    return Boltzman
from collections import Counter
def OccurenceBasePairs(ListPairs, Threshold):
	return [(elem[0],elem[1],Counter(ListPairs)[elem]) for elem in Counter(ListPairs)  if Counter(ListPairs)[elem]>=Threshold]


def ListBasePairsFromStruct(Struct):  # return dic={structure:[liste de pairs de base ],....}
    lista = []
    stack = []
    for i in range(len(Struct)):  # sequence length
        if Struct[i] == '(':  # Opening base pair
            stack.append(i)
        elif Struct[i] == ')':  # Closing base pair
            k = stack.pop()
            lista.append((k, i))
    return lista

def fromPairsToStruct(rna, Pairs):
    structure=["." for i in range(len(rna)-1)]
    for (i,j) in Pairs:
        structure[i]='('
        structure[j]=')'
    return "".join(structure)
def EnergyValuesFromStructure(StructFile, rna):
    Energy = []
    # generate the rnaeval input file
    FromStructFiletoRNAEvalInput(StructFile, "InputRNAeval", rna)
    # launch the RNaeval command
    os.system('RNAeval <' + "InputRNAeval" + '>' + "energyvalues")
    # Parse the RNAevaloutput to extract energy values
    lines = FF.Parsefile("energyvalues")
    for i in xrange(1, len(lines), 2):
        Energy.append(lines[i].split(" ")[1][1:-2])
    return Energy

import scipy
# Calclate eucledian distance between 2 probability dot plot matrix
def Eucledian_distance(B,lenseq):
	dist=scipy.zeros([len(B),len(B)])
	for i in range(len(B)):
		for j in range(i+1,len(B)):
				for x in range(lenseq):
					for y in range(x+1,lenseq):
						#Just upper right half of the dotplot should be considered
						dist[i][j]+=math.pow(B[i][x][y]-B[j][x][y],2)
                                dist[j][i]=dist[i][j]
				#print "i",i,"j",j,math.sqrt(dist[i][j])
	return dist

# generate dot plot base pairs
def DotplotRnaFold(dir,PathConstrainteFile,PathConstrainteFileShape):
    FF.CreateFold(dir)
    for filename in FF.GetListFile(PathConstrainteFile,'.fa'):
        name=os.path.splitext(filename)[0]
        Input= PathConstrainteFile+"/"+filename
        output=dir+'/'+name
        #print "command is", 'RNAfold -p -d2  -C --enforceConstraint <'+Input+  '>'+ output
        os.system('RNAfold --noLP -p -d2  -C --enforceConstraint <'+Input+  '>'+ output)
        #Sp.call("RNAfold -p -d2  -C ", stdin=Input, stdout=output, shell=True)
        # redirect files to the specific folder dir
        shutil.move(name+'_dp.ps', dir+"/"+name+'_dp.ps')
        shutil.move(name+'_ss.ps', dir+"/"+name+'_ss.ps')
    for filename in FF.GetListFile(PathConstrainteFileShape,'.fa'):
        name=os.path.splitext(filename)[0]
        Input= PathConstrainteFileShape+"/"+filename
        ShapeFile=PathConstrainteFileShape+"/"+name+'Shape.txt'
        output=dir+'/'+name
        os.system('RNAfold -p -d2  -C ' + '<' +Input+  ' --shape '+ ShapeFile + '>'+ output)
        shutil.move(name+'_dp.ps', dir+"/"+name+'_dp.ps')
        shutil.move(name+'_ss.ps', dir+"/"+name+'_ss.ps')


#def extract matrix values
def Writeproba( dir,Matrixproba,constraintes,rna):
    FF.CreateFold(Matrixproba)
    for file in constraintes:
        PSPath =dir+"/"+file+"_dp.ps"
        bpm = loadDotPlotPS(PSPath,"RNAfold")
        dp = DotPlot(rna,bpm)
        with open(Matrixproba +file.split('.')[0]+".proba","w") as o:
            for (i,j)in bpm.keys():
                o.write("%i\t%i\t%.6f\n"%(i,j,bpm[(i,j)]))
        o.close()

#**********************************************************Dot plot class
def loadDotPlotPS(path,tag):
    positions=[]
    if(tag=="RNAfold"):
        positions=[0,1,2,3]
    if(tag=="RNAalifold"):
        positions=[3,4,5,6]
    res = {}
    outTmp = open(path)
    for l in outTmp:
        data = l[:-1].split()
        if len(data) == positions[3]+1 and data[positions[3]]=="ubox":
            i = int(data[positions[0]])-1
            j = int(data[positions[1]])-1
            p = math.pow(float(data[positions[2]]),2.)
            res[i,j] = p
    outTmp.close()
    return res

from os.path import isfile, join
from os import listdir
def Load_Probabilities(mypath):
    Dic = {}
    B = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda: 0)))
    Beta = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f.endswith('.proba')]
    for i, file in enumerate(Beta):
        Dic[i] = file.split(".proba")[0]
        lines = FF.Parsefile(mypath + "/" + file)
        for it in range(len(lines)):
            lines[it] = lines[it].split("\t")
            B[i][int(lines[it][0])][int(lines[it][1])] = float(lines[it][2])

    return Dic, B
class DotPlot:
    """Class for holding/producing base-pair probability matrices"""
    def __init__(self, rna , bpm = None):
        self.rna = rna[:]
        if bpm is None:
            # we will avoid this case to be sure that rnafold from the min works well
            self.bpm = self.runRNAFold()
        else:
            self.bpm = bpm
    def getSeq(self):
        return self.rna

    def getBPProb(self,i,j):
        if (i,j) in self.bpm:
            return self.bpm[i,j]
        else:
            return 0.

    def getUnpairedProbs(self):
        res = [1. for i in self.rna]
        for i,j in self.bpm:
            res[i] -= self.bpm[i,j]
            res[j] -= self.bpm[i,j]
        return res