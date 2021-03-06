import pickle as pl
import matplotlib.pyplot as plt1
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn import manifold
from sklearn import cluster, datasets
import plotly
plotly.tools.set_credentials_file(username='afamine', api_key='2znk0yep6s')
import plotly.graph_objs as go
import plotly.plotly as py
from mpl_toolkits.mplot3d import axes3d # To resove the problem of  projection='3d'
from itertools import cycle
import scipy, numpy as np,os
import conf, FileFunctions
from sklearn import manifold
from collections import defaultdict

import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

# Visualize the distribution of the structures into clusters, visualize the MFES
# The 2D plan is the distance matrix [sklearn is called to plot the distribution], the Z axis is the Boltzmannprobabilty
# The MFEs structures are filtred by the end of this process
def ThreeD_MatDistance_Boltzmann(MatDist, Klust, Boltzmannprobabilty, numberofsruct, constrainte,MFEs):
    #clusters=defaultdict(aa)
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    fig_handle = plt1.figure()
    ax1 = fig_handle.add_subplot(111, projection='3d')
    # plot 3 d
    elem = [i for i in range(len(MatDist))]
    D = scipy.zeros([len(MatDist), len(MatDist)])
    for i in range(len(MatDist)):
        for j in range(len(MatDist)):
            D[i][j] = MatDist[i][j]
    adist = np.array(D)
    amax = np.amax(adist)
    adist /= amax
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    coords = results.embedding_
    plt1.subplots_adjust(bottom=0.1)
    for k, col in zip(range(len(Klust)), colors):
        pos = -1
        for i in Klust[k]:
            #print i
            pos += 1
            if i - 1 < numberofsruct * (len(constrainte) - 1):
                #print "passe test"
                ConditionNumber = int((i - 1) / numberofsruct)
                StructureNumber = i - 1 - ConditionNumber * numberofsruct
                #print "ppppppppppp",ConditionNumber,
                ax1.scatter(coords[i - 1, 0], coords[i - 1, 1],
                            Boltzmannprobabilty[constrainte[ConditionNumber]][StructureNumber], c=col, marker='.')
            else:
                #print i-1,"uhuhuh"
                ConditionNumber = len(constrainte) - 1
                StructureNumber = i - 1 - ConditionNumber * numberofsruct
                #print StructureNumber,"gg"
                #print "Cluster number", k, "contains MFE", MFEs[StructureNumber]
                ax1.text(coords[i - 1, 0], coords[i - 1, 1],
                         Boltzmannprobabilty[constrainte[ConditionNumber]][StructureNumber],
                         '*%s' % (MFEs[StructureNumber]), color=col)
                #print Klust[k][pos], "avant avant??"
                #Klust[k][pos]=" "
                #del Klust[k][pos]
                #Klust[k]= {v for v in Klust[k].items() if v=="M"}
                #print "step", Klust
                #pos-=1
                #clusters[k][pos] = ""
                #print Klust[k][pos],"ca fonctionne pas??"
    ax1.set_xlabel('MDS axis1')
    ax1.set_ylabel('MDS axis2')
    ax1.set_zlabel('Boltzmann probability')
    ax1.set_title(' Secondary Structures Multidimensional Scaling ')
    fig_handle.savefig('output_plots/MDS_Structures_MFES.svg')
    #plt1.show()
    #FileFunctions.PickleVariable(clusters, os.path.join(conf.PickledData, "Clusters_WIthoutMFES.pkl"))
    return 0

#To eminitae the MFEs structures
def FilterClusters(Klust, lenconst):
    for C in range(len(Klust)):
        #to eliminate MFEs
        Klust[C] = [v for  v in Klust[C] if v<lenconst]
        # if the cluster becomes empty , delete it
        if Klust[C]==' ':
            del(Klust[C])
    return Klust

# To resolve the problem of pickling a defaultdict!!
def dd():
    return 0
def aa():
    return defaultdict(dd)

def threedcentoids(MatDist, Centroids_Energies, ListDiameters):
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    fig_handle = plt1.figure()
    ax = fig_handle.add_subplot(111, projection='3d')

    D = scipy.zeros([len(MatDist), len(MatDist)])
    for i in range(len(MatDist)):
        for j in range(len(MatDist)):
            D[i][j] = MatDist[i][j]
    # print "hhhhh",len(MatDist)
    adist = np.array(D)
    amax = np.amax(adist)
    adist /= amax
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    coords = results.embedding_
    plt1.subplots_adjust(bottom=0.1)
    for k, col in zip(range(len(MatDist)), colors):
        ax.scatter(coords[k, 0], coords[k, 1], Centroids_Energies[k], c=col, marker='.')
        # ax.add_patch(mpatches.Circle((coords[k, 0], coords[k, 1]),ListDiameters/2,color=col,edgecolor="black"))
        ax.text(coords[k, 0], coords[k, 1], Centroids_Energies[k], 'C%s' % (k + 1), color=col)

    ax.set_xlabel('MDS axis1 (Base pairs distance)')
    ax.set_ylabel('MDS axis2')
    ax.set_zlabel('Boltzmann energy Centroids')
    ax.set_title('Clusters distances with centoid s Boltzmann energies')
    # pl.dump(fig_handle,file('sinus.pickle','wb'))
    fig_handle.savefig('centroids_distribution.svg')
    plt1.show()


def plotDistanceClusters(D, clusters, coloro,title):
    # D=[ [0, 1,2,4], [1, 0,3,2],[2,3,0,5],[4,2,5,0] ]
    Dic = {}
    for elem in clusters:
        Dic[elem] = elem
    adist = np.array(D)
    amax = np.amax(adist)
    adist /= amax
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    coords = results.embedding_
    # plot results
    fig = plt.figure()
    plt.subplots_adjust(bottom=0.1)
    plt.scatter(coords[:, 0], coords[:, 1], marker='o')
    for label, x, y in zip(Dic.values(), coords[:, 0], coords[:, 1]):
        plt.annotate(
            label,
            xy=(x, y), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.2', fc=coloro, alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    fig.savefig('Distance_clusters_' + title + '.png')

'''
def plotClustercBECard(clusternumber,cBE,cardinal,xlabelo,ylabelo,output):
	Labels=clusternumber
    X=cBE
	Y=cardinal
 	fig = plt.figure()
	fig.suptitle('Pareto front for clusters', fontsize=14, fontweight='bold')
	ax = fig.add_subplot(111)
	fig.subplots_adjust(top=0.85)
	ax.set_title('Cluster distribution')
	ax.set_xlabel(xlabelo)
	ax.set_ylabel(ylabelo)
    ax.set_ylim(bottom=0)
    plt.axis([0,max(X)+1 ,0, max(Y)+np.mean(Y)])
    ax.grid(True)
	for i in Labels:
		ax.text(X[i]+0.2,Y[i], i+1,fontsize=10,horizontalalignment='center',color='b')
    plt.plot(X,Y, 'r*')

'''


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
def plotPareto(paretoPoints,dominatedPoints):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	dp = np.array(list(dominatedPoints))
	pp = np.array(list(paretoPoints))
	print(pp.shape,dp.shape)
	ax.scatter(dp[:,0],dp[:,1],dp[:,2])
	ax.scatter(pp[:,0],pp[:,1],pp[:,2],color='red')

	import matplotlib.tri as mtri
	triang = mtri.Triangulation(pp[:,0],pp[:,1])
	ax.plot_trisurf(triang,pp[:,2],color='mediumvioletred')
	plt.show()


from matplotlib import pyplot as PLT
from matplotlib import cm as CM
from matplotlib import mlab as ML
import numpy as NP


def plotPairs(lista, n):
    fig = PLT.figure()
    x = [elem[0] for elem in lista]
    y = [elem[1] for elem in lista]
    z = [elem[2] / float(n) for elem in lista]
    gridsize = 60
    PLT.hexbin(x, y, C=z, gridsize=gridsize, cmap=CM.jet, bins=None)
    PLT.axis([min(x) - 1, max(x) + 1, min(y) - 1, max(y) + 1])
    cb = PLT.colorbar()
    cb.set_label('Probability value  in all  optimal centroids')
    PLT.show()


from matplotlib.pyplot import figure, show
import numpy
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def plotClustercBECard(clusternumber,cBE,cardinal,xlabelo,ylabelo,output):

	Labels=clusternumber
        X=cBE
	Y=cardinal
 	fig = plt.figure()
	fig.suptitle('Pareto front for clusters', fontsize=14, fontweight='bold')
	ax = fig.add_subplot(111)
	fig.subplots_adjust(top=0.85)
	ax.set_title('Cluster distribution')
	ax.set_xlabel(xlabelo)
	ax.set_ylabel(ylabelo)
        ax.set_ylim(bottom=0)
        plt.axis([0,max(X)+1 ,0, max(Y)+np.mean(Y)])
        ax.grid(True)
	for i in Labels:
		ax.text(X[i]+0.2,Y[i], i+1,fontsize=10,horizontalalignment='center',color='b')
        plt.plot(X,Y, 'r*')
        fig.savefig(output)


import StructureFunctions as SF
def plotsecodnarystructures(rnaString, lista, lista2, n, reactivities):
    fig = plt.figure()
    fg, ax = plt.subplots(1, 1)
    import matplotlib as mp

    min_val = 0
    max_val = 1
    my_cmap = cm.get_cmap('Greys')
    norm = matplotlib.colors.Normalize(min_val, max_val)  #
    x = [elem[0] for elem in lista]
    y = [elem[1] for elem in lista]
    z = [elem[2] / float(n) for elem in lista]

    x2 = [elem[0] for elem in lista2]
    y2 = [elem[1] for elem in lista2]
    z2 = [elem[2] / float(n) for elem in lista2]

    rna = list(rnaString)
    cc = ["black" for i in range(len(rna))]
    for elem in range(len(rna)):
        if float(reactivities[elem]) < 0.2:
            cc[elem] = "black"  # "#00509d"#blue
        if float(reactivities[elem]) >= 0.2 and float(reactivities[elem]) < 0.4:
            cc[elem] = "gray"  # "#00c200"#green
        if float(reactivities[elem]) >= 0.4 and float(reactivities[elem]) < 0.7:
            cc[elem] = "mediumvioletred"  # "#f28f00"#yellow
        if float(reactivities[elem]) >= 0.7:
            cc[elem] = "plum"  # "#f20000"#red

    p0 = Rectangle((0, 0), 1, 1, fc="black")
    p1 = Rectangle((0, 0), 1, 1, fc="gray")
    p2 = Rectangle((0, 0), 1, 1, fc="mediumvioletred")
    p3 = Rectangle((0, 0), 1, 1, fc="plum")
    ax.legend([p0, p1, p2, p3], ["Reactivity <0.2", "0.2< <0.4", "0.4<  <0.7", ">0.7"])

    pac = [mpatches.Arc([x[i] + 0.5 + (y[i] - x[i] - 1) / float(2), 0], y[i] - x[i] - 1, (y[i] - x[i] - 1) / float(2),
                        angle=0, theta1=0, theta2=180, color=my_cmap(norm(z[i])), linewidth=1) for i in
           range(len(x))]  # linestyle='dotted', linestyle='dashed
    pac2 = [mpatches.Arc([x2[i] + 0.5 + (y2[i] - x2[i] - 1) / float(2), 0], y2[i] - x2[i] - 1,
                         (y2[i] - x2[i] - 1) / float(2), angle=0, theta1=180, theta2=360, color=my_cmap(norm(z[i])),
                         linewidth=1) for i in range(len(x2))]
    for arc, arc2 in zip(pac, pac2):
        ax.add_patch(arc)
        ax.add_patch(arc2)
    cmmapable = cm.ScalarMappable(norm, my_cmap)
    cmmapable.set_array(range(min_val, max_val))
    colorbar(cmmapable, fraction=0.046, pad=0.04, ticks=[0, 0.5, 1])

    fontProp = mp.font_manager.FontProperties(family="monospace", style="normal", weight="bold", size="8")
    ax.axis([0, max(x) + 20, -max(y) / 3, max(y) / 3])
    for i in range(len(rna)):

        nuc = rna[i]
        ax.add_patch(
            mpatches.Circle((i + 0.5, 0), 0.5, color=cc[i], edgecolor="black"))  # circle at center (x,y), radius 0.5
        ax.annotate(nuc, (i + 0.5, 0), color='white', weight='bold', fontsize=6, ha='center', va='center')
        ax.annotate(i + 1, (i + 0.5, -1), color='black', weight='bold', fontsize=6, ha='center', va='center')
    ax.set_aspect("equal")
    ax.get_yaxis().set_visible(False)
    fg.canvas.draw()

    ax.set_title('Combined optimal centroids ')

    plt.show()
    plt.savefig("res.eps", format='eps', dpi=1000)
    fig.savefig('Arcs_structures.svg')

def plotClusteringDistribution(lenconstraint, Folder_name, Lenrna):
        D = scipy.zeros([lenconstraint, lenconstraint])
        Dic, B = SF.Load_Probabilities(Folder_name)
        #print Dic
        # calculate the Eucledian distance between different matrix
        D = SF.Eucledian_distance(B, Lenrna)
        print "Distance", D
        # Clustering process with th plot
        adist = np.array(D)
        amax = np.amax(adist)
        adist /= amax
        mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
        results = mds.fit(adist)
        coords = results.embedding_
        # plot results
        fig = plt.figure()
        plt.subplots_adjust(bottom=0.1)
        plt.scatter(coords[:, 0], coords[:, 1], marker='o')
        for label, x, y in zip(Dic.values(), coords[:, 0], coords[:, 1]):
            plt.annotate(
                label,
                xy=(x, y), xytext=(0, 20),
                textcoords='offset points', ha='left', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.5),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        # plt.show()
        fig.savefig('Eucledian_distance_dot_plot_Matrix.png')

import FileFunctions as FF

from pyx import *
import StructureFunctions as SF
def Drawvarna(File1, Listclusters, CentroidStructure, numberofsruct, rna, Centroids_Energies,FileShape):
    c = canvas.canvas()
    i = 0
    with open(File1, "w") as OPt:
        lista = []
        for elem in Listclusters:
            i = i + 1
            filename = File1 + numberofsruct + "_" + str(elem)
            OPt.write("%s\n" % (CentroidStructure[elem]))
            # Heatmap
            lista += SF.ListBasePairsFromStruct(CentroidStructure[elem])

            # Varna call
            with open(filename, "w")as outt:
                outt.write(">HIV\n%s\n%s" % (rna, CentroidStructure[elem]))
            print " HIV  structure Centroid of the ", File1, elem, "with the energy value of", Centroids_Energies[
                elem - 1], CentroidStructure[elem]
            drawStructure(filename,FileShape , filename + ".eps")
            c.insert(epsfile.epsfile(0, i * 40, filename + ".eps"))
    return c, lista
def drawStructure(Sequence, Shapefile, outfile):
    lines = ""
    f = FF.Parsefile(Shapefile)
    lines = ";".join(["%.3f" % float(line.strip().split('\t')[1]) for line in f])
    cmd = 'java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -i ' + Sequence + ' -colorMap ' + '"' + lines + '"' + ' -colorMapStyle ' + '"-5.00:#4747B6,0.00:#4747FF,0.40:#1CFF47,0.70:#FF4747,2.41:#FFFF00"' + ' -algorithm line -o ' + outfile + " > /dev/null"
    os.system(cmd)