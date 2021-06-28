#!/Users/keli/anaconda3/bin/python

import sys
import numpy as np
import XYZReader
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

class Connection:

    def __init__(self, boxsize, Otype = '63', Htype = '64'):
 
        self.dimensions = np.array([boxsize, boxsize, boxsize])
        self.Htype = Htype
        self.Otype = Otype

    def pbc_distance(self, coordsHydro, coordsOxy):
        """
        Get the distance between all H and all O in one frame
        """
        delta = np.abs(coordsHydro-coordsOxy[:,None]) # broadcasting 
        # periodic boundary conditions
        delta = np.where(delta > 0.5 * self.dimensions, delta - self.dimensions, delta)

        return np.sqrt((delta ** 2).sum(axis=-1))

    def sortedDisAndIndex(self, coords, atom_names):

#       NumOfH = len(self.coords[0][self.atom_names[0][:,0] == self.Htype]) 
        Ocoords = coords[tuple([atom_names[:,0] == self.Otype])]
        Hcoords = coords[tuple([atom_names[:,0] == self.Htype])]
        HcoordsIndex =np.where(atom_names[:,0] == self.Htype)
        OcoordsIndex = np.where(atom_names[:,0] == self.Otype)
        disMatrix = self.pbc_distance(Hcoords, Ocoords)
#           disMatrix =  [[H1-O1  H2-O1 H3-O1]
#                         [H1-O2  H2-O2 H3-O2]  
#                         [H1-O3  H2-O3 H3-O3]]
#           size: No. Oxy * No. H 
        indexSort = np.argsort(disMatrix, axis=0) # column are sorted
        disSort = np.sort(disMatrix, axis=0)

        HbondInfo = np.where(disSort[1] <= 2.5, True, False)

        HbondIndex = [idx if info else -1 for (info,idx) in zip(HbondInfo, OcoordsIndex[0][indexSort[1]])]
#        print(HcoordsIndex[0][1::2]-2)
        connectionInfo = zip(HcoordsIndex[0], HbondIndex)

        # Oxygen index: OcoordsIndex[indexSort[1]]
#        print(connectionInfo)
        return indexSort, disSort, OcoordsIndex[0], connectionInfo

#   def hbondNetwork(self, waterIndexRange, connectionInfo):
#
#        waterClusterDic = {}
#        while waterIndexRange != []:
#            start = waterIndexRange[0]
#            waterCluster = set()
#            waterCluster.add(start)
#            waterIndexRange.remove(start)
#            findFlag = True
#            n = 0
#            while findFlag:
#                for p in connectionInfo:
#                    if (start + 1 == p[0] and p[1] != -1) or (start + 2 == p[0] and p[1] != -1):
#                        waterCluster.add(p[1])
#                        start = p[1]
#                        try:
#                            waterIndexRange.remove(p[1])
#                        except ValueError:
#                            pass
#                        break
#                    else:
#                        findFlag = False
#                continue
#
#           waterClusterDic[n] = len(waterCluster)
#        return waterClusterDic
#        pass
                


def main():
    
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: HBondNetwork.py <file.xyz> boxsize \n")
        sys.exit(1)
    filename = sys.argv[1]
    boxsize = float(sys.argv[2])
    traj = XYZReader.XYZReader(filename)
    coords = traj.atomC
    nFrames = traj.nFrames
    atom_names = traj.atom_names
    connect = Connection(boxsize)
    size_distribution = []
    largest_list = []
 #   with open("hbondNetwork.txt", "w") as f:
 #       f.write("Source,Target\n")
    for i in range(nFrames):

        Hbond_list = []

            #f.write("Frame No.{}\n".format(i))

        (_, _, _, connectionInfo) = connect.sortedDisAndIndex(coords[i], atom_names[i])
    
        for p in connectionInfo:
            if p[1] != -1:
                Hbond_list.append([int(p[0]/3), int(p[1]/3)])
        Hbond_list = np.array(Hbond_list)
        node = list(set(Hbond_list[:,0]))
        G = nx.Graph()
        for p in Hbond_list:
            G.add_edge(p[0],p[1])
        for n in node:
            G.add_node(n)
        largest_list.append(len(max(nx.connected_components(G), key=len)))
        for size in nx.connected_components(G):
            size_distribution.append(len(size))
    #print sorted(size_distribution)
    sns.distplot(size_distribution,hist = False, kde = True)
                    #f.write("index{0},index{1}\n".format(p[0]/3 ,p[1]/3))
    plt.show()
    print("The largest HBond network size is {} +- {}".format(np.array(largest_list).mean(),np.std(largest_list)))
    np.save('water_cluster', np.array(largest_list))

if __name__ == '__main__':
    main()
