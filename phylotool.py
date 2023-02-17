#! /usr/bin/env python3

import argparse
import os
import sys
import pickle
import glob
import itertools
from Bio import Phylo

def readTaxonomy(tf):
    """ Function to read in the tab-separated RDP taxonomy database file and
    returns a dictionary with ID:TAX key:value pairs in string format. """
    taxDict = {}
    with open(tf, 'r') as f:
        for l in f.readlines():
            idx, tax = l.split("\t")
            taxDict[str(idx)] = str(tax)
    return(taxDict)

def cleanDatabaseDict(rdpD):
    """ Removes all entries without species or genus information from
    the RDP database dictionary """
    keys = []
    for key in rdpD:
        taxonomy = rdpD[key].split(";")
        species = taxonomy[-1]
        genus = taxonomy[-2]
        if  not genus.startswith("g__"):
            genus = "g__" + species.split("_")[2]
        if not species.startswith("s__"):
            keys.append(key)
    for key in keys:
        del rdpD[key]
    return(rdpD)

def collectGenusInfo(tr, rdpD):
    """ Collects all the genuses and its member species (clade) from the tree
    into a dictionary in the following format:Key: genus, Value:
    list(members) """
    gen = {}
    for clade in tr.get_terminals():
        idNum = clade.name.split("_")[0]
        taxonomy = rdpD[idNum].split(";")
        genus = taxonomy[-2].strip()
        if genus not in gen:
            gen[genus] = []
        gen[genus].append(clade)
    return(gen)


def getMinGenusDist(tr, gen):
    """ Calculates the minimum genus distance between the members of the family.
    If a family contains less than 2 member genuses, use the maximum species
    variation as the criteria instead. Input is a genus dictionary for a family
    and a tree file. """
    min_genus_dist = {}
    if len(list(gen)) < 2:
        genus = list(gen)[0]
        distL = []
        for sp_pair in itertools.combinations(gen[genus], 2):
            sp1, sp2 = sp_pair
            distL.append(tr.distance(sp1, sp2))
        return(max(distL))
    for genuspair in itertools.combinations(list(gen), 2):
        genus1, genus2 = genuspair
        distL = []
        for sp_pair in itertools.product(gen[genus1], gen[genus2]):
            sp1, sp2 = sp_pair
            dist = tr.distance(sp1, sp2)
            distL.append(dist)
        min_genus_dist[genuspair] = min(distL) # Get min from the species variation
    minimum_dist = min(min_genus_dist.values()) # Return the min of genus variation
    return(minimum_dist)


def terminalNeighborDists(self):
    """Return a list of tuples between adjacent consensus terminals.
    The tuple includes both terminal names and the distance inbetween
    """
    def generate_pairs(self):
        pairs = itertools.tee(self)
        next(pairs[1])
        return(zip(pairs[0], pairs[1]))
    dist_list = [(i[0].name, i[1].name, self.distance(*i)) for i in
                 generate_pairs(self.find_clades(terminal=True))]
                 #if "_C" in i[0].name
                 #]
    return(dist_list)


def clusterTreeNodes(rdp_taxonomy, \
                     tree, \
                     min_family_dist):

    """
    DESCRIPTION:
    This function clusters nodes in the given tree using \
    thresholds from min_family_dist.

    PARAMETERS:
    tree - A phylogenetic tree built from in-house \
    tool output
    
    OUTPUT:
    List of clusters, where each cluster is list \
    of clades beloning to the cluster

    """

    cladeClusters = []

    cladeSet = set(tree.find_clades())

    for terminalNode in tree.get_terminals():

        cluster = []
        for clade in cladeSet:
            if clade.name == None:
                continue

            dist = tree.distance(terminalNode, clade)
            idNum = clade.name.split("_")[0]
            taxonomy = rdp_taxonomy[idNum].split(";")
            familyTag = taxonomy[-3].strip().strip().split(" ")[0]
            if not familyTag.startswith("f__"):
                print ("Can not cluster " + taxonomy[-1] + ". No family info.")
                continue

            family = familyTag[3:]

            if dist < min_family_dist[family]:
                cluster.append(clade)

        cladeSet = cladeSet - set(cluster)
        if len(cluster) > 0:
            cladeClusters.append(cluster)

    return cladeClusters


def findRef(rdp_taxonomy, \
            tree, \
            min_family_dist):

    cladeSet = set([x for x in tree.find_clades() if x.name is not None])

    D = []
    for node_a in cladeSet:
        if node_a.name.endswith("_C"):

            l = []
            for node_b in cladeSet:
                dist = tree.distance(node_a, node_b)
                l.append((node_b, dist))

            for t in sorted(l, key=lambda k: k[1]):
                node_b, dist = t
                if not node_b.name.endswith("_C"):
                    D.append((node_a, node_b, dist))
                    break

    for t in sorted(D, key=lambda k: k[2]):
        node_a, node_b, dist = t

        idNum = node_a.name.split("_")[0]
        taxonomy = rdp_taxonomy[idNum].split(";")
        familyTag = taxonomy[-3].strip().strip().split(" ")[0]
        threshold = -1
        thresholdtxt = ""
        if familyTag.startswith("f__"):
            family = familyTag[3:]
            if family in min_family_dist:
                threshold = min_family_dist[family] 
                if dist < threshold:
                    thresholdtxt = ", distance is within the threshold " + str(threshold)
                else:
                    thresholdtxt = ", distance is NOT within the threshold " + str(threshold)



        print(node_a.name + " mapped to " + node_b.name + " with max distance: " + str(dist) + thresholdtxt)
                            



    return
    



if __name__=="__main__":

    parser = argparse.ArgumentParser(description='This is a script for ASA')

    parser.add_argument('--input-tree', 
                        action='store', 
                        dest='tree', 
                        required=True, 
                        help='Newick tree file')

    parser.add_argument('--RDP-taxonomy', 
                        action='store', 
                        dest='rdp_taxonomy', 
                        required=True, 
                        help='RDP taxonomy file.')

    parser.add_argument('--RDP-tree-dir', 
                        action='store', 
                        dest='rdp_tree_dir', 
                        required=False, 
                        help='RDP tree directory.',
                        default=False)

    args = parser.parse_args()
    tree = Phylo.read(args.tree, 'newick')
    rdp_taxonomy = cleanDatabaseDict(readTaxonomy(args.rdp_taxonomy))

    min_family_dist = {}
    if not os.path.isfile("min_family_dist.pickle"):
        print ("BUILDING min dists...")
        rdp_tree_dir = args.rdp_tree_dir

        for fil in glob.glob(rdp_tree_dir +'*.tree'):
            print ("Processing a tree " + fil)
            base_name = os.path.basename(fil)
            family = os.path.splitext(base_name)[0]
            fil = Phylo.read(fil, 'newick')
            genera = collectGenusInfo(fil, rdp_taxonomy)
            min_family_dist[family] = getMinGenusDist(fil, genera)
        
        if len(min_family_dist) == 0:
            print ("ERROR: No distances from trees at: " + rdp_tree_dir)
            sys.exit(1)

        with open("min_family_dist.pickle", "wb") as fil:
            pickle.dump(min_family_dist, fil)
    else:
        with open("min_family_dist.pickle", "rb") as fil: 
            print ("LOADING min dists...")
            min_family_dist = pickle.load(fil)
            print ("LOADING min dists, DONE.")


    findRef(rdp_taxonomy, tree, min_family_dist)

    print("done")

    sys.exit(0)
