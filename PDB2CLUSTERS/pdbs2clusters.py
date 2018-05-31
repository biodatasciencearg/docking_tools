#!/usr/bin/env python
import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
import os, time
#from Bio import PDB
from Bio.PDB.PDBList import PDBList
from scipy.spatial.distance import squareform, pdist
import pymol
import numpy as np
import scipy
pymol.finish_launching()

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f s' % (f.func_name, (time2-time1))
        return ret
    return wrap

#Make pdblist from file.
def file2lst (filename):
	#return list
	rl=[]
	#chains sublist
	chains = []
	#ids sublist
	ids=[]
	with open(filename) as file:
    		for line in file:
			column = line.split() 
        		#line = line.strip() #or some other preprocessing
        		ids.append(column[0].lower()) #storing everything in memory!
			chains.append(column[1].upper())
	return (ids , chains) 
# Download a list of pdbs into a data directory.
def download_pdblst(pdbslst):
	pdblist = PDBList()
	for pdb in pdbslst:
		pdblist.retrieve_pdb_file(pdb, file_format="pdb", pdir='./data/')

# Align two Proteins 
# chain1 belogns to referece while the chain2 to a mobile selection. 
def Align (reference, mobile, amask, chain1, chain2):
	if not amask:
		amask=''
	else:
		amask=' and ' + amask
	# Load Structures
	#load mobile
	pymol.cmd.load('./data/pdb' + mobile + '.ent'  , mobile)
	#load reference.
	pymol.cmd.load('./data/pdb' + reference + '.ent'  , reference)
        pymol.cmd.do("cealign %s, %s" %(reference + " and chain " + chain1  + amask, mobile + " and chain " + chain2 ))
        time.sleep(1)
	#Save  structures.
	# save(file, selection, state (0 default), format)
	pymol.cmd.save("./data/%s.pdb" %(mobile + "_" + chain2 + '_aln'), mobile + " and chain " + chain2, 0, 'pdb')
	pymol.cmd.save("./data/%s.pdb" %(reference + "_" + chain1 + '_aln'), reference + " and chain " + chain1, 0, 'pdb')
	pymol.cmd.delete('all')

# Structural alignment of pdb list.
def AlignLst (pdbslst,amask):
	codes=pdbslst[0]
	chains=pdbslst[1]
	for code, chain  in zip(codes[1:], chains[1:]):
		# align two pdb files
		print(codes[0], code, amask,chains[0],chain)
		Align(codes[0], code, amask,chains[0],chain)
	print("\n\nINFO:All the pdbs files were saved in the data folder with the extension *_aln.pdb\n")
# get coordinates of water molecules in a list o pdbs.	
def get_waters(pdbslst,return_mask):
	codes=pdbslst[0]
        chains=pdbslst[1]
	from pymol import stored
	#load all pdb files
	for code, chain in zip(codes, chains):
		# Load Structures
        	pymol.cmd.load('./data/' + code + "_" + chain + '_aln.pdb'  , code + "_" + chain)
	#if return mask is all return all sites.
	if return_mask=="all":
		selection_mask = 'solvent'
        	pymol.cmd.select('waters', selection_mask)
	else:
		# reference residue.
		reference=codes[0] + "_" + chain + ' and ' + str(return_mask)
		pymol.cmd.select('S1', reference)
		pymol.cmd.select('S2', 'solvent')		
		selection_mask = 'S2 within 5 of S1'
		pymol.cmd.select('waters', selection_mask)
		time.sleep(2)
	pos = []
	pymol.cmd.iterate_state(1, 'waters', 'pos.append([x,y,z])', space={'pos': pos})
	time.sleep(2)
	pymol.cmd.save("./data/waters.pdb", 'waters', 0, 'pdb')
	pymol.cmd.delete('all')
	pymol.cmd.quit()
	return pos
#auxiliar fuction, to get indexes from cluster object.	
def cluster_indices(cluster_assignments):
    n = cluster_assignments.max()
    indices = []
    for cluster_number in range(1, n + 1):
        indices.append(np.where(cluster_assignments == cluster_number)[0])
    return indices


# compute geometric center.
def centroid(XYZ,cluster_indices):
	#coordinates of water molecules of cluster n.
	clusterXYZ=[]
	#Iterate over indexes.
	for index in cluster_indices:
		clusterXYZ.append(XYZ[index])
	return np.mean(clusterXYZ, axis=0)
#Function to get the clusters. Arguments, XYZ array and distance threshold.
def Aglomerative_clustering(XYZ,t,occupancy_threshold):
	print('Starting Clustering Process using Aglomerative method...')
	import scipy.cluster.hierarchy as hac
	clusters_assignments=hac.fclusterdata(XYZ,t,criterion='distance',method='single')
        #create a list o indexes.
        indices = cluster_indices(clusters_assignments)
	#output list
        out=[]
        #iterate over clusters.
        for k, ind in enumerate(indices):
                geom_center=centroid(XYZ,ind)
                #R90(XYZ,ind) Funcion a ser desarrolada.
                occupancy=len(ind)
                if int(occupancy) >= int(occupancy_threshold):
                        info_cluster= (geom_center[0],geom_center[1],geom_center[2],occupancy)
                        out.append(info_cluster)
        #sort the results by occupancy.
        out.sort(key=lambda x: x[3], reverse=True)
	return out
# another type of clustering method.
def DBSCAN_clustering(XYZ,epsilon,minsamples):
	occupancy_threshold=minsamples
	X = np.empty(shape=[0, 3])
	coords=[]	
	for coord in XYZ:
		coordenadas= np.array([coord[0], coord[1],coord[2]])
		coords.append(coordenadas)
	X=np.append(X, coords, axis=0)
	from sklearn.cluster import DBSCAN
	from sklearn import metrics
	from sklearn.preprocessing import StandardScaler
	print('Starting Clustering Process using DBSCAN method...')
	print'We found {} points ...'.format(len(X))
	db=DBSCAN(eps=float(epsilon), min_samples=int(minsamples)).fit(X)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
	#Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

        print('Estimated number of clusters: %d' % n_clusters_)
        #print("Silhouette Coefficient: %0.3f"
	#% metrics.silhouette_score(X, labels))
        unique_labels = set(labels)
	out=[]
	for k in unique_labels:
		K=int(k)
        	#Solo itero sobre los clusters mayores que -1.
                if K is not -1:
                        clusterNumb=K+1
                        class_member_mask = (labels == k)
                        xy = X[class_member_mask & core_samples_mask]
                        # Calculo el centroide para este mismo.
                        centroid=np.mean(xy, axis=0)
                        # occupancy
                        occupancy= len(xy)
			if int(occupancy) >= int(occupancy_threshold):
				info_cluster= (centroid[0],centroid[1],centroid[2],occupancy) 
				out.append(info_cluster)
	#sort the results by occupancy.
        out.sort(key=lambda x: x[3], reverse=True)
	return out
# New clustering method.
def HDBSCAN_clustering(XYZ,min_clust_size,occupancy_threshold):
	import hdbscan
	clusterer = hdbscan.HDBSCAN(metric='euclidean')
	X = np.empty(shape=[0, 3])
        coords=[]
        for coord in XYZ:
		coordenadas= np.array([coord[0], coord[1],coord[2]])
		coords.append(coordenadas)
	X=np.append(X, coords, axis=0)
	print('Starting Clustering Process using HDBSCAN method...')
	print'We found {} points ...'.format(len(X))
	clusterer.fit(X)
	        #create a list o indexes.
        indices = cluster_indices(clusterer.labels_)
        #output list
        out=[]
        #iterate over clusters.
        for k, ind in enumerate(indices):
                geom_center=centroid(XYZ,ind)
                #R90(XYZ,ind) Funcion a ser desarrolada.
                occupancy=len(ind)
                if int(occupancy) >= int(occupancy_threshold):
                        info_cluster= (geom_center[0],geom_center[1],geom_center[2],occupancy)
                        out.append(info_cluster)
        #sort the results by occupancy.
        out.sort(key=lambda x: x[3], reverse=True)
        return out



	




	 
# compute r90. (Funcion a ser desarrollada)
def R90(XYZ,cluster_indices):
	#get center.
	center = centroid(XYZ,cluster_indices)
	center_vect = [(center[0],center[1],center[2])]
	# get coodinates.
	clusterXYZ= []
        #Iterate over indexes.
        for index in cluster_indices:
		vector=(XYZ[index][0],XYZ[index][1],XYZ[index][2])
                clusterXYZ.append(vector)
	distance_matrix= scipy.spatial.distance.cdist(clusterXYZ,center_vect)
	print(distance_matrix)






#Main fuction to get a water clusters.
#Arguments filename with pdb codes, download var to set if download or not pdb files, aling variable to set if align the files and t is the threshold is a distance criterion to make the clusters.
@timing
def get_clusters (filename, download, align, amask, clustering, method, distance_threshold, epsilon, minsamples, min_cluster_size, return_mask, return_minimal_occupancy):
	#Generate pdb list from file.
	pdbslst=file2lst(filename)
	#Download pdbs
	if download:
		download_pdblst(pdbslst[0])
	#align pdbs.
	if align:
		#align pdbs.
		AlignLst(pdbslst,amask)
	if not clustering:
		print('\nDone. You chose not to do the clustering analysis\n\n')
		sys.exit()
	# Get Water Molecules.	
	print('	Loading Water Molecules.\n')
	XYZ=get_waters(pdbslst,return_mask)
	#Each number is the cluster number to which the water molecule belongs.
	print('\nDone.\nClustering calculation ...')
	#Make the clustering
	#choose method.
	if method=='Agglomerative':
		# Agglomerative Clustering. http://www.saedsayad.com/clustering_hierarchical.htm.
        	# We choose Single Linkage.
		clusters=Aglomerative_clustering(XYZ,distance_threshold,return_minimal_occupancy)
	elif method=='DBSCAN':
		clusters=DBSCAN_clustering(XYZ,epsilon, minsamples)
	elif method=='HDBSCAN':
		clusters=HDBSCAN_clustering(XYZ,min_cluster_size,return_minimal_occupancy)
	
	#clusters counter.
	cluster_number=0	
	out_final= []
	#iterate over clustes a add format.
	for element in clusters:
		#increment  cluster number.
		cluster_number+=1
		info_cluster= '%-5s %10.3f %10.3f %10.3f %10i %10i'  % ('O' + str(cluster_number),element[0],element[1],element[2],cluster_number,element[3])
		out_final.append(info_cluster)
	out_header=['%i' % (cluster_number)]
	out_header.append('ATOM        x       y       z        ClusterNumber   occupancy')
	# Open a file to save the output.
	outfile = open('cluster_cg.xyz', 'w')
	outfile.writelines(["%s\n" % item  for item in out_header])
	outfile.writelines(["%s\n" % item  for item in out_final])
	outfile.close()
	print('\nDone. The calculation of the sites has ended successfully. They were saved in the file clusters_cg.xyz\n\n')


def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line
def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError # evil ValueError that doesn't tell you what the wrong value was
	
def parsefile(filename):
	#set defaul parameters.
	pdb_list=None
	download = True
	align = True
	align_mask = None
	clustering=True
	method='HDBSCAN'
	distance_threshold=2.5
	epsilon=0.80
	minsamples=2
	min_cluster_size=2
	return_mask='all'
	return_minimal_occupancy=2
	#iterate over lines.	
	with open(filename) as inf:
                for line in nonblank_lines(inf):
                        if not line.startswith("#"):
				parameters=line.split("=")
				#setting new parameters.
				if parameters[0]=="pdb_list":
					pdb_list=parameters[1]
				elif parameters[0]=="download":
					download = str_to_bool(parameters[1])
				elif parameters[0]=="align":
					align = str_to_bool(parameters[1])
				elif parameters[0]=="align_mask":
					align_mask=parameters[1]
				elif parameters[0]=="clustering":
					clustering=str_to_bool(parameters[1])
				elif parameters[0]=="method":
					method=parameters[1]
				elif parameters[0]=="distance_threshold":
					distance_threshold=parameters[1]
				elif parameters[0]=="epsilon":
					epsilon=parameters[1]
				elif parameters[0]=="minsamples":
					minsamples=parameters[1]
				elif parameters[0]=="min_cluster_size":
					min_cluster_size=parameters[1]
				elif parameters[0]=="return_mask":
					return_mask=parameters[1]
				elif parameters[0]=="return_minimal_occupancy":
					return_minimal_occupancy=parameters[1]
	return (pdb_list, download, align, align_mask, clustering, method, distance_threshold,epsilon,minsamples,min_cluster_size,return_mask,return_minimal_occupancy)

##############
#####MAIN##### 
##############



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: ./pdbs2clusters.py -l filename"
        print
        print "    Description of command..."
        print "         -f    input filename"
        print "    parameters:"
	print "        <list_files> pdbid list \n"
        print "        <download> (default=True) \n"
	print "        <align> pdb files (default=True) \n"
	print "        <align_mask> optional mask to align (default=None) \n"
	print "        <clustering> optional parameter to make clustering (default=True) \n"
	print "        <method> Clustering methods: aglomerative, HDBSCAN or DBSCAN methods (default=HDBSCAN)\n"
	print "        <return_minimal_occupancy>  (default=2 molecules) \n"	
	print "        <return_mask> returns  water sites surrounding the active site defined by pymol mask (default=all)"
        print "							Developed by Elias Lopez (2017) contact:eliaslopez@qb.fcen.uba.ar\n"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'f:')
    except getopt.GetoptError, msg:
        print 'pdbs2clusters: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-l: filename
    filename =  None
    #-o Download or not.
 
    #'l:vo:d:A:CKU:B:R:MFI:Zgs'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-f', '--f'):
            filename = a
	if o in ('-h', '--'):
            usage()
            sys.exit()
    if not  filename:
        print 'pdbs2clusters: filename with input parameters must be specified.'
        usage()
        sys.exit()
    else:
		# parse input file
		pdb_list, download, align, amask, clustering, method, distance_threshold, epsilon, minsamples, min_cluster_size, return_mask, return_minimal_occupancy = parsefile(filename)
		if not pdb_list:
			print 'pdbs2clusters: pdb_list  with pdb codes must be specified.'
		        usage()
        		sys.exit()
		#
		#print(pdb_list, download, align, amask, clustering, method, distance_threshold, epsilon, minsamples, min_cluster_size, return_mask, return_minimal_occupancy)
		get_clusters(pdb_list, download, align, amask, clustering, method, distance_threshold,epsilon,minsamples,min_cluster_size,return_mask,return_minimal_occupancy) 

