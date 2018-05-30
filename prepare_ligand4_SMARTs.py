#!/usr/bin/python2.7 

import sys, os
import pybel
from itertools import chain
from os.path import basename
ob = pybel.ob

#Hits carboxylic acid, ester, ketone, aldehyde, carbonic acid/ester,anhydride, carbamic acid/ester, acyl halide, amide.
# defino nuevo tipo de oxigeno del carbonilo o cetona como "OC" antes llamado "OA"
smarts1 = "[OX1]=*"
#	
# 	defino nuevo tipo de nitrogeno dador puro como  "ND"  antes llamado "N"
smarts2 = "[N&H3,H2,H1]"
#	hydroxyl (includes alcohol, phenol)
#	define a new atom derived from ethanol as "OE"  before called  "OA"
smarts3 = "[OX2H]"


################################################################
################################################################
#	Define a dictionary with  Atomtype_New:[Atomtype_original,SMART]
SMARTS = {
"OC":["OA",smarts1],
"ND":["N",smarts2],
"OE":["OA",smarts3],
}
################################################################
#	Defino el procedimiento que busca los atomos a modificar con SMARTS y devuelve una lista con los seriales.
def find_smarts(data):
		#itero sobre las moleculas cargadas
        for mol in pybel.readfile("pdbqt",data):
			#creo una lista de smarts vacia para usar.
			output = []
			#itero sobre los smarts definidos anteriormente.
			for key in SMARTS:
				smarts = pybel.Smarts(SMARTS[key][1])          # define SMARTS query
				sfind = []
				sfind = list(chain.from_iterable(smarts.findall(mol)))                   # Matches an SMARTS
				#	verifico que no sea vacio sfind
				if sfind:
					temp = [[sfind] + [SMARTS[key][0]] + [key]]  # format output
					output = output + temp
			return output 
################################################################
################################################################
################################################################





################################################################
####################Clase que parsea el pdb#####################
################################################################

class PDBQT():
    def __init__(self, line):
        self._parse_common(line)     # used here (PDB) and in PDBQT
        self._parse_specific(line)   # differs in PDBQT

    def getline(self):
        txt = self._print_common()    # no \n; PDB + PDBQT
        txt += self._print_specific() # differs in PDBQT
        return txt

    def _parse_common(self, line):
        """Common to PDB and PDBQT formats"""
        self.keyword     = line      [ 0: 6]     # ATOM or HETATM
        self.serial      = int(line  [ 6:11])    # atom id
        #                            [11:12]
        self.name        = line      [12:16]     # atom name
        self.altLoc      = line      [16:17]     # Alternate location
        self.resName     = line      [17:20]     # Residue name
        #                            [20:21] 
        self.chain       = line      [21:22]     # chain
        self.resNum      = int(line  [22:26])    # Residue number
        self.icode       = line      [26:27]     # ???
        #                            [27:30]
        self.x           = float(line[30:38])    # X
        self.y           = float(line[38:46])    # Y
        self.z           = float(line[46:54])    # Z
        self.occupancy   = float(line[54:60])    # Occupancy
        self.bfact       = float(line[60:66])    # Temperature factor

    def _parse_specific(self, line):
        """ PDBQT characters [68:79] """
        self.charge      = float(line[68:76])   # Charge
        self.atype       = line      [77:79]    # Atom type
        self.atype = self.atype.strip().upper()

    def _print_common(self):
        """ Characters [0:68]"""
        linestr = ''
        linestr += '%6s' % (self.keyword)
        linestr += '%5d' % (self.serial)
        linestr += ' ' 
        linestr += '%4s' % (self.name)
        linestr += '%1s' % (self.altLoc) 
        linestr += '%3s' % (self.resName)
        linestr += ' ' 
        linestr += '%1s' % (self.chain)
        linestr += '%4d' % (self.resNum)
        linestr += '%1s' % (self.icode)
        linestr += ' ' * 3 
        linestr += '%8.3f' % (self.x)
        linestr += '%8.3f' % (self.y)
        linestr += '%8.3f' % (self.z)
        linestr += '%6.2f' % (self.occupancy)
        linestr += '%6.2f' % (self.bfact)
        return linestr

    def _print_specific(self):
        """ PDBQT characters [68:79] """
        linestr =  ' ' * 2                      # [66:68]
        linestr += '%8.3f' % (self.charge)      # [68:76]
        linestr += ' ' * 1                      # [76:77]
        linestr += '%-2s' % (self.atype)       # [77:79]
        #linestr += '\n'
        return linestr






################################################################
###########Prodecimiento  que modifica el pdbqt#################
################################################################
def mod_pdbqt(filename, outfilename):
	"""Creates a list of PDBQT atom objects"""
	#abro el archivo a modificar
	f = open(filename)
	#defino el nombre de salida
	salida = outfilename
	#en caso de no haber dado archivo de salida seteamos el default
	if salida is None:
		salida = os.path.splitext(basename(filename))[0]+"_biase.pdbqt"
	#abro el archivo de salida
	out = open(salida, 'w')
	#Escribo el header de salida referenciando la molecula original.
	out.write("REMARK  Name = " + filename + "\n"),
	#itero sobre las lineas del archivo.
	for line in f:
		if line.startswith('ATOM  ') or line.startswith('HETATM'):
			atom = PDBQT(line)
			#Busco SMARTS sobre un archivo. 
			smarts_list = find_smarts(filename)
			#print smarts_list
			#chequeo que la lista no este vacia.
			if smarts_list:
				#itero sobre los smarts que debo cambiar
				for s in smarts_list:
					#lista de seriales a buscar
					serials = s[0]
					#tipo de atomo a cambiar
					oldtype = s[1]
					#nuevo tipo de atomo a reemplazar
					newtype = s[2]
					#itero sobre los posibles seriales.
					for serialnumber in serials:
						#selecciono el serial que ademas tiene el tipo de atomo a reemplazar
						if atom.serial == serialnumber and atom.atype == oldtype:
							#seteo el nuevo tipo de atomo.
							atom.atype = newtype
			#imprimo las lineas hayan o no sido modificadas.
			out.write(str(atom.getline()) + "\n"),
		else:
			#imprimo las lineas relacionadas con los REMARKS o torsiones ENDS etc.
			out.write(str(line)),
	# cierro el archivo finalmente.
	f.close
	out.close()




################################################################
#####Ejecuto el programa tomando como archivo de entrada el##### 
#####primer argumento por consola###############################
#print "REMARK  Name = ", sys.argv[1]
################################################################



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: prepare_ligand4_SMARTs.py -l filename"
        print
        print "    Description of command..."
        print "         -l     ligand_filename (.pdbq format)"
        print "    Optional parameters:"
        print "        [-o pdbqt_filename] (default output filename is ligand_filename_stem + _biase + .pdbqt)\n\n\n"
        print "							Develop by Elias Lopez (2015) contact:eliaslopez at qb.fcen.uba.ar\n"


    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'l:o:')
    except getopt.GetoptError, msg:
        print 'prepare_ligand4_SMARTs.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-l: ligand
    ligand_filename =  None
    # optional parameters
    #-o outputfilename
    outputfilename = None

    #'l:vo:d:A:CKU:B:R:MFI:Zgs'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-l', '--l'):
            ligand_filename = a
        if o in ('-o', '--o'):
            outputfilename = a
        if o in ('-h', '--'):
            usage()
            sys.exit()
    if not  ligand_filename:
        print 'prepare_ligand4: ligand filename must be specified.'
        usage()
        sys.exit()
    else:
		# ejecuto el programa propiamente dicho.
		mod_pdbqt(ligand_filename, outputfilename) 

