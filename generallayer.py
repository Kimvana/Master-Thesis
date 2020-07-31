

## ------------- Dont edit anything between here and next line! ----------


import numpy as np
import inspect

def initialize(ion): ## create different variables storing the ion's xyz data
	## two challenges here: (1) this function has to create new variables, 
	## (2) use the name of input variable as a string 
	## the readme contains the 'original' version, instead of this function
	callers_global_vars=inspect.currentframe().f_back.f_globals.items()
	string=[key for key, value in callers_global_vars if value is ion][0]
	globals()[string+"_coords"] = np.genfromtxt(Molecule_folder+ion+Molecule_handler, skip_header=2, usecols=(1,2,3))
	globals()[string+"_coords"]/=10
	globals()[string+"_types"] = np.genfromtxt(Molecule_folder+ion+Molecule_handler, skip_header=2, usecols=0, dtype=str)
	globals()[string+"_coords"].flags.writeable = False
	globals()[string+"_types"].flags.writeable = False
	globals()[string+"_atomcount"] = globals()[string+"_types"].shape[0]
	
	## convert 'only atom type' to list of numbered types ready for .gro formation
	globals()[string+"_labeldict"] = {}
	globals()[string+"_labels"] = []
	for atom in globals()[string+"_types"]:
		if atom in globals()[string+"_labeldict"]: globals()[string+"_labeldict"][atom]+=1
		else: globals()[string+"_labeldict"][atom]=1
		globals()[string+"_labels"].append(atom+str(globals()[string+"_labeldict"][atom]))

def write_atom(atomtype):
	global Unit_Counter
	global Atom_Counter
	Output_gro.write(str(Unit_Counter).rjust(5)+"Unit "+atomtype.ljust(5)+str(Atom_Counter).rjust(5)+f"{xcoord:.3f}".rjust(8)+f"{ycoord:.3f}".rjust(8)+f"{zcoord:.3f}".rjust(8)+"  0.0000  0.0000  0.0000\n")
	Output_xyz.write(atomtype+" "+f"{10*xcoord:.3f}".rjust(8)+f"{10*ycoord:.3f}".rjust(8)+f"{10*zcoord:.3f}".rjust(8)+"\n")
	Atom_Counter += 1
	atomtypeS=atomtype+"S"
	Output_gro.write(str(Unit_Counter).rjust(5)+"Unit "+atomtypeS.ljust(5)+str(Atom_Counter).rjust(5)+f"{xcoord:.3f}".rjust(8)+f"{ycoord:.3f}".rjust(8)+f"{zcoord:.3f}".rjust(8)+"  0.0000  0.0000  0.0000\n")
	Atom_Counter += 1
	Unit_Counter += 1
	
def write_molecule(moltype):
	global mol_atom_count
	global Atom_Counter
	global xcoord
	global ycoord
	if xcoord > BoxX: xcoord-=BoxX
	if xcoord < 0: xcoord+=BoxX
	if ycoord > BoxY: ycoord-=BoxY
	if ycoord < 0: ycoord+=BoxY
	Output_gro.write(str(Unit_Counter).rjust(5)+moltype.ljust(5)+Label.ljust(5)+str(Atom_Counter).rjust(5)+f"{xcoord:.3f}".rjust(8)+f"{ycoord:.3f}".rjust(8)+f"{zcoord:.3f}".rjust(8)+"  0.0000  0.0000  0.0000\n")
	Output_xyz.write(Label[0]+" "+f"{10*xcoord:.3f}".rjust(8)+f"{10*ycoord:.3f}".rjust(8)+f"{10*zcoord:.3f}".rjust(8)+"\n")
	Atom_Counter += 1
	mol_atom_count += 1

	
## -------------- Dont edit anything above this line! --------------


## names of the ions to place in the structure
## note: the A cation and ligand are assumed to be organic molecules, and
## should already be oriented/moved correctly. In the .xyz file, 0,0,0 should
## be the location of whatever is supposed to be in the centre of the unit
## cell. In case of the ligand, the length of the tail is assumed to point in
## the positive z direction in the xyz file.
## The B cation and X anion are assumed to be atomic ions
## If no ligand is required (in case of 3D structures), anything can be
## entered, including an empty string.
Aion = "FA"  ## FA or MA or some atom
Bion = "Pd"
Xion = "I"
Ligand = "" ## PEA or nBA or ALA or PGA

## Parameters for special cases. Both should be True by default, but can be 
## swithed off. 
Aion_is_molecule = True ## to write the Aion in atom of molecule format
Crystal_2D = False ## whether to create a 2D layered structure (True), or 3D; without ligands (False)

NumXUnitCells = 4
NumYUnitCells = 4
NumZUnitCells = 3 ## complete unit cells per layer
Layers = 2 ## amount of layers

## what the files are called, and where to find them.
Molecule_folder = "../Optimized molecules/" ## example: "../Optimized molecules/"
Output_folder = "../gromacs input files/" ## example: "../gromacs input files/"
Molecule_handler = "_corr.xyz" ## example: "_corr.xyz"

CellXsize = 0.64  ## size in nm
CellYsize = 0.64
CellZsize = 0.64
Gap = 2-CellZsize ## A-A (or X-X) distance around ligand gap. When defined as
##					 [number]-CellZsize, [number] is the B-B distance.



## -------------- Dont edit anything below this line! --------------



## generate output name/files
if Crystal_2D:
	Output_name = Aion+Bion+Xion+"3-"+Ligand+"_"+str(Layers)+"_"+str(NumXUnitCells)+"-"+str(NumYUnitCells)+"-"+str(NumZUnitCells)
else:
	Output_name = Aion+Bion+Xion+"3-3D"+"_"+str(NumXUnitCells)+"-"+str(NumYUnitCells)+"-"+str(NumZUnitCells*Layers)
Output_gro = open(Output_folder+Output_name+".gro","w")
Output_xyz = open(Output_folder+Output_name+".xyz","w")

## determine box parameters
BoxX = NumXUnitCells * CellXsize
BoxY = NumYUnitCells * CellYsize
Layer_height = (NumZUnitCells+Crystal_2D)*CellZsize + Gap * Crystal_2D
Box_height = Layer_height * Layers

## print the properties of the chosen structure
print("\nChosen structure:\n   "+str(Layers)+" layer(s) of "+str(NumXUnitCells)+"*"+str(NumYUnitCells)+"*"+str(NumZUnitCells)+" unit cells")
if Crystal_2D: print("   "+Aion+Bion+Xion+"3 with "+Ligand+" as the A cation site ligand")
else: print("   "+Aion+Bion+Xion+"3 ")
print("   Size of single unit cell: "+str(CellXsize)+"*"+str(CellYsize)+"*"+str(CellZsize)+" nm")

## create ion_coords, ion_types, ion_atomcount, ion_labels
initialize(Aion)
if Crystal_2D: initialize(Ligand)

## calculations for amounts of atoms
NumUnitCells_stack = NumXUnitCells * NumYUnitCells * Layers
NumUnitCells = NumUnitCells_stack * NumZUnitCells

Aion_atomcount_total = NumUnitCells * Aion_atomcount
Bion_atomcount_total = NumUnitCells + NumUnitCells_stack * Crystal_2D
Xion_atomcount_total = NumUnitCells * 3 + NumUnitCells_stack * 4 * Crystal_2D
if Crystal_2D: Ligand_atomcount_total = NumUnitCells_stack * 2 * Ligand_atomcount * Crystal_2D
else: Ligand_atomcount_total = 0

Total_atomcount = Aion_atomcount_total + Bion_atomcount_total*2 + Xion_atomcount_total*2 + Ligand_atomcount_total
xyz_atomcount = Aion_atomcount_total + Bion_atomcount_total + Xion_atomcount_total + Ligand_atomcount_total

print("\nTotal amount of atoms: ", Total_atomcount)
print("   Atoms in "+Aion+": "+str(Aion_atomcount_total)+" in "+str(NumUnitCells)+" molecules (= "+str(Aion_atomcount)+" atoms each)")
print("   Atoms in "+Bion+": "+str(Bion_atomcount_total))
print("   Atoms in "+Xion+": "+str(Xion_atomcount_total))
if Crystal_2D: print("   Atoms in "+Ligand+": "+str(Ligand_atomcount_total)+" in "+str(NumUnitCells_stack*2)+" molecules (= "+str(Ligand_atomcount)+" atoms each)")

if Crystal_2D: Output_text = Aion+Bion+Xion+"3 with "+Ligand+" as ligand; "+str(Layers)+" layer(s) of "+str(NumXUnitCells)+"*"+str(NumYUnitCells)+"*"+str(NumZUnitCells)+" cells"
else: Output_text = Aion+Bion+Xion+"3; "+str(NumXUnitCells)+"*"+str(NumYUnitCells)+"*"+str(NumZUnitCells*Layers)+" cells"
Output_gro.write(Output_text+"\n"+str(Total_atomcount)+"\n")
Output_xyz.write(str(xyz_atomcount)+"\n"+Output_text+"\n")

Atom_Counter = 1
Unit_Counter = 1

for Layer in range(Layers): ## place all B cations
	for h in range(NumXUnitCells):
		for k in range(NumYUnitCells):
			for l in range(NumZUnitCells+Crystal_2D):
				xcoord = h*CellXsize
				ycoord = k*CellYsize
				zcoord = l*CellZsize + Layer * Layer_height
				write_atom(Bion)

for Layer in range(Layers): ## place all X anions
	for h in range(NumXUnitCells):
		for k in range(NumYUnitCells):
			for l in range(NumZUnitCells+Crystal_2D):
				xcoord = (h+0.5)*CellXsize
				ycoord = k*CellYsize
				zcoord = l*CellZsize + Layer * Layer_height
				write_atom(Xion)
				
				xcoord = h*CellXsize
				ycoord = (k+0.5)*CellYsize
				zcoord = l*CellZsize + Layer * Layer_height
				write_atom(Xion)
				
				xcoord = h*CellXsize
				ycoord = k*CellYsize
				zcoord = (l+0.5)*CellZsize + Layer * Layer_height
				write_atom(Xion)
				
				if l == 0 and Crystal_2D:
					xcoord = h*CellXsize
					ycoord = k*CellYsize
					zcoord = (l-0.5)*CellZsize + (Layer+1)*Layer_height
					write_atom(Xion)

for Layer in range(Layers): ## place all A cations
	for h in range(NumXUnitCells):
		for k in range(NumYUnitCells):
			for l in range(NumZUnitCells):
				if Aion_is_molecule is True:
					mol_atom_count=0
					for atom in range(Aion_atomcount):
						xcoord = (h+0.5)*CellXsize + Aion_coords[atom,0]
						ycoord = (k+0.5)*CellYsize + Aion_coords[atom,1]
						zcoord = (l+0.5)*CellZsize + Aion_coords[atom,2] + Layer * Layer_height
						Label = Aion_labels[mol_atom_count]
						write_molecule(Aion)
					Unit_Counter +=1
				else:
					xcoord = (h+0.5)*CellXsize
					ycoord = (k+0.5)*CellYsize
					zcoord = (l+0.5)*CellZsize + Layer * Layer_height
					write_atom(Aion)

if Crystal_2D: ## place all Ligands
	for Layer in range(Layers): 
		for h in range(NumXUnitCells):
			for k in range(NumYUnitCells):
				for up_down in range(2):
					mol_atom_count=0
					for atom in range(Ligand_atomcount):
						xcoord = (h+0.5)*CellXsize + Ligand_coords[atom,0]
						ycoord = (k+0.5)*CellYsize + Ligand_coords[atom,1]
						zcoord = (0.5-up_down)*CellZsize + NumZUnitCells*(1-up_down)*CellZsize + Ligand_coords[atom,2]*2*(0.5-up_down) + (Layer+up_down)*((NumZUnitCells+1)*CellZsize + Gap)
						Label = Ligand_labels[mol_atom_count]
						write_molecule(Ligand)
					Unit_Counter +=1

Output_xyz.close()
Output_gro.write(str(round(BoxX,2))+" "+str(round(BoxY,2))+" "+str(round(Box_height,2))+"\n")
Output_gro.close()
