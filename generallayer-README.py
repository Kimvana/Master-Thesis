

## instead of 'def initialize(Ligand)', the following lines of text could 
## be used. The variables 'Molecule_folder', 'Ligand' and 'Molecule handler'
## are all defined in the main file before calling this function.

## the function loads in the coordinates and atom types from an external
## .xyz file. This file contains a single (vacuum geometry optimized) molecule.
## As mentioned in the main file, these molecules are assumed to be oriented
## and translated correctly. The 0,0,0 point in the .xyz file is taken as the
## part that goes into the centre of the unit cell. The orientation in the 
## .xyz file is the orientation used in the generated structure. Especially
## for the A2 cations (ligand for clarity in code), it is important to keep
## in mind that the protruding tail should be positioned in the positive
## z-direction in these .xyz files.


Ligand_coords=np.genfromtxt(Molecule_folder+Ligand+Molecule_handler, skip_header=2, usecols=(1,2,3))
Ligand_types=np.genfromtxt(Molecule_folder+Ligand+Molecule_handler, skip_header=2, usecols=0, dtype=str)
Ligand_coords.flags.writeable = False
Ligand_types.flags.writeable = False
Ligand_atomcount=Ligand_types.shape[0]

Ligand_labeldict={}
Ligand_labels=[]
for atom in Ligand_types:
	if atom in Ligand_labeldict: Ligand_labeldict[atom]+=1
	else: Ligand_labeldict[atom]=1
	Ligand_labels.append(atom+str(Ligand_labeldict[atom]))
	
	
## The code contains 2 options for the Aion, and 4 options for the ligand.
## In principle, the Aion options could also be inserted as ligand, but this
## does not make any physical sense. 
## These lists can be appended relatively easily: One only has to create the
## .xyz files to the aforementioned specifications, and add them to the file
## location. Of course, the code makes Gromacs input structures. Gromacs also
## needs .itp files, which have to be created separately.


## Special case parameters. In case one wants to make structures with an
## atomic Aion (like CsSnI3), Aion_is_molecule must be switched to false, and
## the desired atomtype should be defined in Aion = [atomtype]. Crystal_2D 
## generates layered structures when set to True, and 'bulk' 3D materials 
## when set to False.
## In case of Crystal_2D, the value is used as a number: True functions like
## 1 in calculations, False like 0.
Aion_is_molecule = True
Crystal_2D = True


## these are to define file locations. In the example, it looks for the 
## following location when promted to look for FA:
## .../Optimized molecules/FA_corr.xyz
Molecule_folder = ""


## generate output name/files.

## The names are designed such that the name contains all required info.
## Example names: FAPbI3-PEA_2_4-4-3, FASnI3-3D_4-4-4


## calculations for amounts of atoms

## NumUnitCells is the total amount of complete unit cells. In other words,
## it is the amount of Aion molecules present in the structure.
## Considering the 2D formula (Ligand)_2(Aion)_n(Bion)_(n+1)(Xion)_(3n+4),
## NumUnitCells is the factor n, NumUnitCells_stack is used for the additional
## 1 and 4.


## Place all B cations

## all are placed at 0,0,0 of every unit cell. 
## if Crystal_2D = True, for every layer, an extra set is required. To place
## these, l loops over 'NumZUnitCells + Crystal_2D'


## Place all X anions

## all are placed at (0.5, 0 ,0), (0, 0.5, 0) and (0, 0, 0.5).
## if Crystal_2D = True, for every Layer, 4 extra's are required. Three of
## these are already placed using the 'standard' loops (that is why l loops
## over 'NumZUnitCells + Crystal_2D'), the fourth is placed with the extra
## section.


## Place all A cations

## all are placed at (0.5, 0.5, 0.5)


## Place all Ligands
 
## Very similar to place A cations. However, it does not loop over l, as every
## layer only needs a pair of ligands. As it is a pair, it loops over up_down.
## This up_down factor is also used to translate (and mirror) the ligands in
## the correct way to build into the structure.
## This mirroring should be noted! If the Ligand is chiral, one half of the
## interface gets one chirality, the other half the other. Should be fixable
## by also mirroring the x OR y coordinate, but was not considered relevant
## during this project.




