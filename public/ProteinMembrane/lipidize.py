import MDAnalysis
import numpy as np
from random import shuffle
import re
import sys
import os.path as path

def norm(x):
    return x / np.sqrt(x.dot(x))

def flatten(x): 
    return [item for sublist in x for item in sublist]

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def pdbRemark(filename):
    remark= []; lipids = []; fract = []
    file = open(filename, "r")
    for line in file:
        if re.search('REMARK 547', line):
              remark.append(line.strip().split(' ')[2:])

    for line in remark:
        for l,f in zip(line[0::2], line[1::2]):
            lipids.append(l.lower())
            fract.append(num(f))

    file.close()

    return lipids,fract

def replaceDum(dummypdb,nconformer=50,thickscale=0.992):
    
    # Structure with dummy atoms
    uprot = MDAnalysis.Universe(dummypdb)
   
    # Parse the input pdb for a remark 547 (LIP) line
    lipidpdb, composition = pdbRemark(dummypdb)

    if lipidpdb == [] or composition == []\
	    or (len(lipidpdb) != len(composition)):
        raise compositionException(
		"Malformed or no REMARK 547 found,"\
		"lipids or compositions are not specified") 

    # Load in the lipids we need
    ulipid = []
    for lip in lipidpdb:
        ulipid.append(MDAnalysis.Universe('./pdbs/'+lip+'.pdb.gz'))
    
    # Selection for specified lipids
    sel_lipid = [];
    nlipid = len(ulipid)
    for u in ulipid:
        sel = u.select_atoms('all')
        sel_lipid.append(sel)
        
    # Create selections for dummy atoms, select dummys
    # within a distance of the protein if necessary
    dum_n = uprot.select_atoms('resname DUM and name N and (not around 5 protein)')
    dum_o = uprot.select_atoms('resname DUM and name O and (not around 5 protein)')
    dum_no = dum_n + dum_o
  
    # Number of dummy atoms
    ndum =  dum_no.n_atoms
    ndumn = dum_n.n_atoms
    ndumo = dum_o.n_atoms
        
    # We need to specify some form of composition / ratio
    # of the lipids. By Default, assume a balanced amount
    # of the specified lipids
    # TODO: Make robust when ndum * species_frac != integer
    composition = np.array(composition);
    if not composition.size:
        composition = (np.array([1.]*nlipid)) / nlipid
    composition *= ndum

    # Create a list of lipids to be merged
    tomerge = []
    for l,c in zip(sel_lipid,composition):
        tomerge.append([l]*np.int(c))
        
    # Flatten the list and shuffle so placement on the
    # surface is random
    tomerge = flatten(tomerge)
    shuffle(tomerge)

    # Select all the atoms for lipid and create a new
    # universe that has ndum lipid molecules
    mergemol = MDAnalysis.Merge(*tomerge)
    sel = mergemol.select_atoms('all')
    nlipidatoms = sel.n_atoms
    
    # Cleanup the lipid numbering, segments etc..
    resid = 0; resids = []
    for lipid in tomerge:
        resids.append([resid]*lipid.n_atoms)
        resid += 1
    sel.set_resids(flatten(resids))
            
    # Array to store coordinates of lipid molecules
    newcoords = np.empty([nlipidatoms,3])
    start = 0; end = tomerge[0].n_atoms
      
    # Generate the new popc coordinates for each 'O' dummy atom
    np.random.shuffle(dum_o.positions)
    for pos, lipid in zip(dum_o.positions,tomerge[0:ndumo]):
        
        # Select a random conformer from the lipid database
        tsid = np.random.randint(low=0, high=nconformer)
        lipid.universe.trajectory[tsid]
        
        # Store the lipid position
        init_lipid_pos = lipid.get_positions()
        
        # Align the principal axis with the dummy atom vector
        lipid.align_principal_axis(0,norm(-1*pos))
       
        # Rotate about the principal axes some random amount to
        # give the illusion (MAGIC!) of disorder
        e1, e2, e3 = lipid.principal_axes()
        
        lipid.rotateby(np.random.uniform(low=-180, high=180,size=1),
                       norm(pos),[0, 0, 0])
        lipid.rotateby(np.random.uniform(low=-20, high=20,size=1),
                       norm(e2),[0, 0, 0])
        lipid.rotateby(np.random.uniform(low=-20, high=20,size=1),
                       norm(e3),[0, 0, 0])
        
        # Translate to the dummy atom position, but add some fuzzyness
        end = start + lipid.n_atoms
        newcoords[start:end] = pos * thickscale + lipid.get_positions() +\
                np.random.uniform(low=-2.0,high=2.0,size=[1,3])
        
	# Coordinate array offsets
        start += lipid.n_atoms
               
        # Return the lipid to its initial position before moving
        lipid.positions = init_lipid_pos
        
    # Generate the new popc coordinates for each 'N' dummy atom
    np.random.shuffle(dum_n.positions)
    for pos, lipid in zip(dum_n.positions,tomerge[ndumo:ndum]):
        tsid = np.random.randint(low=0, high=nconformer)
        lipid.universe.trajectory[tsid]
        init_lipid_pos = lipid.get_positions()
        lipid.align_principal_axis(0,norm(pos))
        e1, e2, e3 = lipid.principal_axes()
        lipid.rotateby(np.random.uniform(low=-180, high=180, size=1),
                       norm(pos),[0, 0, 0])
        lipid.rotateby(np.random.uniform(low=-20, high=20,size=1),
                       norm(e2),[0, 0, 0])
        lipid.rotateby(np.random.uniform(low=-20, high=20,size=1),
                       norm(e3),[0, 0, 0])
        end = start + lipid.n_atoms
        newcoords[start:end] = pos * (1.0 / thickscale) + lipid.get_positions() +\
                np.random.uniform(low=-2.0,high=2.0,size=[1,3])
        start += lipid.n_atoms
        lipid.positions = init_lipid_pos
        
    # Apply the coordinates to the selection
    sel.set_positions(newcoords)
   
    ## Merge the protein and the lipds
    prot   =  uprot.select_atoms('protein')
    lipids =  mergemol.select_atoms('all') 
    lipprot = MDAnalysis.Merge(prot,lipids)

    # Write out the pdb
    lipprot.atoms.write(path.splitext(dummypdb)[0]+'lipids.pdb')


class compositionException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

if __name__ == "__main__":
	replaceDum(sys.argv[1])

