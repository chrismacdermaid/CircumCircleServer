import MDAnalysis
import numpy.linalg

u = MDAnalysis.Universe('cer160.pdb.gz')
sel = u.select_atoms("all")

print sel.n_atoms
