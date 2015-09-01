"""
-----------
Setup stage
-----------

*** Identify residues in the molecule - playground/rotamer/identify_residue.py
*** Find dihedrals and store atom identities. - playground/rotamer/find_dihedrals.py
* Choose backbone-dependent, local environment etc. to determine what sort of library data to generate.
* Read in relevant part of rotamer library data: store dictionary for each residue's rotamer states.

----------
Move stage
----------

*** Read in coordinates from GMIN - playground/amber/coords_io.py
* Select a set of residues to move.
* Measure rotamer state of residues and their neighbours (if appropriate).
* Select a new conformation from the rotamer library, or a random configuration, with probability given by
  Good-Turing frequency estimation.
*** Write coordinates for GMIN - playground/amber/coords_io.py

------
Future
------

- Decomposed energies for dihedrals wrt other sidechain dihedrals.

"""

