noe-way-jose
============

Probabilistic NOE module for PHAISTOS
 
 
 - Options:


    energy-noe = 1                                                  # Activate energy-noe [number of occurrences]
    energy-noe-debug = 0                                            # Debug level
    energy-noe-weight = 1                                           # Weight used when summing energy terms
    energy-noe-weight-constant = 1                                  # Apply an extra weight factor (default is 1.0)
    energy-noe-active-restraints = 1                                # Number of active restraints (default is 1)
    energy-noe-seamless = false                                     # Seamless restraint switching -- experimental! (default is False)
    energy-noe-upl-filename = ""             # CYANA UPL formatted list of NOE contacts

 - Input format is the CYANA UPL format, i.e.:


    178 GLU  H     179 GLU  H       4.52
    179 GLU  H     180 LYS  H       4.66
    28  TYR  H      28 TYR  QD      3.613
    28  TYR  H     176 ILE  QD1     5.532
    27  VAL  QG1    28 TYR  H       4.128
    27  VAL  QG2    28 TYR  H       4.12
    28  TYR  H     171 LEU  H       3.829


 - Only N out of all restraints are active at a time.

 - The potential is the ROSETTA bounded potential explained here:

    https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/de/d50/constraint_file.html
    http://www.pnas.org/content/suppl/2012/06/25/1203013109.DCSupplemental/Appendix.pdf

 - Force constant simpy multiplies the potential by a scaling factor.

 - Seamless restraint switching adds a bias, which connects potential surfaces when switching from one set of restraints to another. This effectively forces switching of restraints. Otherwise, they are sampled. This can be used to escape local minima in simulations.

 - When used as an observable, all NOEs are used, regardless of the specified number of active restraints. Additionally, then number of restraints violations are reported The output format is [XXXX.XX,YY,ZZ] where XXXX.XX is the energy, YY is the number of NOE restraint violations less than 1 angstrom, and ZZ is the number of NOE restraint violations greater than 1 angstrom.
