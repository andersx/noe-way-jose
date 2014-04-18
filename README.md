noe-way-jose
============

Probabilistic NOE module for PHAISTOS
 
    --energy-noe [=arg(= )]                 Activate energy-noe [number of occurrences]
    --energy-noe-debug arg (=0)             Debug level
    --energy-noe-weight arg (=1)            Weight used when summing energy terms
    --energy-noe-force-constant arg (=1)    The harmonic force constant (default is 1.0)
    --energy-noe-active-restraints arg (=1) Number of active restraints (default is 1)
    --energy-noe-seamless arg (=0)          Seamless restraint switching -- experimental! (default is False)
    --energy-noe-contact-map-filename arg   List of NOE contacts.


 - The potential is the ROSETTA bounded potential explained here:

    https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/de/d50/constraint_file.html
    http://www.pnas.org/content/suppl/2012/06/25/1203013109.DCSupplemental/Appendix.pdf

 - Force constant simpy multiplies the potential by a scaling factor.

 - Seamless restraint switching adds a bias, which connects potential surfaces when switching from one set of restraints to another. This effectively forces switching of restraints.

 - Input file format is the same as the energy-contact-map. E.g.

    19 CA 56 CA 7.0

 - E.g. "resid type resid type distance", where  resid, runs from 1 to N.
