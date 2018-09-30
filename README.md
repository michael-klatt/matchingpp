# matchingpp
**An R Package for the Stable Matching of Point Patterns**  
Author: Michael Andreas Klatt (<software@mklatt.org>)  
License: GNU GPLv3  
Version: 0.1.0  

Two point patterns (e.g. a stationarized lattice, a Poisson point
process, or data from csv files) are matched to each other via the
stable mutual nearest neighbor matching (on the torus in arbitrary
dimensions). Typically this results in a thinnig of one of the two point
patterns.

**Based on the paper**  
Michael Andreas Klatt, Günter Last and D. Yogeshwaran.
*Hyperuniform and rigid stable matchings* (2018)

Installation
------------
To install the package from github:    

                       devtools::install_github("michael-klatt/matchingpp")
To test the package:  

                       library('matchingpp')
                       example('MatchingPoissonToLattice')
Note that this package depends on the package RANN.

R Functions
-----------
**MatchingPoissonToLattice**
matches a stationarized lattice to a Poisson point process.  
More details are avaiable via the command  

                       help(MatchingPoissonToLattice)

**MatchingInfileToLattice**
matches a stationarized lattice to a point pattern whose points are
stored in a csv file, where each row stores the coordinates of one
point.  
More details are avaiable via the command  

                       help(MatchingInfileToLattice)

**MatchingInfileToInfile**
matches two point patterns whose points are stored in two separate csv
files, where each row stores the coordinates of one point.  
More details are avaiable via the command  

                       help(MatchingInfileToInfile)

**MatchingPP**
is the core function to match two point patterns.  
More details are avaiable via the command  

                       help(MatchingPP)

