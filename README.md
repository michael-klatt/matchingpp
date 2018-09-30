# matchingpp
**Stable Matching of Point Patterns**  
An R package for a mutual nearest neighbor matching of two point
patterns on the torus in arbitrary dimension

Author: Michael Andreas Klatt  
Maintainer: Michael Andreas Klatt (<software@mklatt.org>)  
License: GNU GPLv3  
Version: 0.1.0

Two point patterns (e.g. a stationarized lattice, a Poisson point
process, or data from csv files) are matched to each other via the
stable mutual nearest neighbor matching. Typically this results in a
thinnig of one of the two point patterns.

MatchingPoissonToLattice
------------------------
Matches a stationarized lattice to a Poisson point process
More details are avaiable via the command
	help(MatchingPoissonToLattice)

MatchingInfileToLattice
-----------------------
Matches a stationarized lattice to a point pattern whose points are stored in a csv file, where each row stores the coordinates of one point.
More details are avaiable via the command
	help(MatchingInfileToLattice)

MatchingInfileToInfile
----------------------
Matches two point patterns whose points are stored in two separate csv files, where each row stores the coordinates of one point.
More details are avaiable via the command
	help(MatchingInfileToInfile)

MatchingPP
----------
Core function to match two point patterns.
More details are avaiable via the command
	help(MatchingPP)

Installation
------------
                       Install from GitHub :  devtools::install_github("michael-klatt/matchingpp")
This package depends on the package RANN.
