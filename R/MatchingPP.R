#' Stable matching of point patterns
#'
#' Mutual nearest neighbor matching of two point patterns on the torus.

#' imports RANN::n2

# STABLE MATCHING: USER FUNCTIONS ---------------------------------------------
#' Stable matching of a randomized lattice to a Poisson point process
#'
#' @name MatchingPoissonToLattice
#' @usage MatchingPoissonToLattice(dimension, linear.system.size, intensity, fileout=FALSE, filename="example-matchedpp.dat", unthinned.fileout=FALSE, unthinned.filename="example-ppp.dat", verbose=TRUE)
#' @description Stable matching of a randomized lattice to a Poisson point process
#'
#' @param dimension Dimension d of the simulation box
#' @param linear.system.size Linear system size L of the simulation box [0,L)^d
#' @param intensity Intensity of the Poisson point process (mean number of points in[0,1)^d)
#' @param fileout Output of matched points to file if set to TRUE (default=FALSE)
#' @param filename Name of output file (default="example-matchedpp.dat")
#' @param unthinned.fileout Output of unthinned Poisson point process to file if set to TRUE (default=FALSE)
#' @param unthinned.filename Name of output file for the unthinned Poisson point process (default="example-ppp.dat")
#' @param verbose Output of progress if set to TRUE (default=TRUE)
#' @return Matrix that stores all matched points
#'         of the Poisson point process (in columns 1 to dimension)
#'         and of the lattice (in columns dimension+1 to end).
#'         Each row contains one pair of matched points.
#' @examples
#' out.1 <- MatchingPoissonToLattice(dimension=3,linear.system.size=4.0,intensity=25.0)
#' out.2 <- MatchingPoissonToLattice(2,10.0,2.5,fileout=TRUE,filename="example-2D-matched-possion-to-Z.dat")
#' @export
MatchingPoissonToLattice <- function(dimension, linear.system.size, intensity,
                                     fileout=FALSE,
				     filename="example-matchedpp.dat",
                                     unthinned.fileout=FALSE,
                                     unthinned.filename="example-ppp.dat",
                                     verbose=TRUE) {
  # Screen message ------------------------------------------------------------
  if(verbose){
    cat(paste0("Matching a Poisson point process to Z+U\n"))
  }

  # Renaming variables for convenience ----------------------------------------
  dim <- dimension
  L <- linear.system.size

  # Poisson point process -----------------------------------------------------
  p1.N <- rpois(1,intensity*L^dim)
  p1.coords <- runif(dim*p1.N,0,L)
  p1.coords <- matrix(p1.coords, nrow = p1.N)
  if(unthinned.fileout) {
    write.table(p1.coords, file=unthinned.filename,
                row.names=FALSE, col.names=FALSE)
  }

  # Randomized lattice --------------------------------------------------------
  # Unit intensity
  p2.l <- 1              # Lattice constant (linear size of unit cell)
  # Random shift of lattice
  p2.offset <- runif(dim,0,p2.l)
  # Number of lattice points
  p2.lN <- ceiling(L-p2.offset) # Lattice: linear number of lattice points
  p2.N <- prod(p2.lN)           # Total number of points
  # Sample shifted lattice
  p2.coords <- Lattice(dim,p2.lN,p2.l)
  p2.coords <- sweep(p2.coords,2,-p2.offset)

  # MATCHING THE TWO POINT PATTERNS -------------------------------------------
  matched.pp <- MatchingPP(L, p1.coords, p2.coords,
                           fout = fileout,
                           fname = filename,
                           verb = verbose)

  return(matched.pp)
}

#' Stable matching of points from an infile to a randomized lattice
#'
#' @name MatchingInfileToLattice
#' @usage MatchingInfileToLattice(dimension, linear.system.size, infile, fileout=FALSE, filename="example-matchedpp.dat", verbose=TRUE)
#' @description Stable matching of points from an infile to a randomized lattice
#'
#' @param dimension Dimension d of the simulation box
#' @param linear.system.size Linear system size L of the simulation box [0,L)^d
#' @param infile Name of file that contains a matrix, where each row stores the coordinates of one point
#' @param fileout Output of matched points to file if set to TRUE (default=FALSE)
#' @param filename Name of output file (default="example-matchedpp.dat")
#' @param verbose Output of progress if set to TRUE (default=TRUE)
#' @return Matrix that stores all matched points
#'         of the point pattern from infile (in columns 1 to dimension)
#'         and of the lattice (in columns dimension+1 to end).
#'         Each row contains one pair of matched points.
#' @examples
#' tmp.infile <- "example-ppp.dat"
#' linear.sys.size <- 4.0
#' n.pts <- rpois(1,1000.0*linear.sys.size^2)
#' coords.pois <- matrix(runif(2*n.pts,0,linear.sys.size), nrow = n.pts)
#' write.table(coords.pois,file=tmp.infile,row.names=FALSE, col.names=FALSE)
#' out <- MatchingInfileToLattice(dimension=2, linear.sys.size, tmp.infile)
#' @export
MatchingInfileToLattice <- function(dimension, linear.system.size, infile,
                                    fileout=FALSE,
				    filename="example-matchedpp.dat",
                                    verbose=TRUE) {
  # Screen message ------------------------------------------------------------
  if(verbose){
    cat(paste0("Matching ",infile," to Z+U\n"))
  }

  # Renaming variables for convenience ----------------------------------------
  dim <- dimension
  L <- linear.system.size

  # Read in data and extract parameters ---------------------------------------
  # Coordinates in infile
  p1.coords <- read.csv(file=infile, header=FALSE, sep=' ')
  # Number of points in infile
  p1.N <- dim(p1.coords)[1]
  # Check dimensions
  stopifnot(dim == dim(p1.coords)[2])

  # Randomized lattice --------------------------------------------------------
  # Unit intensity
  p2.l <- 1              # Lattice constant (linear size of unit cell)
  # Random shift of lattice
  p2.offset <- runif(dim,0,p2.l)
  # Number of lattice points
  p2.lN <- ceiling(L-p2.offset) # Lattice: linear number of lattice points
  p2.N <- prod(p2.lN)           # Total number of points
  # Sample shifted lattice
  p2.coords <- Lattice(dim,p2.lN,p2.l)
  p2.coords <- sweep(p2.coords,2,-p2.offset)

  # MATCHING THE TWO POINT PATTERNS -------------------------------------------
  matched.pp <- MatchingPP(L, p1.coords, p2.coords,
                           fout = fileout,
                           fname = filename,
                           verb = verbose)

  return(matched.pp)
}

#' Stable matching of point patterns from two infiles
#'
#' @name MatchingInfileToInfile
#' @usage MatchingInfileToInfile(linear.system.size, infile.1, infile.2, fileout=FALSE, filename="example-matchedpp.dat", verbose=TRUE)
#' @description Stable matching of point patterns from two infiles
#'
#' @param linear.system.size Linear system size L of the simulation box [0,L)^d
#' @param infile.1 Name of a file that contains a matrix, where each row stores the coordinates of one point
#' @param infile.2 Name of another file that contains a matrix, where each row stores the coordinates of one point
#' @param fileout Output of matched points to file if set to TRUE (default=FALSE)
#' @param filename Name of output file (default="example-matchedpp.dat")
#' @param verbose Output of progress if set to TRUE (default=TRUE)
#' @return Matrix that stores all matched points
#'         of the point pattern from infile.1 (in columns 1 to dimension)
#'         and of the point pattern from infile.2 (in columns dimension+1 to end).
#'         Each row contains one pair of matched points.
#' @examples
#' tmp.bino <- "example-bpp.dat"
#' tmp.pois <- "example-ppp.dat"
#' linear.sys.size <- 4.0
#' n.pts <- 25
#' coords.bino <- matrix(runif(2*n.pts,0,linear.sys.size), nrow = n.pts)
#' n.pts <- rpois(1,n.pts)
#' coords.pois <- matrix(runif(2*n.pts,0,linear.sys.size), nrow = n.pts)
#' write.table(coords.bino,file=tmp.bino,row.names=FALSE, col.names=FALSE)
#' write.table(coords.pois,file=tmp.pois,row.names=FALSE, col.names=FALSE)
#' out <- MatchingInfileToInfile(linear.sys.size, tmp.pois, tmp.bino)
#' @export
MatchingInfileToInfile <- function(linear.system.size, infile.1, infile.2,
                                   fileout=FALSE,
				   filename="example-matchedpp.dat",
                                   verbose=TRUE) {
  # Screen message ------------------------------------------------------------
  if(verbose){
    cat(paste0("Matching ",infile.1," to ",infile.2,"\n"))
  }

  # Renaming variable for convenience -----------------------------------------
  L <- linear.system.size

  # Read in data and extract parameters ---------------------------------------
  # Coordinates in infile.1
  p1.coords <- read.csv(file=infile.1, header=FALSE, sep=' ')
  # Number of points in infile.1
  p1.N <- dim(p1.coords)[1]
  # Coordinates in infile.2
  p2.coords <- read.csv(file=infile.2, header=FALSE, sep=' ')
  # Number of points in infile.2
  p2.N <- dim(p2.coords)[1]
  # Check dimensions
  stopifnot(dim(p1.coords)[2] == dim(p2.coords)[2])
  dim <- dim(p1.coords)[2]

  # MATCHING THE TWO POINT PATTERNS -------------------------------------------
  matched.pp <- MatchingPP(L, p1.coords, p2.coords,
                           fout = fileout,
                           fname = filename,
                           verb = verbose)

  return(matched.pp)
}

# STABLE MATCHING: MAIN FUNCTION ----------------------------------------------
#' Stable matching of point patterns
#'
#' @name MatchingPP
#' @usage MatchingPP(linear.system.size, coords.1, coords.2, fout=FALSE, fname="example-matchedpp.dat", verb=TRUE, nn2.treetype = "kd", nn2.eps=0)
#' @description Stable matching of two point patterns in a cubic simulation box [0,L)^d:
#'              It is the core function of the package MatchingPP.
#'
#'              The nearest neighbor search uses the function nn2 from the package RANN.
#'
#' @param linear.system.size Linear system size L of the simulation box [0,L)^d
#' @param coords.1 A point pattern stored in a matrix, where each row stores the coordinates of one point
#' @param coords.2 Another point pattern stored in a matrix, where each row stores the coordinates of one point
#' @param fout Output of matched points to file if set to TRUE (default=FALSE)
#' @param fname Name of output file (default="example-matchedpp.dat")
#' @param verb Output of progress if set to TRUE (default=TRUE)
#' @param nn2.treetype Optional argument of the function nn2 from the package RANN, which is used for the nearest neighbor search: "Either the standard kd tree or a bd (box-decomposition, AMNSW98) tree which may perform better for larger point sets" (default="kd")
#' @param nn2.eps Optional argument of the function nn2 from the package RANN, which is used for the nearest neighbor search: "error bound: default of 0.0 implies exact nearest neighbour search" (default=0)
#' @return Matrix that stores all matched points
#'         of the 1st point pattern (in columns 1 to dimension)
#'         and of the 2nd point pattern (in columns dimension+1 to end).
#'         Each row contains one pair of matched points.
#' @examples
#' linear.sys.size <- 4.0
#' n.pts <- 25
#' coords.bino <- matrix(runif(2*n.pts,0,linear.sys.size), nrow = n.pts)
#' n.pts <- rpois(1,n.pts)
#' coords.pois <- matrix(runif(2*n.pts,0,linear.sys.size), nrow = n.pts)
#' out <- MatchingPP(linear.sys.size, coords.pois, coords.bino)
#' @export
MatchingPP <- function(linear.system.size, coords.1, coords.2,
                       fout=FALSE, fname="example-matchedpp.dat", verb=TRUE,
                       nn2.treetype = "kd", nn2.eps=0) {
  # Check dimensions
  stopifnot(dim(coords.1)[2] == dim(coords.2)[2])
  dimension <- dim(coords.1)[2]

  matched.all <- matrix(0, nrow = 0, ncol = 2*dimension)

  # Check which point patterns will contain unmatched points
  # Number of points in point patterns
  N1 <- nrow(coords.1)
  N2 <- nrow(coords.2)
  swap_pp <- FALSE
  if(N2 < N1){
    swap_pp <- TRUE
    tmp_coords <- coords.1
    coords.1 <- coords.2
    coords.2 <- tmp_coords
    N1 <- nrow(coords.1)
    N2 <- nrow(coords.2)
  }

  # Match mutual nearest neighbors until all points are matched
  repeat{
    if(verb) {
      cat(paste0("Remaining points to be matched: ", nrow(coords.1),"\n"))
    }

    # Number of points in point patterns
    N1 <- nrow(coords.1)
    N2 <- nrow(coords.2)

    # Prepare periodic copies for neighbor search on torus
    p1.pbc <- PeriodicCopies(dimension,coords.1,linear.system.size)
    p2.pbc <- PeriodicCopies(dimension,coords.2,linear.system.size)

    # Nearest neighbor search (parameters)
    # treetype = "kd" or "bd"
    # eps = error bound (accepting deviations in nearest neighbor search)
    # searchtype = c("standard","priority","radius"),radius=0.0
    nn2.searchtype <- "priority"
    nn.of.p2.in.p1 <- nn2(p1.pbc, query = coords.2, k = 1,
                          treetype=nn2.treetype, searchtype=nn2.searchtype, eps=nn2.eps)
    nn.of.p1.in.p2 <- nn2(p2.pbc, query = coords.1, k = 1,
                          treetype=nn2.treetype, searchtype=nn2.searchtype, eps=nn2.eps)

    # Find mutual nearest neighbors -- representative in [0,L)^d
    of.p1 <- apply(nn.of.p1.in.p2$nn.idx,1,function(x) (x-1)%%N2 + 1)
    of.p2 <- apply(nn.of.p2.in.p1$nn.idx,1,function(x) (x-1)%%N1 + 1)

    idx <- 1:N1
    matched.in.p1 <- idx[of.p2[of.p1] == 1:N1]
    matched.in.p2 <- of.p1[matched.in.p1]

    # save matched points
    if(swap_pp){
      matched <- cbind(coords.2[matched.in.p2,, drop=FALSE],
  		                 coords.1[matched.in.p1,, drop=FALSE])
      matched.all <- rbind(matched.all,matched)
    }
    else{
      matched <- cbind(coords.1[matched.in.p1,, drop=FALSE],
  		                 coords.2[matched.in.p2,, drop=FALSE])
      matched.all <- rbind(matched.all,matched)
    }

    if(nrow(coords.1) == nrow(matched)) {
      if(verb){cat("Remaining points to be matched: 0\n")}
      break;
    }

    # Remove the matched points from the point sets
    coords.1 <- coords.1[-matched.in.p1,, drop=FALSE]
    coords.2 <- coords.2[-matched.in.p2,, drop=FALSE]
  }

  # Output of pairs of matched points
  if(fout) {
    write.table(matched.all, file = fname,
                row.names = FALSE, col.names = FALSE)
  }

  colnames(matched.all) <- NULL
  return(as.matrix(matched.all))
}

# AUXILIARY FUNCTIONS ---------------------------------------------------------
#' Coordinates of a lattice
#'
#' @name Lattice
#' @usage Lattice(dimension,lin.N,lin.len)
#' @description Creates a patch of the lattice a*Z^d, where d is the dimension (auxiliary function)
#'
#' @param dimension Dimension of the simulation box
#' @param lin.N Vector of linear numbers of lattice points
#' @param lin.len Lattice constant "a" (nearest neighbor distance)
#' @return Matrix R^{prod(lin.N) x dimension}, each row contains coordinates of one cell center
Lattice <- function(dimension,lin.N,lin.len) {

  coords <- c()
  i=1
  while(i<=dimension){
    coords <- c(coords, list((0:(lin.N[i]-1))*lin.len))
    i <- i+1
  }
  return(as.matrix(expand.grid(coords)))
  # Inspired by stackoverflow.com/questions/9422945
}

#' Create periodic copies of points in a simulation box
#'
#' @name PeriodicCopies
#' @usage PeriodicCopies(dimension,coordinates,L)
#' @description Create 3^d periodic copies of points in [0,L)^d (auxiliary function)
#'
#' @param dimension Dimension of the simulation box
#' @param coordinates Matrix that stores the coordinates of the points
#' @param L Linear size of cubic simulation box
#' @return Matrix that stores all periodic copies
PeriodicCopies <- function(dimension,coordinates,L) {

  pbc.lattice <- Lattice(dimension,rep(3,dimension),L) - L

  p1.pbc.list <- list()
  for (i in 1:3^dimension){
    p1.pbc.list[[i]] <- sweep(coordinates,2,-pbc.lattice[i,])
  }
  return(do.call(rbind,p1.pbc.list))
}

#' Distances between matched points on a torus
#'
#' @name DistancesOnTorus
#' @usage DistancesOnTorus(dimension, linear.system.size, coords.1, coords.2)
#' @description Output of distances on torus between points in two patterns (auxiliary function)
#'
#' @param dimension Dimension of the simulation box
#' @param linear.system.size Linear system size L of the simulation box
#' @param coords.1 Matrix with coordinates of points in 1st pattern (each row stores the coordinates of one point)
#' @param coords.2 Matrix with coordinates of points in 2nd pattern (each row stores the coordinates of one point)
#' @return Vector of distances, each row corresponds to the distances of the
#'  corresponding rows in coords.1 and coords.2
DistancesOnTorus <- function(dimension, linear.system.size, coords.1, coords.2){

  abs.diff <- abs(coords.1-coords.2)
  dist.pbc <- apply(pmin(abs.diff,linear.system.size-abs.diff), 1,
                    function(x) sqrt(sum(x^2)))
  return(dist.pbc)
}
