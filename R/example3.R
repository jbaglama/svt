# PAPER:
#  J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm for 
#  Computing a Partial Singular Value Decomposition Satisfying a Given 
#  Threshold", submitted for publication 2024.
# 
# EXAMPLE 3 FROM PAPER:
# Demonstrates how the energy percentages for the image 'tiger' from the R- 
# package rsvd are computed and obtained. The rsvd reference:
#                              Erichson, N. B., Voronin, S., Brunton, 
#                              S. L., & Kutz, J. N., Randomized Matrix 
#                              Decompositions Using R. Journal of Statistical 
#                              Software, 89, 1-48.  
#              
# TO RUN DEMO IN R:
# set the working directory - see below   
#  > source("example3.R")
#  > example3()
#
# REQUIRED SOFTWARE:  
# svt_irlba.R
# https://github.com/jbaglama/svt/
#
# DATE LAST MODIFIED: 
# 7/8/24
#
# LANGUAGE:
# R version: 4.3.2 
# 
# AUTHORS: 
# James Baglama            email: jbaglama@uri.edu
# Jonathan Chávez-Casillas email: jchavezc@uri.edu
# Vasilije Perovic         email: perovic@uri.edu
# ---------------------------------------------------
example3 <- 
  function(){
    
# ** IMPORTANT **  
# Working directory where script and files are located
# ** NEEDS TO BE SET BY USER ** for example
# > workingdir <- getwd()
# > setwd(workingdir)
# Alternatively, the user will have to run the 
# 'example43()' function into the R-environment.
# -----------------------------------------------------

# Needed Package names
# ---------------------
packages <- c("irlba","Matrix","methods","R.matlab","SparseM","rsvd")

# Installing packages not yet installed.
# --------------------------------------
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
# -----------------
invisible(lapply(packages, library, character.only = TRUE))

# load svt_irlba script
# ----------------------
source("svt_irlba.R")

# load the tiger data set from package rsvd- it is loaded directly into a matrix. 
# tiger is a 1600 x 1200 matrix.
# ------------------------------------------------------------------------------
data("tiger", package = "rsvd")
m <- nrow(tiger)
n <- ncol(tiger)

# Computing the Frobenius norm of tiger matrix.
# ---------------------------------------------
nrmse <- norm(tiger,"F")^2

# Setting tolerance and psvdmax.
# ------------------------------
tol <- 1e-5
psvdmax <- min(m,n)

# 1. Call the program svt_irlba with energy = 0.9854
# --------------------------------------------------
energy1 <- 0.9854
tstart <- proc.time()
psvd <- svt_irlba(tiger,tol = tol, energy = energy1,psvdmax = psvdmax)
irlba_time1 <- proc.time() - tstart
maxsingval1 <- max(psvd$d)
minsingval1 <- min(psvd$d)
numsingval1 <- length(psvd$d)
nrmse1 <- NULL

# Compute the norm of the reconstruction error 
# ---------------------------------------------
if (!is.null(psvd$d)){
  tiger.re1 <- psvd$u %*% diag(psvd$d) %*% t(psvd$v)
  nrmse1 <- sqrt(sum((tiger - tiger.re1) ** 2) / nrmse)
  erra <- norm(t(t(psvd$u) %*% tiger) - psvd$v %*% diag(psvd$d),"2")
  errb <- norm( tiger %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
  errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
  errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
  errVU1 <- sqrt(errV^2+errU^2)
  err1 <- sqrt(erra^2+errb^2)
}

# 2. Call the program svt_irlba with energy = 0.99
# --------------------------------------------------
energy2 <- 0.99
tstart <- proc.time()
psvd <- svt_irlba(tiger,tol = tol,energy = energy2,psvdmax = min(m,n),psvd0=psvd)
irlba_time2 <- proc.time() - tstart
maxsingval2 <- max(psvd$d)
minsingval2 <- min(psvd$d)
numsingval2 <- length(psvd$d)
nrmse2 <- NULL

# Computing the norm of the reconstruction error.
# -----------------------------------------------
if (!is.null(psvd$d)){
  tiger.re2 <- psvd$u %*% diag(psvd$d) %*% t(psvd$v)
  nrmse2 <- sqrt(sum((tiger - tiger.re2) ** 2) / nrmse)
  erra <- norm(t(t(psvd$u) %*% tiger) - psvd$v %*% diag(psvd$d),"2")
  errb <- norm( tiger %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
  errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
  errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
  errVU2 <- sqrt(errV^2+errU^2)
  err2 <- sqrt(erra^2+errb^2)
}

# 3. Calling the program rsvd with default values and k = 100.
# ------------------------------------------------------------
tstart <- proc.time()
psvd <- rsvd(tiger,k=100)
rsvd_time3 <- proc.time() - tstart
maxsingval3 <- max(psvd$d)
minsingval3 <- min(psvd$d)
numsingval3 <- length(psvd$d)
nrmse3 <- NULL

# Computing the norm of the reconstruction error.
# -----------------------------------------------
if (!is.null(psvd$d)){
  tiger.re3 <- psvd$u %*% diag(psvd$d) %*% t(psvd$v)
  nrmse3 <- sqrt(sum((tiger - tiger.re3) ** 2) / nrmse)
  erra3 <- norm(t(t(psvd$u) %*% tiger) - psvd$v %*% diag(psvd$d),"2")
  errb3 <- norm( tiger %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
  errV3 <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
  errU3 <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
  errVU3 <- sqrt(errV3^2+errU3^2)
  err3 <- sqrt(erra3^2+errb3^2)
}

# 4. Calling the program rsvd with q = 3 and k = 100.
# ----------------------------------------------------
q4 <- 3
tstart <- proc.time()
psvd <- rsvd(tiger,k=100,q=q4)
rsvd_time4 <- proc.time() - tstart
maxsingval4 <- max(psvd$d)
minsingval4 <- min(psvd$d)
numsingval4 <- length(psvd$d)
nrmse4 <- NULL

# Computing the norm of the reconstruction error.
# -----------------------------------------------
if (!is.null(psvd$d)){
  tiger.re4 <- psvd$u %*% diag(psvd$d) %*% t(psvd$v)
  nrmse4 <- sqrt(sum((tiger - tiger.re4) ** 2) / nrmse)
  erra4 <- norm(t(t(psvd$u) %*% tiger) - psvd$v %*% diag(psvd$d),"2")
  errb4 <- norm( tiger %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
  errV4 <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
  errU4 <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
  errVU4 <- sqrt(errV4^2+errU4^2)
  err4 <- sqrt(erra4^2+errb4^2)
}

# 5. Calling the program rsvd with q = 15 and k = 100.
# ----------------------------------------------------
q5 <- 15
tstart <- proc.time()
psvd <- rsvd(tiger,k=100,q=q5)
rsvd_time5 <- proc.time() - tstart
maxsingval5 <- max(psvd$d)
minsingval5 <- min(psvd$d)
numsingval5 <- length(psvd$d)
nrmse5 <- NULL

# Computing the norm of the reconstruction error.
# -----------------------------------------------
if (!is.null(psvd$d)){
  tiger.re5 <- psvd$u %*% diag(psvd$d) %*% t(psvd$v)
  nrmse5 <- sqrt(sum((tiger - tiger.re5) ** 2) / nrmse)
  erra5 <- norm(t(t(psvd$u) %*% tiger) - psvd$v %*% diag(psvd$d),"2")
  errb5 <- norm( tiger %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
  errV5 <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
  errU5 <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
  errVU5 <- sqrt(errV5^2+errU5^2)
  err5 <- sqrt(erra5^2+errb5^2)
}
# Displaying the results.
# -----------------------
cat(sprintf(" \n"))
cat(sprintf("Example 3.3: Compute energy percentages for the image tiger \n"))
cat(sprintf("*******************************************************************\n"))
cat(sprintf("Matrix tiger (%d x %d): \n",m,n))
cat(sprintf("SVT_IRLBA\n"))
cat(sprintf("------------------------\n"))
cat(sprintf("energy percentage %0.5g  \n",energy1))
cat(sprintf("   max. singval =  %0.5g\n",maxsingval1))
cat(sprintf("   min. singval =  %0.5g\n",minsingval1))
cat(sprintf("   # singvals = %0.5g\n",numsingval1))
cat(sprintf("   norm reconstruction error = %0.5g\n",nrmse1))
cat(sprintf("   cputime    = %0.5g secs\n",irlba_time1[["elapsed"]]))
cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err1))
cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU1))
cat(sprintf("------------------------\n"))
cat(sprintf("energy percentage %0.5g  \n",energy2))
cat(sprintf("   max. singval =  %0.5g\n",maxsingval2))
cat(sprintf("   min. singval =  %0.5g\n",minsingval2))
cat(sprintf("   # singvals = %0.5g\n",numsingval2))
cat(sprintf("   norm reconstruction error = %0.5g\n",nrmse2))
cat(sprintf("   cputime    = %0.5g secs\n",irlba_time2[["elapsed"]]))
cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err2))
cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU2))
cat(sprintf(" \n"))
cat(sprintf("RSVD\n"))
cat(sprintf("------------------------\n"))
cat(sprintf("   q  default value \n"))
cat(sprintf("   max. singval =  %0.5g\n",maxsingval3))
cat(sprintf("   min. singval =  %0.5g\n",minsingval3))
cat(sprintf("   # singvals = %0.5g\n",numsingval3))
cat(sprintf("   norm reconstruction error = %0.5g\n",nrmse3))
cat(sprintf("   cputime    = %0.5g secs\n",rsvd_time3[["elapsed"]]))
cat(sprintf("   ||AV - US|| = %0.5g\n",errb3))
cat(sprintf("   ||A^TU - VS|| = %0.5g\n",erra3))
cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err3))
cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU3))
cat(sprintf("------------------------\n"))
cat(sprintf("   q =  %0.5g\n",q4))
cat(sprintf("   max. singval =  %0.5g\n",maxsingval4))
cat(sprintf("   min. singval =  %0.5g\n",minsingval4))
cat(sprintf("   # singvals = %0.5g\n",numsingval4))
cat(sprintf("   norm reconstruction error = %0.5g\n",nrmse4))
cat(sprintf("   cputime    = %0.5g secs\n",rsvd_time4[["elapsed"]]))
cat(sprintf("   ||AV - US|| = %0.5g\n",errb4))
cat(sprintf("   ||A^TU - VS|| = %0.5g\n",erra4))
cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err4))
cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU4))
cat(sprintf("------------------------\n"))
cat(sprintf("   q =  %0.5g\n",q5))
cat(sprintf("   max. singval =  %0.5g\n",maxsingval5))
cat(sprintf("   min. singval =  %0.5g\n",minsingval5))
cat(sprintf("   # singvals = %0.5g\n",numsingval5))
cat(sprintf("   norm reconstruction error = %0.5g\n",nrmse5))
cat(sprintf("   cputime    = %0.5g secs\n",rsvd_time5[["elapsed"]]))
cat(sprintf("   ||AV - US|| = %0.5g\n",errb5))
cat(sprintf("   ||A^TU - VS|| = %0.5g\n",erra5))
cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err5))
cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU5))
cat(sprintf("------------------------\n"))

# Uncomment below to preview tiger images
# ---------------------------------------
# par(mfrow=c(2,2))
# tiger_image <- 255*tiger
# tiger_image.re1 <- 255*tiger.re1
# tiger_image.re2 <- 255*tiger.re2
# i1 <- image(tiger_image, col = gray(0:255 / 255),xaxt="n",yaxt="n",main="Original Image")
# i2 <- image(tiger_image.re1, col = gray(0:255 / 255),xaxt="n",yaxt="n",main=paste("Energy = ",energy1*100,"%"))
# i3 <- image(tiger_image.re2, col = gray(0:255 / 255),xaxt="n",yaxt="n",main=paste("Energy = ",energy2*100,"%"))
}  