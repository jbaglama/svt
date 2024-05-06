#  DEMO FUNCTION FOR PAPER:
#  J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm for 
#  Computing a Partial Singular Value Decomposition Satisfying a Given 
#  Threshold", submitted for publication 2024.
# 
# TO RUN DEMO IN R:
# Set the working directory - see below   
#  > source("svt_demo.R")
#  > svt_demo()
#
# REQUIRED SOFTWARE:  
# svt_irlba.R   
# https://github.com/jbaglama/svt
#
# DATE LAST MODIFIED: 
# 5/1/24
#
# LANGUAGE:
# R version: 4.3.2 
# 
# AUTHORS: 
# James Baglama            email: jbaglama@uri.edu
# Jonathan Chávez-Casillas email: jchavezc@uri.edu
# Vasilije Perovic         email: perovic@uri.edu
# ---------------------------------------------------
svt_demo <- 
  function(){

# ** IMPORTANT **  
# Working directory where script and files are located
# ** NEEDS TO BE SET BY USER **
# for example    
# > workingdir <- getwd()
# > setwd(workingdir)     
#-----------------------------------------------------
  
# Needed Package names
#---------------------
packages <- c("irlba","Matrix","methods","R.matlab","SparseM")

# Installing packages not yet installed.
#---------------------------------------
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Loading packages.
#------------------
invisible(lapply(packages, library, character.only = TRUE))

# Loading svt_irlba.R script
#---------------------------
source("svt_irlba.R")

# Setting a seed to reproduce examples.
#--------------------------------------
set.seed(1)

cat(sprintf("**********************************************************\n"))
cat(sprintf("DEMO:\n"))                                                              
cat(sprintf(" psvd <- svt_irlba(A,PARAMETERS)\n")) 
cat(sprintf(" A hybrid function for computing all singular triplets above a given\n")) 
cat(sprintf(" threshold. The function repeatedly calls a Partial Singular Value  \n")) 
cat(sprintf(" Decomposition (PSVD) method (svds or irlba) and a block SVD power  \n")) 
cat(sprintf(" method to compute a PSVD of a matrix A with all the singular values\n")) 
cat(sprintf(" above a threshold (sigma).\n"))                 
cat(sprintf(" \n")) 
cat(sprintf(" Demo runs examples with two matrices: illc1033 and bibd_20_10 \n"))
cat(sprintf(" Required matrix download from https://sparse.tamu.edu/ \n"))
cat(sprintf(" Demo will try to download matrices if not already in path \n"))
cat(sprintf(" \n"))                                                         
cat(sprintf("PAPER: \n"))                                                         
cat(sprintf("  J. Baglama, J.Chavez-Casillas and V. Perovic, A Hybrid Algorithm \n"))
cat(sprintf("  for Computing a Partial Singular Value Decomposition Satisfying a \n"))
cat(sprintf("  Given Threshold, submitted for publication 2024.                 \n"))
cat(sprintf("**********************************************************\n"))

# The following two matrices from the SuiteSparse Matrix Collection 
# https://sparse.tamu.edu/ will be used for this demo. Matrices from 
# SuiteSparse Matrix must be loaded in the path or the function will try 
# to download them. 
# illc1033  (1033 x 320) 
#   - largest singval ~ 2.144 and smallest singval ~ 0.00011353
#   - Large # of repeated singvals ~ 1 starting around singval(106)
# bibd_20_10  (190 x 184,756) -> A = A + LR'; L & R are randn(*,10) matrices
#   - largest singval ~ 7,133  and smallest singval ~ 113.37 
#   - Large # of repeated singvals ~ 113.446 starting around singval(31)
# 
# --------------------------------------------------------------------------

# illc1033  (1033 x 320)
# ---------------------- 
name <-  "illc1033.mat"
if (!file.exists(name)){
  
   # Save url of matrix in using MATLAB format
   # -----------------------------------------
   url <- "https://suitesparse-collection-website.herokuapp.com/mat/HB/illc1033.mat"
    
   # Download to working directory with given filename.
   # --------------------------------------------------
   download.file(url, paste0(workingdir,"/",name))
}
if (!file.exists(name))
  cat(sprintf("FAILED getting %s from web and does not exist in path\n",name))
 
if (file.exists(name)){
  
   # Loading in the file into the workspace - complete problem set and pull out the matrix only.
   # -------------------------------------------------------------------------------------------
   A <- readMat(name)
   A <- A$Problem[2][[1]]
  
   # A is a 1033 x 320 matrix. The code svt_irlba runs faster when m < n. 
   # Therefore, work with the transpose of A. This is not required, but will 
   # speed up the overall cpu time. 
   # ----------------------------------------------------------------------
   A <- t(A) 
   m <- nrow(A)
   n <- ncol(A)
 
   # Initializing the parameters.
   # ----------------------------
   k <- 6               # Initial number of singular triplets
   sigma <- 1.5         # Threshold 
   tol <- 1e-8          # Tolerance used for convergence in irlba or svds
   incre <- 5           # Initial increment
   pwrsvd <- 0          # Number of iters of a block SVD power method
   psvdmax <- min(m,n)  # Maximum "size" of the output PSVD [U,S,V]
   display <- 0         # Displays some iters diagnostic
   kmax <- min(m,n)     # Maximum value k can reach
   cmmf <- FALSE        # Use irlba mult parameter
  
   # Calling the program svt_irlba.
   # ------------------------------
   tstart <- proc.time()
   psvd <- svt_irlba(A,sigma = sigma,tol = tol,psvdmax = psvdmax,kmax = kmax,
                    pwrsvd = pwrsvd,incre = incre,k = k,cmmf = cmmf)
   irlba_time <- proc.time() - tstart
   
   # Computing the errors.
   # ---------------------
   if (!is.null(psvd$d)){
      err1 <- norm(t(t(psvd$u) %*% A) - psvd$v %*% diag(psvd$d),"2")
      err2 <- norm( A %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
       err <- sqrt(err1^2+err2^2)
      errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
      errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
     errVU <- sqrt(errV^2+errU^2)
   }
   
   # Displaying the results.
   # -----------------------
   cat(sprintf(" \n"))
   cat(sprintf("Example 1: Compute all singvals above sigma = %0.5g\n",sigma))
   cat(sprintf("*******************************************************************\n"))
   cat(sprintf("Matrix A (%d x %d): %s\n",m,n,name))
   cat(sprintf(" Parameters:\n"))
   cat(sprintf("   sigma = %0.5g\n",sigma))
   cat(sprintf("   tol = %0.5g\n ",tol))
   cat(sprintf("   psvdmax = %d\n",psvdmax))
   cat(sprintf("   kmax = %d\n",kmax))
   cat(sprintf("   pwrsvd = %d\n",pwrsvd))
   cat(sprintf("   incre = %d\n",incre))
   cat(sprintf("   k = %d\n",k))
   cat(sprintf("-------------------\n"))
   cat(sprintf("SVT_IRLBA\n")) 
   cat(sprintf("   FLAG = %d\n",psvd$flag))
   cat(sprintf("   max. singval =  %0.5g\n",max(psvd$d)))
   cat(sprintf("   min. singval (>=%0.5g)  =  %0.5g\n",sigma,min(psvd$d)))
   cat(sprintf("   # singvals = %0.5g\n",length(psvd$d)))
   cat(sprintf("   cputime    = %0.5g secs\n",irlba_time[["elapsed"]]))
   if (!is.null(psvd$d)){
      cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err))
      cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU))
   }
   cat(sprintf("**********************************************************\n"))
   
   # Continuing to compute more singular values.
   # -------------------------------------------
   sigma <- 0.9 
   
   # Calling the program svt_irlba.
   # ------------------------------
   tstart <- proc.time()
   psvd <- svt_irlba(A,sigma = sigma,tol = tol,psvdmax = psvdmax,kmax = kmax,
                     pwrsvd = pwrsvd,incre = incre,k = k,cmmf = cmmf,psvd0=psvd)
   irlba_time <- proc.time() - tstart
   
   # Computing the errors.
   # ---------------------
   if (!is.null(psvd$d)){
     err1 <- norm(t(t(psvd$u) %*% A) - psvd$v %*% diag(psvd$d),"2")
     err2 <- norm( A %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
     err <- sqrt(err1^2+err2^2)
     errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
     errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
     errVU <- sqrt(errV^2+errU^2)
   }
   
   # Displaying the results.
   # -----------------------
   cat(sprintf(" \n"))
   cat(sprintf("Example 2: Continue Example 1: Compute all singvals above sigma =  %0.5g\n",sigma))
   cat(sprintf("*******************************************************************\n"))
   cat(sprintf("Matrix A (%d x %d): %s\n",m,n,name))
   cat(sprintf(" Parameters:\n"))
   cat(sprintf("   sigma = %0.5g\n",sigma))
   cat(sprintf("   tol = %0.5g\n",tol))
   cat(sprintf("   psvdmax = %d\n",psvdmax))
   cat(sprintf("   kmax = %d",kmax))
   cat(sprintf("   pwrsvd = %d\n",pwrsvd))
   cat(sprintf("   incre = %d\n",incre))
   cat(sprintf("   k = %d\n",k))
   cat(sprintf("-------------------\n"))
   cat(sprintf("SVT_IRLBA\n")) 
   cat(sprintf("   FLAG = %d\n",psvd$flag))
   cat(sprintf("   max. singval =  %0.5g\n",max(psvd$d)))
   cat(sprintf("   min. singval (>=%0.5g)  =  %0.5g\n",sigma,min(psvd$d)))
   cat(sprintf("   # singvals = %0.5g\n",length(psvd$d)))
   cat(sprintf("   cputime    = %0.5g secs\n",irlba_time[["elapsed"]]))
   if (!is.null(psvd$d)){
     cat(sprintf("  sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err))
     cat(sprintf("  sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU))
   }
   cat(sprintf("**********************************************************\n"))
   
   # Changing the tolerance and threshold sigma.
   # -------------------------------------------
   tol <- 1e-4
   sigma <- 1.2
 
   # Calling the program svt_irlba.
   # ------------------------------
   tstart <- proc.time()
   psvd <- svt_irlba(A,sigma = sigma,tol = tol,psvdmax = psvdmax,kmax = kmax,
                     pwrsvd = pwrsvd,incre = incre,k = k,cmmf = cmmf)
   irlba_time <- proc.time() - tstart
   
   # Computing the errors.
   # ---------------------
   if (!is.null(psvd$d)){
     err1 <- norm(t(t(psvd$u) %*% A) - psvd$v %*% diag(psvd$d),"2")
     err2 <- norm( A %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
     err <- sqrt(err1^2+err2^2)
     errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
     errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
     errVU <- sqrt(errV^2+errU^2)
   }
   
   # Displaying the results.
   # -----------------------
   cat(sprintf(" \n"))
   cat(sprintf("Example 3: Compute all singvals above sigma =  %0.5g\n",sigma))
   cat(sprintf("*******************************************************************\n"))
   cat(sprintf("Matrix A (%d x %d): %s\n",m,n,name))
   cat(sprintf(" Parameters:\n"))
   cat(sprintf("   sigma = %0.5g\n",sigma))
   cat(sprintf("   tol = %0.5g\n",tol))
   cat(sprintf("   psvdmax = %d\n",psvdmax))
   cat(sprintf("   kmax = %d\n",kmax))
   cat(sprintf("   pwrsvd = %d\n",pwrsvd))
   cat(sprintf("   incre = %d\n",incre))
   cat(sprintf("   k = %d\n",k))
   cat(sprintf("-------------------\n"))
   cat(sprintf("SVT_IRLBA\n")) 
   cat(sprintf("   FLAG = %d\n",psvd$flag))
   cat(sprintf("   max. singval =  %0.5g\n",max(psvd$d)))
   cat(sprintf("   min. singval (>=%0.5g)  =  %0.5g\n",sigma,min(psvd$d)))
   cat(sprintf("   # singvals = %0.5g\n",length(psvd$d)))
   cat(sprintf("   cputime    = %0.5g secs\n",irlba_time[["elapsed"]]))
   if (!is.null(psvd$d)){
     cat(sprintf("  sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err))
     cat(sprintf("  sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU))
   }
   cat(sprintf("**********************************************************\n"))
   
   # Continuing to compute more singular values.
   # -------------------------------------------
   sigma <- 0.9 
   pwrsvd <- 1
   
   # Calling the program svt_irlba.
   # ------------------------------
   tstart <- proc.time()
   psvd <- svt_irlba(A,sigma = sigma,tol = tol,psvdmax = psvdmax,kmax = kmax,
                     pwrsvd = pwrsvd,incre = incre,k = k,cmmf = cmmf)
   irlba_time <- proc.time() - tstart
   
   # Computing the errors.
   # ---------------------
   if (!is.null(psvd$d)){
     err1 <- norm(t(t(psvd$u) %*% A) - psvd$v %*% diag(psvd$d),"2")
     err2 <- norm( A %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
     err <- sqrt(err1^2+err2^2)
     errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
     errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
     errVU <- sqrt(errV^2+errU^2)
   }
   
   # Displaying the results.
   # -----------------------
   cat(sprintf(" \n"))
   cat(sprintf("Example 4: Continue Example 3 with power svd iteration to compute all singvals above sigma =  %0.5g\n",sigma))
   cat(sprintf("*******************************************************************\n"))
   cat(sprintf("Matrix A (%d x %d): %s\n",m,n,name))
   cat(sprintf(" Parameters:\n"))
   cat(sprintf("   sigma = %0.5g\n",sigma))
   cat(sprintf("   tol = %0.5g\n",tol))
   cat(sprintf("   psvdmax = %d\n",psvdmax))
   cat(sprintf("   kmax = %d\n",kmax))
   cat(sprintf("   pwrsvd = %d\n",pwrsvd))
   cat(sprintf("   incre = %d\n",incre))
   cat(sprintf("   k = %d\n",k))
   cat(sprintf("-------------------\n"))
   cat(sprintf("SVT_IRLBA\n")) 
   cat(sprintf("   FLAG = %d\n",psvd$flag))
   cat(sprintf("   max. singval =  %0.5g\n",max(psvd$d)))
   cat(sprintf("   min. singval (>=%0.5g)  =  %0.5g\n",sigma,min(psvd$d)))
   cat(sprintf("   # singvals = %0.5g\n",length(psvd$d)))
   cat(sprintf("   cputime    = %0.5g secs\n",irlba_time[["elapsed"]]))
   if (!is.null(psvd$d)){
     cat(sprintf("  sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err))
     cat(sprintf("  sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU))
   }
   cat(sprintf("**********************************************************\n"))
}
   
# bibd_20_10  Size: 190 x 184,756
# ------------------------------ 
name <-  "bibd_20_10.mat"
if (!file.exists(name)){
  
  # Saving url of matrix in using MATLAB format.
  # --------------------------------------------
  url <- "https://suitesparse-collection-website.herokuapp.com/mat/JGD_BIBD/bibd_20_10.mat"
  
  # Downloading to working directory with given filename.
  # -----------------------------------------------------
  download.file(url, paste0(workingdir,"/",name))
}
if (!file.exists(name))
  cat(sprintf("FAILED getting %s from web and does not exist in path\n",name))   
   

if (file.exists(name)){
  
  # Loading the file into R - complete problem set and pull out the matrix only.
  # ----------------------------------------------------------------------------
  A <- readMat(name)
  A <- A$Problem[3][[1]]
  m <- nrow(A) 
  n <- ncol(A) 
  L <- matrix(rnorm(m*10),nrow=m) 
  R <- matrix(rnorm(n*10),nrow=n) 
  LR <- tcrossprod(L,R) 
  A  <- A + LR
  
  # Initializing the  parameters.
  # -----------------------------
  k <- 6               # Initial number of singular triplets
  sigma <- 4000        # Threshold 
  tol <- 1e-8          # Tolerance used for convergence in irlba or svds
  incre <- 5           # Initial increment
  pwrsvd <- 0          # Number of iters of a block SVD power method
  psvdmax <- min(m,n)  # Maximum "size" of the output PSVD [U,S,V]
  display <- 0         # Displays some iters diagnostic
  kmax <- min(m,n)     # Maximum value k can reach
  cmmf <- FALSE        # Use irlba mult parameter
  
  # Calling the program svt_irlba.
  # ------------------------------
  tstart <- proc.time()
  psvd <- svt_irlba(A,sigma = sigma,tol = tol,psvdmax = psvdmax,kmax = kmax,
                    pwrsvd = pwrsvd,incre = incre,k = k,cmmf = cmmf)
  irlba_time <- proc.time() - tstart
  
  # Computing the errors.
  # ---------------------
  if (!is.null(psvd$d)){
    err1 <- norm(t(t(psvd$u) %*% A) - psvd$v %*% diag(psvd$d),"2")
    err2 <- norm( A %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
    err <- sqrt(err1^2+err2^2)
    errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
    errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
    errVU <- sqrt(errV^2+errU^2)
  }
  
  # Displaying the results.
  # -----------------------
  cat(sprintf(" \n"))
  cat(sprintf("Example 5: Compute all singvals above sigma = %0.5g\n",sigma))
  cat(sprintf("*******************************************************************\n"))
  cat(sprintf("Matrix A (%d x %d): %s\n",m,n,name))
  cat(sprintf(" Parameters:\n"))
  cat(sprintf("   sigma = %0.5g\n",sigma))
  cat(sprintf("   tol = %0.5g\n",tol))
  cat(sprintf("   psvdmax = %d\n",psvdmax))
  cat(sprintf("   kmax = %d\n",kmax))
  cat(sprintf("   pwrsvd = %d\n",pwrsvd))
  cat(sprintf("   incre = %d\n",incre))
  cat(sprintf("   k = %d\n",k))
  cat(sprintf("-------------------\n"))
  cat(sprintf("SVT_IRLBA\n")) 
  cat(sprintf("   FLAG = %d\n",psvd$flag))
  cat(sprintf("   max. singval =  %0.5g\n",max(psvd$d)))
  cat(sprintf("   min. singval (>=%0.5g)  =  %0.5g\n",sigma,min(psvd$d)))
  cat(sprintf("   # singvals = %0.5g\n",length(psvd$d)))
  cat(sprintf("   cputime    = %0.5g secs\n",irlba_time[["elapsed"]]))
  if (!is.null(psvd$d)){
    cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err))
    cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU))
  }
  cat(sprintf("**********************************************************\n"))
  
  # Continuing to compute more singular values.
  # -------------------------------------------
  sigma <- 200  
  
  # Calling the program svt_irlba.
  # ------------------------------
  tstart <- proc.time()
  psvd <- svt_irlba(A,sigma = sigma,tol = tol,psvdmax = psvdmax,kmax = kmax,
                    pwrsvd = pwrsvd,incre = incre,k = k,cmmf = cmmf,psvd0=psvd)
  irlba_time <- proc.time() - tstart
  
  # Computing the errors.
  # ---------------------
  if (!is.null(psvd$d)){
    err1 <- norm(t(t(psvd$u) %*% A) - psvd$v %*% diag(psvd$d),"2")
    err2 <- norm( A %*% psvd$v  - psvd$u %*% diag(psvd$d),"2")
    err <- sqrt(err1^2+err2^2)
    errV <-  norm(t(psvd$v) %*% psvd$v - diag(ncol(psvd$v)))
    errU <- norm(t(psvd$u) %*% psvd$u - diag(ncol(psvd$u)))
    errVU <- sqrt(errV^2+errU^2)
  }
  
  # Displaying the results.
  # -----------------------
  cat(sprintf(" \n"))
  cat(sprintf("Example 6: Continue Example 5: Compute all singvals above sigma =  %0.5g\n",sigma))
  cat(sprintf("*******************************************************************\n"))
  cat(sprintf("Matrix A (%d x %d): %s\n",m,n,name))
  cat(sprintf(" Parameters:\n"))
  cat(sprintf("   sigma = %0.5g\n",sigma))
  cat(sprintf("   tol = %0.5g\n",tol))
  cat(sprintf("   psvdmax = %d\n",psvdmax))
  cat(sprintf("   kmax = %d\n",kmax))
  cat(sprintf("   pwrsvd = %d\n",pwrsvd))
  cat(sprintf("   incre = %d\n",incre))
  cat(sprintf("   k = %d\n",k))
  cat(sprintf("-------------------\n"))
  cat(sprintf("SVT_IRLBA\n")) 
  cat(sprintf("   FLAG = %d\n",psvd$flag))
  cat(sprintf("   max. singval =  %0.5g\n",max(psvd$d)))
  cat(sprintf("   min. singval (>=%0.5g)  =  %0.5g\n",sigma,min(psvd$d)))
  cat(sprintf("   # singvals = %0.5g\n",length(psvd$d)))
  cat(sprintf("   cputime    = %0.5g secs\n",irlba_time[["elapsed"]]))
  if (!is.null(psvd$d)){
    cat(sprintf("   sqrt(||AV - US||^2+||A^TU - VS||^2) = %0.5g\n",err))
    cat(sprintf("   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n",errVU))
  }
  cat(sprintf("**********************************************************\n"))
}
return(invisible())
}