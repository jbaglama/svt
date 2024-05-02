#  DESCRIPTION:
#  Implements Algorithm 1 in reference 1.
#  A hybrid function for computing all singular triplets above a given threshold. 
#  The function repeatedly calls a Partial Singular Value Decomposition (PSVD)
#  method (irlba) and a block SVD power method to compute a PSVDof a matrix A. 
#  If a user threshold value (sigma) is specified, the function computes all the
#  singular values above such threshold (sigma). However, if an energy percentage
#  (<=1) (energy) is specified, this function computes  the energy percentage. 
#  That is, computes an r-rank approximation A_r of A where 
#  (||A_r||_F/||A||_F)^2 >= energy. If neither are given then the function 
#  returns a k-PSVD with k=6 by default. 
#
# *** Required libraries (and dependencies): irlba, Matrix, and methods  ***
#
# INPUT:
#    A - An m x n numeric sparse real matrix (matrix-product input is not available)
#    PARAMETERS:
#        sigma  - Singular value threshold (>= 0). If missing, the function returns a k-PSVD.
#         (or)
#       energy  - Energy percentage (decimal <= 1) - r-rank approximation A_r of  
#                 A where (||A_r||_F/||A||_F)^2 >= energy. 
#                 It cannot be combined with a sigma value. (default: NULL)
#        psvd0  - List of an already computed PSVD of A. (default: NULL)
#                     A %*% psvd0$v = psvd0$u %*% diag(psvd0$d)
#                  t(A) %*% psvd0$u = psvd0$v %*% diag(psvd0$d)
#        tol   -  Tolerance used for convergence in the irlba routine. (default: 1e-8)
#        k     -  Initial number of singular triplets. (default: 6)
#        incre -  Increment added to k. (internally doubled each iteration) (default: 5)
#        kmax  -  Maximum value k can reach. (default: min(0.1*min(m,n),100))
#         p0   -  Starting vector for irlba. (default: p = rnorm(n))
#     psvdmax  -  Maximum dimension of the output of PSVD. The output 
#                 psvd will contain the input psvd if given.
#                 (default: max(min(100+size(S0),min(n,m)),k))
#                 *NOTE* This function will allocate memory for the full matrices 
#                        and will run out of memory (return system error) on 
#                        initialization for very large A and psvdmax value. 
#      pwrsvd  -  If set to an integer > 0 then, on each iteration, it will perform 
#                 pwrsvd iterations of a block SVD power method with the
#                 output from irlba. If set to 0 only one iteration is performed when 
#                 required. (e.g. loss of orthogonality of basis vectors)(default: 0)
#     display  -  If set to 1 displays iteration diagnostic information. (default: 0)
#      cmmf   -   Logic variable, if TRUE use the internal custom matrix multiplication 
#                 function (cmmf) to compute the matrix-product deflation. If FALSE, the
#                 routine will use the irlba mult parameter - see irlba documentation
#                 for details. (default: FALSE) 
#
#  OUTPUT:
#   Returns a list with entries:
#      u  -  Matrix of left singular vectors. 
#      d  -  Vector of singular values in decreasing order.
#      v  -  Matrix of right singular vectors.
#            *NOTE* - The output u,d,v will include input values from psvd if given.
#   flag  - 0 successful output - either threshold (sigma) met or energy 
#             percentage satisfied.  
#           1 if irlba iteration fails to compute any singular triplets - output 
#             last values for u,d,v.
#           2 psvdmax was reached before threshold (sigma) was met or energy 
#             percentage was satisfied - it outputs the last values for u,d,v.
#           3 if no singular valuess above threshold (sigma) exist - it outputs
#             NULL for u,d,v.
#
#  EXAMPLES:
#
#      1. Compute all singular triplets with singular values exceeding 5.1:
#          > psvd <- svt_irlba(A,sigma=5.1)
#
#      2. Compute all singular triplets with singular values exceeding 5.1:
#         given an initial PSVD psvd:
#          > psvd <- svt_irlba(A,sigma=5.1,psvd0=psvd)
#
#      3. Compute all singular triplets with singular values exceeding 1.1,
#         within tolerance of 1e-10 and provide a maximum of 20 singular triplets:
#          > psvd <- svt_irlba(A,sigma=1.1,tol=1e-10,psvdmax=20)
#
#      4. Compute the top 6 singular triplets and then continue:
#          > psvd <- svt_irlba(A)  
#         Assume that based on the output the desired threshold is 100 and set
#         psvdmax, the number of singular triplets to 20: 
#          > psvd <- svt_irlba(A,sigma=100,psvd0=psvd,psvdmax=20)
#         The user can check if any multiple singular values have been missed by
#         calling the function again with the same parameters and threshold, but
#         with the output from the previous call
#
#      5. Compute the energy percentage 0.9:
#          > psvd <- svt_irlba(A,energy = 0.9)
#    
#  DATE LAST MODIFIED: 
#  5/1/24
#  
#  AUTHORS: 
#   James Baglama            email: jbaglama@uri.edu
#   Jonathan Chávez-Casillas email: jchavezc@uri.edu
#   Vasilije Perovic         email: perovic@uri.edu
#
# REFERENCES:
#  1. J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm for 
#     Computing a Partial Singular Value Decomposition Satisfying a Given 
#     Threshold", submitted for publication 2024.
#  2. J. Baglama and V. Perovic, "Explicit Deflation in Golub–Kahan–Lanczos 
#     Bidiagonalization Methods, ETNA, Vol. 58, (2023), pp. 164–176.
#  3. B.W. Lewis, J.Baglama, L. Reichel, "The irlba Package", (2021)
#     https://cran.r-project.org/web/packages/irlba/
# --------------------------------------------------------------------------
svt_irlba <- 
function(A,               # m x n real matrix - complex currently not supported  
         sigma = Inf,     # Threshold 
         energy = NULL,   # Energy percentage (decimal <= 1)
         psvd0 = NULL,    # List of an already computed PSVD
         tol = 1e-8,      # Tolerance used for convergence in the irlba routine
         psvdmax = NULL,  # Maximum size of computed PSVD - default value set later 
         display = 0,     # Displays iteration diagnostic information
         p0 = NULL,       # Starting vector for irlba
         kmax = NULL,     # Maximum value k can reach over all iterations - default value set later
         pwrsvd = 0,      # Number of iterations of a simple block SVD power method
         incre = 5,       # Initial increment (internally doubled each iteration) to reach threshold
         k = 6,           # Number of initial singular triplets used in irlba
         cmmf = FALSE     # Use the internal custom matrix multiplication function (cmmf)
    ){

# Libraries irlba, Matrix, and methods are required.
# -------------------------------------------------- 
require(irlba)
require(Matrix)
require(methods)

# Predefined a two-norm function for vectors. 
# -------------------------------------------
norm2 <- function(x) drop(sqrt(crossprod(x)))  

# NOTE:
# Parameters not available by user for input into irlba:
# CENTER, SCALE, SHIFT, MAXIT, NV, NU, WORK, REORTH, V, RIGHT_ONLY,VERBOSE, 
# MULT, FASTPATH, SVTOL, SMALLEST - many of these are set by this function.
# Parameters (e.g. CENTER, SHIFT,SCALE) may become available in future versions.
# For details:  B.W. Lewis, J.Baglama, L. Reichel, "The irlba Package", (2021)
#     https://cran.r-project.org/web/packages/irlba/
# -----------------------------------------------------------------------------

# Setting clock time for display option.
# --------------------------------------
tstart <- proc.time()

# Checking input matrix A.
#-------------------------
 if (is.complex(A)) 
   stop("svt_irlba is not set up for complex matrices")
 if (typeof(A) == "integer") A <- A + 0.0
 typeA <- NULL
 if (is(A,'Matrix')) typeA <- "Matrix"
 if (is(A,'matrix')) typeA <- "matrix" 
 if (is.null(typeA)) 
   stop("input matrix must be numeric or sparse matrix")
 m <- nrow(A)
 n <- ncol(A)
 interchange <- 0
 
# Checking for some errors of parameters and set some internal variables.
#------------------------------------------------------------------------
 min_mn <- min(m,n)
 if (is.null(kmax)) 
   kmax <- max(floor(min(0.1*min_mn,100)),k) 
 if (!is.numeric(kmax)  || kmax < 0 || kmax > min_mn || kmax != floor(kmax))
   stop(" incorrect value for kmax  ")
 if (kmax > min_mn - 1) kmax <- min_mn - 2
 if (!is.numeric(k) || k <= 0 || k != floor(k))
   stop(" incorrect value for k ")
 k_org <- k
 if (!is.numeric(incre) || incre <= 0 || incre != floor(incre)) 
   stop(" incorrect value for incre ")
 if (!is.numeric(pwrsvd) || pwrsvd < 0 || pwrsvd != floor(pwrsvd))
   stop(" incorrect value for pwrsvd ")
 if (!is.numeric(display) || (display != 0 & display != 1))
   stop(" incorrect value for display")
 if (!is.numeric(sigma) || sigma < 0)
   stop(" incorrect value for sigma ")
 if (!is.null(energy)){
    if (!is.numeric(energy) || energy <= 0 || energy > 1)
       stop("Incorrect value for energy - 0 < energy <= 1.")
    if (sigma != Inf)
       stop("Cannot have a value for sigma and energy.")
    FronormA <- norm(A,"F")^2
 }
 
 if (is.null(p0)) p0 <- rnorm(n)
 if (length(p0) != n)
   stop(" incorrect size for p0 ")
 
 # Checking tolerance and set internal values.
 # -------------------------------------------
 if (tol < 0) stop(" tol must be non-negative ")
 if (tol < .Machine$double.eps)
    tol <- .Machine$double.eps
 eps12 <- .Machine$double.eps^0.5
 eps34 <- .Machine$double.eps^0.75
 Smax <- 1
 Slen <- 0
 
 # Checking the sizes of the input psvd:  
 # This assumes A is m x n and that the PSVD has form: 
 #    A %*% psvd0$v = psvd0$u %*% diag(psvd0$d)
 # t(A) %*% psvd0$u = psvd0$v %*% diag(psvd0$d)
 # ---------------------------------------------------
 if (is.list(psvd0) && !is.null(psvd0$v) && !is.null(psvd0$u) && !is.null(psvd0$d)){
   mu <- nrow(psvd0$u) 
   ku <- ncol(psvd0$u)
   nv <- nrow(psvd0$v) 
   kv <- ncol(psvd0$v)
   Slen <- length(psvd0$d)
   if (nv != n || mu != m || kv != ku || Slen != kv) 
     stop(" incorrect dimensions of input - psvd0$v, psvd0$u, psvd0$d ")
   if (is.unsorted(psvd0$d)){
     ss <- sort(psvd0$d,decreasing = TRUE,index.return=TRUE)
     psvd0$d <- ss$x
     psvd0$u <- psvd0$u[,ss$ix]
     psvd0$v <- psvd0$v[,ss$ix]
   }
   Smax <- psvd0$d[1]
   if (psvd0$d[Slen] < sigma && sigma != Inf)
     return(list(d=psvd0$d, v=psvd0$v, u=psvd0$u, flag=0))
 }
 
 # Setting the initial value for psvdmax.
 # --------------------------------------
 S0len <- Slen
 SlenL <- 1:Slen
 if (is.null(psvdmax)) 
   psvdmax <- min(100,min_mn-Slen)
 else
   psvdmax <- min(psvdmax,min_mn-Slen)
 if (k>= min_mn - Slen)
   stop(" incorrect value for k ")
 psvdmax <- max(psvdmax,k)
 if (psvdmax <= 0){
   warning("Input dim of psvd is >= min(m,n)) ",immediate.=TRUE)
   if (is.list(psvd0)) return(list(d=psvd0$d, v=psvd0$v, u=psvd0$u, flag=2))
   return(list(d=NULL, v=NULL, u=NULL, flag = 3))
 }

 # Define the matrix-vector product deflation functions. 
 # -----------------------------------------------------
 # Define the custom matrix multiplication function (cmmf).
 if (cmmf) {
   # Matrix vector product with (I - UU')A.
   if (m <= n) {
     if (is(A,'Matrix')){
        setClass("AUV_class", contains=typeA, slots=c(A_mat=typeA,U="matrix",V="matrix"))
     }
     else{
       setClass("AUV_class", contains=typeA, slots=c(A_mat=typeA,U="matrix",V="matrix",Dim="numeric",dim="numeric"))
     }
      setMethod("%*%", signature(x="AUV_class", y="numeric"),
             function(x ,y){
             z <-  drop(x@A_mat %*% y)
             z - (x@U[,SlenL] %*% drop(z %*% x@U[,SlenL]))}) # z = (I - UU')*A*y
      setMethod("%*%", signature(x="numeric", y="AUV_class"),
           function(x ,y){
             z <- x - drop(y@U[,SlenL] %*% drop(x %*% y@U[,SlenL]))
             z %*% y@A_mat}) # z = A'(I-UU')*x
   }
   else{
      # Matrix vector product with A(I - VV'). 
      interchange <- 1
      pwrsvd <- max(1,pwrsvd)
      if (is(A,'Matrix')){
          setClass("AUV_class", contains=typeA, slots=c(A_mat=typeA,U="matrix",V="matrix"))
      }
      else{
        setClass("AUV_class", contains=typeA, slots=c(A_mat=typeA,U="matrix",V="matrix",Dim="numeric",dim="numeric"))
      }
      setMethod("%*%", signature(x="AUV_class", y="numeric"),
                function(x ,y){
                  z <- y - drop(x@V[,SlenL] %*% drop(y %*% x@V[,SlenL]))
                  x@A_mat %*% z}) # z = A(I - VV')*y
      setMethod("%*%", signature(x="numeric", y="AUV_class"),
                function(x ,y){
                  z <- drop(x %*% y@A_mat)
                  z - (y@V[,SlenL] %*% drop(z %*% y@V[,SlenL]))})# z = (I - VV')A'x
    }
    C <- new("AUV_class", 
                 A_mat = A, 
                 V = matrix(0.0, n, psvdmax+Slen),
                 U = matrix(0.0, m, psvdmax+Slen), Dim = dim(A))
    if (is(A,'matrix')){
       C@dim <- C@Dim
       # Checking the dimensions of C.
       if (any(dim(C) != dim(A)))
       stop("Error setting up custom matrix multiplication - try cmmf = FALSE")
    }
 } 
 # Using the irlba mult (Deprecated) parameter option.  
 else {
   # Matrix vector product with (I - UU')A.
   if (m <= n) {
      Amult <- function(x, y){
        # Checking if y is a  vector -> z = (I - UU')*A*y.
        if (is.vector(y)) {
           z <-  drop(x %*% y)
           return(z - (C@U[,SlenL] %*% drop(z %*% C@U[,SlenL])))
        }
        # else y is the matrix -> z = A'(I-UU')*x.
        z <- x - drop(C@U[,SlenL] %*% drop(x %*% C@U[,SlenL]))
        z %*% y
      }
   }
   else{
     # Matrix vector product with A(I - VV'). 
     interchange <- 1
     pwrsvd <- max(1,pwrsvd)
     Amult <- function(x, y){
       # Checking if y is a  vector -> z = A(I - VV')*y.
       if (is.vector(y)) {
         z <- y - drop(C@V[,SlenL] %*% drop(y %*% C@V[,SlenL]))
         return(x %*% z)
       }
       # else y is the matrix -> z = (I - VV')A'x.
       z <- drop(x %*% y)
       z - (C@V[,SlenL] %*% drop(z %*% C@V[,SlenL]))# z = (I - VV')A'x
     }
   }
   setClass("UV_class", slots=c(U="matrix",V="matrix"))
   C <- new("UV_class", 
            V = matrix(0.0, n, psvdmax+Slen),
            U = matrix(0.0, m, psvdmax+Slen))
  }
 
 # Allocating memory for V, U, and S. NOTE: If psvdmax+Slen is too large R will 
 # return an out of memory error here. If a PSVD is inputted, set the values in
 # V, U, S.
 # -----------------------------------------------------------------------------
 S <- matrix(0.0, psvdmax+Slen,1)
 if (Slen > 0){
   C@V[,1:Slen] <- psvd0$v
   C@U[,1:Slen] <- psvd0$u
     S[1:Slen]  <- psvd0$d
 }

# ---------------------#
# BEGIN: MAIN PROGRAM  #
# ---------------------#
 maxit <- 100          # Initial max # of iterations for irlba.
 work  <- max(3*k,20)  # Initial max dimension of subspace for irlba.
 iter <- 1             # Iteration count.

 # Checking to see if a power SVD iteration should be performed with input psvd.
 # -----------------------------------------------------------------------------
 if (Slen > 0){
   if (pwrsvd){
     if (interchange){
        psvd0 <- blocksvdpower(A,C@V[,1:Slen],pwrsvd,interchange)
     }
     else{
        psvd0 <- blocksvdpower(A,C@U[,1:Slen],pwrsvd,interchange)
     }
     Slen <- length(psvd0$d)
     C@U[,1:Slen] <- psvd0$u
     C@V[,1:Slen] <- psvd0$v
     S[1:Slen] <-  psvd0$d
   }
   p0 <- p0 - drop(C@V[,1:Slen] %*% drop(p0 %*% C@V[,1:Slen]))
 } 
 p0 <- p0/norm2(p0)
 
# Main while loop.
# ----------------
 while (iter < min_mn){
   
     # Resetting internal values on each iteration.
     # --------------------------------------------
     svd_iter <- 1 
     svd_fail <- 1
     orth_chk <- 0
     tStart_irlba <- proc.time()
     lksvdpw <- 0 
     blksvdpw_time <- 0
   
     # Calling the irlba function to get the singular triplets, with a 
     # maximum of 2 iterations - examples do not suggest more iterations 
     # are needed.
     # ------------------------------------------------------------------
     while (svd_iter <= 2 && svd_fail){
       if (Slen == 0){
          psvd <- suppressWarnings(irlba(A,nv=k,maxit=svd_iter*maxit,work=svd_iter*work,
                       tol=tol,svtol=Smax*min(eps12,tol),fastpath=FALSE,reorth=TRUE))
       } 
       else {
          if (cmmf) 
             psvd <- suppressWarnings(irlba(C,nv=k,maxit=svd_iter*maxit,work=svd_iter*work,
                        tol=tol,v=p0,svtol=Smax*min(eps12,tol),fastpath=FALSE,reorth=TRUE))
          else 
             psvd <- suppressWarnings(irlba(A,mult=Amult,nv=k,maxit=svd_iter*maxit,work=svd_iter*work,
                       tol=tol,v=p0,svtol= Smax*min(eps12,tol),fastpath=FALSE,reorth=TRUE))
       }
       conv <- 1:k
       svd_fail <-  0
        
       # Implies irlba failed to converge - reached the maximum number of iterations.
       # It is still needed to determine if any singular triplets converged. This 
       # occurs rarely, but needs to be addressed as well.
       # ----------------------------------------------------------------------------
       if (psvd$iter >= svd_iter*maxit){
          err_r <- t(t(psvd$u) %*% A) - psvd$v %*% diag(psvd$d)
          conv <- which(apply(err_r, 2, function(x) norm2(x)) < max(1,Smax)*tol,arr.ind = TRUE)
          svd_fail <- 1
          maxit <- 2*maxit
          if (length(conv) > 0){
             psvd$v <- psvd$v[,conv]
             psvd$u <- psvd$u[,conv]
             psvd$d <- psvd$d[conv]
             break # implies k > 0 singular triplets have been computed
          }
          svd_iter <- svd_iter+1
       }
     } # end while loop for irlba-
     irlba_time <- proc.time() - tStart_irlba
      
     # Checking if irlba failed to compute anymore singular triplets. 
     # --------------------------------------------------------------
     if (svd_fail && length(conv) == 0){
       if (display){
         cat(sprintf("IRLBA failed to return any singular values.\n"))
       }
       if (Slen > 0) return(list(d=S[1:Slen], v=C@V[,1:Slen], u=C@U[,1:Slen], flag=1))
       return(list(d=NULL, v=NULL, u=NULL, flag=3))
     }
        
     # Reset the k value based on the return from irlba
     # --------------------------------------------------
     k <- length(conv)
     comput_k <- k
     Slenk <- Slen+k
     Skindex <- (Slen+1):Slenk
      
     # Checking orthogonality among the "long" vectors - That is, 
     # the vectors not explicit orthogonalized with basis vectors.
     # orth_chk is based on Theorem 5 (page 31)
     # Larsen, Rasmus Munk. "Lanczos bidiagonalization with partial 
     # reorthogonalization." DAIMI Report Series 537 (1998).
     # -------------------------------------------------------------
     if (Slen > 0 && !svd_fail && pwrsvd == 0){
       if (interchange){
          orth_chk <- max(abs(t(psvd$u) %*% C@U[,1:Slen]))
       }
       else{
         orth_chk <- max(abs(t(psvd$v) %*% C@V[,1:Slen]))
       }
     }
     
     # Calling the block svd power method. This improves the orthogonality of 
     # the basis vectors and the overall error. However, this can be 
     # computationally expensive especially if A is large and/or number of 
     # computed singular triplets is large. 
     # If pwrsvd = 0, a determination is set by orth_chk.
     # ----------------------------------------------------------------------
     if (orth_chk > eps12/sqrt(Slen+k) || pwrsvd || svd_fail || psvd$d[k] < Smax*eps12){
       tStart_blk <- proc.time()
       if (interchange){
         C@V[,Skindex] <- psvd$v 
          psvd <- blocksvdpower(A,C@V[,1:Slenk],max(1,pwrsvd),interchange)
       }
       else{
         C@U[,Skindex] <- psvd$u 
         psvd <- blocksvdpower(A,C@U[,1:Slenk],max(1,pwrsvd),interchange)
       }
       Slenk <- length(psvd$d)
       Skindex <- 1:Slenk
       blksvdpw <- 1 
       blksvdpw_time = proc.time() - tStart_blk
     }
     
     # Updating the matrices U, V, S, and starting vector p0.
     # ------------------------------------------------------
     C@U[,Skindex] <- psvd$u
     C@V[,Skindex] <- psvd$v
     S[Skindex] <-  psvd$d
     p0 <- p0 - drop(psvd$v %*% drop(p0 %*% psvd$v))
     p0 <- p0/norm2(p0)
     
     # Updating the number of computed singular triplets.
     #---------------------------------------------------
     Slen <- Slenk
     SlenL <- 1:Slen
     
     # Checking for convergence.
     #--------------------------
     Smin <-  min(S[1:Slen])
     Smax <-  max(max(S[1:Slen]),Smax)
     statement1 <- 0
     statement2 <- 0
     if (is.null(energy)){
        statement1 <-  Smin < sigma
     }
     else{
        statement2 <- (norm2(S[1:Slen])^2/FronormA) >= energy
     }
     statement3 <-  Slen >= (psvdmax+S0len)
     statement4 <-  Slen >= min_mn
     if (display){
        cat(sprintf(" \n"))
        cat(sprintf(" IRLBA: CPU time = %0.5g secs. \n",irlba_time[["elapsed"]]))
        cat(sprintf(" IRLBA: computed singvals : %d   \n",comput_k))
        cat(sprintf(" Current size of PSVD of A: %d \n", Slen))
        cat(sprintf(" Approx. max. singval of A:  %0.5g \n",Smax))
        cat(sprintf(" Approx. min. singval of A: %0.5g \n", Smin))
        if (blksvdpw){
         cat(sprintf(" SVD power method: iterations = %d,CPU time = %0.5g secs.\n",
         max(1,pwrsvd),blksvdpw_time[["elapsed"]]))  
       }
     }
     if (statement1 || statement2 || statement3 || statement4){
        # Converged or max number of singular triplets reached
        # Resize S,V, U for output 
       C@U <- C@U[,1:Slen]
       C@V <- C@V[,1:Slen]
        S <- S[1:Slen]
        flag <- 0
        # perhaps change to all(diff(S) >= 0)
        if (is.unsorted(S)){
          ss <- sort(S,decreasing = TRUE,index.return=TRUE)
          S <- ss$x
          C@U <- C@U[,ss$ix]
          C@V <- C@V[,ss$ix]
        }
        Slen1 <-  min(Slen,min_mn,psvdmax+S0len)
        if (Slen1 < Slen) {
          C@U <- C@U[,1:Slen1]
          C@V <- C@V[,1:Slen1]
          S <- S[1:Slen1]
        }
        # Check if sigma is set to a threshold value
        if (sigma != Inf){
          # Already ordered largest to smallest.
          thres = which(S >= sigma,arr.ind = TRUE)
          C@U <- C@U[,thres]
          C@V <- C@V[,thres]
          S <- S[thres]
        }
        if (statement2){
          for (i in Slen:1){
            FronormS <- norm2(S[1:i])^2
            if (FronormS/FronormA < energy){
              S <- S[1:(i+1)] 
              C@V <- C@V[,1:(i+1)] 
              C@U <- C@U[,1:(i+1)] 
              break
            }
          }
        }
        if (!statement1 && !statement2 && statement3 && !statement4) flag <- 2
        if (length(S) == 0) {
           S <- NULL
           C@V <- NULL
           C@U <- NULL
           flag <- 3 
        }
        if (display){
           if (flag == 0){
             cat(sprintf(" Successful return - threshold met \n"))
           }
           tot_time <- proc.time() - tstart 
           cat(sprintf(" Total time for function svt_irlba: %0.5g secs. \n",tot_time[["elapsed"]]))
        }
        return(list(d=S, v=C@V, u=C@U, flag=flag))
     }
      
     # Adjusting the iteration value of k     ->     on going area of research. 
     # Setting incre = 2*incre is based on the deflation scheme used in the 
     # paper: Li, C., & Zhou, H. "svt: Singular value thresholding in MATLAB", 
     # Journal of statistical software, 81(2), (2017).
     #-------------------------------------------------------------------------
     k <- min(k+incre,psvdmax-Slen+S0len,kmax)
     work  <- max(3*k,20)
     incre <- 2*incre 
     iter <- iter+1
     maxit <- min(maxit+10,min_mn)
  } # end main while loop 
}   # end of svt_irlba
 # -------------------#
 # END: MAIN PROGRAM  #
 # -------------------#
 
# Routine to compute the block svd power method as described in: Bentbib, A.H.
# and  Kanber, A., "Block power method for SVD decomposition", Analele 
# ştiinţifice ale Universităţii Ovidius Constanţa. Seria Matematică, 
# 23(2), (2015) pp.45-58. 
# -----------------------------------------------------------------------------
 blocksvdpower <- function(A, x, maxiter,interchange) {
   if (interchange){
     qrstr <- qr(x)
     v <- qr.Q(qrstr) 
     for (j in 1:maxiter) {
       qrstr <- qr(A %*% v)
       u <- qr.Q(qrstr)
       qrstr <- qr(t(t(u) %*% A))
       R <- qr.R(qrstr)
       v <- qr.Q(qrstr)
     }
     svd_result <- svd(R)
     # Note: required change output order for svd of R.
     return(list(v = v %*%svd_result$u, d = svd_result$d, u = u %*% svd_result$v ))
   }
   qrstr <- qr(x)
   u <- qr.Q(qrstr) 
   for (j in 1:maxiter) {
      qrstr <- qr(t(t(u) %*% A))
          v <- qr.Q(qrstr)
      qrstr <- qr(A %*% v)
          R <- qr.R(qrstr)
          u <- qr.Q(qrstr)
   }
   svd_result <- svd(R)
   return(list(u = u %*%svd_result$u, d = svd_result$d, v = v %*% svd_result$v ))
 }
