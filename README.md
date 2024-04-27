#  **A Hybrid Algorithm for Computing a Partial Singular Value Decomposition Satisfying a Given Threshold**


## Description:

A hybrid singular value wrapper function that repeatedly calls the thick-restarted GKLB routine irlba and a block SVD power method to compute an approximate Partial Singular Value Decomposition (PSVD) of a matrix $A$. If a user threshold value (_sigma_) is given, the function will compute all singular values above such threshold (sigma). However, if an energy percentage ($\leq 1$) (_energy_) is given, the function computes the energy percentage. That is, it computes an $r$-rank approximation of $A$ where $(||S_r||_F/||A||_F)^2 >=$_energy_.  If neither are given then the function returns a $k$-PSVD. 


## Contents of the repository

| File name | Description  |
|-----------|--------------|
|svt_irlba.m | A matlab wrapper using the internal (MATLAB) thick-restarted GKLB routine ***irlba***.|
|svt_svds.m  | A matlab wrapper using the calls the MATLAB function ***svds***.|
|svt_interact_demo.m | This demo will allow the user to explore svt_svds or svt_irlba with varying thresholds<br> and options for 7 different matrices from the SuiteSparse Matrix Collection. <br> **Note:** This demo requires to download the matrices from https://sparse.tamu.edu/|
|svt_demo.m | This demo runs several examples with the two matrices `illc1033` and `bibd_20_10` from the SuiteSparse  <br> Matrix Collection. <br> **Note:** This demo requires to download the matrices from https://sparse.tamu.edu/|
|example41.m| 







#
# INPUT:
#    A - An m x n numeric sparse real matrix (matrix-product input is not available)
#    PARAMETERS:
#        sigma  - Singular value threshold (>= 0). If missing, the function returns a k-PSVD.
#         (or)
#       energy  - Energy percentage (decimal <= 1) - r-rank approximation of  
#                 A where (||psvd$d||_F/||A||_F)^2 >= energy. 
#                 Cannot be combined with a sigma value. (default: NULL)
#        psvd   - List of an already computed PSVD of A. (default: NULL)
#                     A %*% psvd$v = psvd$u %*% diag(psvd$d)
#                  t(A) %*% psvd$u = psvd$v %*% diag(psvd$d)
#        tol   -  Tolerance used for convergence in the irlba routine. (default: 1e-8)
#        k     -  Initial number of singular triplets. (default: 6)
#        incre -  Increment added to k. (internally doubled each iteration) (default: 5)
#        kmax  -  Maximum value k can reach. (default: min(.1*min(m,n),100))
#         p0   -  Starting vector for irlba. (default: p = rnorm(n))
#     psvdmax  -  Maximum dimension of the output ofPSVD psvd. The output 
#                 psvd will contain the input psvd if given.
#                 (default: max(min(100,min(n,m)-size(S0)),k)
#                 *NOTE* This function will allocate memory for the full matrices 
#                        and will run out of memory (return system error) on 
#                        initialization for very large A and psvdmax value. 
#      pwrsvd  -  If set to an integer > 0 then, on each iteration, it will perform 
#                 pwrsvd iterations of a block SVD power method with the
#                 output from irlba. If set to 0 only one iteration is performed when 
#                 required. (e.g. loss of orthogonality of basis vectors)(default: 0)
#     display  -  If set to 1 displays iteration diagnostic information. (default: 0)
#      cmmf   -   Logic variable, if TRUE use the internal custom matrix multiplication 
#                 function (cmmf) to compute the matrix-product deflation. If FALSE, will
#                 use the irlba mult parameter - see irlba documentation for details.  
#                 (default: FALSE)
#
#  OUTPUT:
#   Returns a list with entries:
#      u  -  Matrix of left singular vectors. 
#      d  -  Vector of singular values in decreasing order.
#      v  -  Matrix of right singular vectors.
#            *NOTE* - Output u,d,v will include input values from psvd if given.
#   flag  - 0 successful output - either threshold (sigma) met or energy 
#             percentage satisfied.  
#           1 if irlba iteration fails to compute any singular triplets - output 
#             last values for u,s,v.
#           2 if psvdmax is reached before threshold (sigma) is met or energy 
#             percentage is satisfied- output last values for u,s,v.
#           3 if no singvals above threshold (sigma) exist - output NULL for u,s,v.
#
#  EXAMPLES:
#
#      1. Compute all singular triplets with singular values exceeding 5.1:
#          > psvd <- svt_irlba(A,sigma=5.1)
#
#      2. Compute all singular triplets with singular values exceeding 5.1:
#         given an initial PSVD psvd0:
#          > psvd <- svt_irlba(A,sigma=5.1,psvd=psvd0)
#
#      3. Compute all singular triplets with singular values exceeding 1.1,
#         within tolerance of 1e-10 and provide a maximum of 20 singular triplets:
#          > psvd <- svt_irlba(A,sigma=1.1,tol=1e-10,psvdmax=20)
#
#      4. Compute the top 6 singular triplets and then continue:
#          > psvd0 <- svt_irlba(A)  
#         Assume that based on the output the desired threshold is 100 and set
#         psvdmax, the number of singular triplets to 20: 
#          > psvd <- svt_irlba(A,sigma=100,psvd=psvd0,psvdmax=20)
#         The user can check if any multiple singular values have been missed by
#         calling the function again with the same parameters and threshold, but
#         with the output from the previous call
#
#      5. Compute the energy percentage 0.9:
#          > psvd <- svt_irlba(A,energy = 0.9)
#    
#  DATE LAST MODIFIED: 
#  4/24/24
#  
#  AUTHORS: 
#   James Baglama            email: jbaglama@uri.edu
#   Jonathan Chavez-Casillas email: jchavezc@uri.edu
#   Vasilije Perovic         email: perovic@uri.edu
#
## REFERENCES:
[1.] J. Baglama, J.Chavez-Casillas and V. Perovic, "A Hybrid Algorithm for Computing a Partial Singular Value Decomposition Satisfying a Given Threshold", submitted for publication 2024.

[2.] J. Baglama and V. Perovic, "Explicit Deflation in Golub–Kahan–Lanczos Bidiagonalization Methods, ETNA, Vol. 58, (2023), pp. 164–176.

[3.] B.W. Lewis, J.Baglama, L. Reichel, "The irlba Package", (2021) https://cran.r-project.org/web/packages/irlba/

