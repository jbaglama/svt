#  A Hybrid Algorithm for Computing a Partial Singular Value Decomposition Satisfying a Given Threshold

---
#  Authors: 
-   James Baglama $\qquad\qquad\qquad$     &nbsp;      jbaglama@uri.edu
-   Jonathan Chávez-Casillas $\qquad$&nbsp; jchavezc@uri.edu
-   Vasilije Perovic  $\qquad\qquad\qquad$ &nbsp;&nbsp;      perovic@uri.edu
---

## Description:

Implements Algorithm 1 in reference [[1]](#1).
A hybrid function for computing all singular triplets above a given threshold. The function repeatedly calls a Partial Singular Value Decomposition (PSVD) method (svds or irlba) and a block SVD power method to compute a PSVD  of a matrix $A$. If a user threshold value `sigma` is specified, the function computes all the singular values above such threshold `sigma`. However, if an energy percentage ($\leq 1$) (`energy`) is specified, this function computes  the energy percentage. i.e.  an $r$-rank approximation $A_r$ of $A$ where  $(\Vert A_r\Vert_F/\Vert A\Vert_F)^2 \geq$ `energy`. If neither are given then the function returns a $k$-PSVD with `k=6` by default. 

## Input and Optional Parameters

The necessary inputs for the routine are:

| Input | Version |Description |
|:-------:|:---------:|------------|
|`A`| **Matlab/Octave** |An $m\times n$ numeric real matrix $A$ or a function handle (for Matlab only).|
|`A`| **${\tt R}$**| An $m\times n$ numeric sparse real matrix (matrix-product input is not available).|

Meanwhile, the optional arguments for the function in the distinct languages are:

| Parameters | Version | Description |
|:------------:|:-------------:|---------|
|`sigma`| all | Singular value threshold ($\geq 0$). If missing returns a $k$-PSVD.|
| or | |
|`energy`| all | Energy percentage (decimal $\leq 1$). Provides an $r$-rank approximation of $A$ <br>where  $(\Vert S_r\Vert_F/\Vert A\Vert_F)^2 \geq$ `energy`. It cannot be combined with a sigma value.  <br> **Matlab:** If specified, $A$ cannot be a function handle. (default: `[]`)<br> **${\tt R default:}$** `NULL`|
|`m`| **Matlab** | Number rows of $A$ - required if $A$ is a function handle (default: `[]`)|
|`n`| **Matlab** | Number columns of $A$ - required if $A$ is a function handle (default: `[]`)|
|`U0`| **Matlab/Octave** | Left singular vectors of an already computed PSVD of $A$ (default: `[]`)|
|`V0`| **Matlab/Octave** | Right singular vectors of an already computed PSVD of $A$ (default: `[]`)|
|`S0`| **Matlab/Octave** | Diagonal matrix  of an already computed PSVD of $A$ (default: `[]`)|
|`psvd`|  **${\tt R}$**  | ${\tt List}$ of an already computed PSVD of $A$. (default: `NULL`) <br> It should satisfy that `A %*% psvd$v = psvd$u %*% diag(psvd$d)` and <br> `t(A) %*% psvd$u = psvd$v %*% diag(psvd$d)`|
|`tol`| all | Tolerance used for convergence in the `svds`/`irlba` routine. <br> **Matlab default:** `sqrt(eps)` <br> **${\tt R}$ default:** `1e-8`|
|`k`| all | Initial number of singular triplets (default: `6`) |
|`incre`| all | Increment added to $k$ -internally doubled each iteration. (default: `5`) |
|`kmax`| all | Maximum value $k$ can reach. (default: `min(0.1*min(m,n),100)`)|
|`p0`| all | Starting vector for the `svds`/`irlba` routine. <br> **Matlab default:** `randn(max(n,m),1)` <br> **${\tt R}$ default:** `rnorm(n)`|
|`psvdmax`| all | Maximum dimension of the output PSVD. The output psvd will contain the  <br>input psvd if given. (default: `max(min(100+size(S0),min(n,m))),k)`) <br> <br> **NOTE:** This function will allocate memory for the full matrices $U$ and $V$.  It<br> might run out of memory  (return system error) on initialization for very <br>large $A$ and `psvdmax` value. |
|`pwrsvd`| all | If set to an integer $>0$, each iteration will perform `pwrsvd` iterations of a block<br> SVD power method with the output from  svds. If set to `0` then only one <br>iteration is performed if required (e.g. loss of orthogonality of basis vectors)<br>        (default: `0`)|
|`display`| all | If set to `1` displays iteration diagnostic information. (default: `0`)|
|`cmmf`| **${\tt R}$**| Logic variable, if `TRUE` use the internal custom matrix multiplication function<br> (`cmmf`)  to compute the matrix-product deflation. If `FALSE`, the routine will<br> use the `irlba` mult parameter - see irlba documentation for details.<br>  (default: `FALSE`)

## Output

Below, we describe the general output of the main routines.


### **Matlab/Octave Output:**

| Output | Description |
|:--------:|-------------|
| `U`| Left singular vectors.|
| `S`| Diagonal matrix of singular values sorted in decreasing order.|
| `V`| Right singular vectors. <br> **NOTE:** The output `U`,`S`,`V` will include `U0`,`S0`,`V0` if given.|
|`FLAG`| `0`: successful output - either the threshold (`sigma`) was met <br>  or the energy percentage was satisfied. <br><br>`1`:  the svds iteration failed to compute any singular triplets.<br> Outputs last values for `U`,`S`,`V`. <br><br> `2`:  `psvdmax` was reached before the threshold (`sigma`)<br>  was met or the energy percentage was satisfied.<br> Outputs the last values for `U`,`S`,`V`. <br><br> `3`:  no singular values above the specified threshold (`sigma`) exist. <br> Outputs `U=[]`,`V=[]`,`S=[]`.|


### **${\tt R}$ Output:**

| Output | Description |
|--------|-------------|
| `u`| Matrix of left singular vectors. |
| `d`| Vector of singular values in decreasing order.|
| `v`| Matrix of right singular vectors. <br> **NOTE:** The output `u`,`d`,`v` will include input values from psvd if given.|
|`FLAG`| `0`: successful output - either the threshold (`sigma`) was met <br>  or the energy percentage was satisfied <br><br>`1`:  if `irlba` iteration fails to compute any singular triplets. <br> Outputs last values for `u`,`d`,`v` <br><br> `2`:  `psvdmax` was reached before the threshold (`sigma`)<br>  was met or the energy percentage was satisfied.<br> Outputs the last values for `u`,`d`,`v` <br><br> `3`:  no singular values above the specified threshold (`sigma`) exist. <br> Outputs `u=NULL`,`d=NULL`,`v=NULL`.|

## Contents of the repository

| File name   |   Octave<br> Compatibility? |  Description  |
|:-----------:|:-----------------------------:|---------------|
|Matlab/svt_irlba.m | *Yes*| A matlab hybrid function that calls the internal MATLAB thick-restarted <br>GKLB routine `irlba`.|
|Matlab/svt_svds.m  | *No*| A matlab hybrid function that calls the internal MATLAB function `svds`.|
|Matlab/svt_interact_demo.m | *Yes*| This demo will allow the user to explore svt_svds or svt_irlba with varying <br>thresholds and options for 7 different matrices from the SuiteSparse<br> Matrix Collection. <br> **Note:** This requires to download the matrices from https://sparse.tamu.edu/|
|Matlab/svt_demo.m | *Yes*| This demo runs several examples with the two matrices `illc1033` and<br> `bibd_20_10` from the SuiteSparse Matrix Collection. <br> **Note:** This requires to download the matrices from https://sparse.tamu.edu/|
|Matlab/example41.m| *Yes*| Reproduces the Example 4.1 in [[1]](#1) |
|Matlab/example42.m| *No*| Reproduces the Example 4.2 in [[1]](#1) |
|R/svt_irlba.R| *NA*| An ${\tt R}$ hybrid function that calls the thick-restarted GKLB routine `irlba`. <br> This function can be though as the ${\tt R}$ version of `svt_irlba.m`.|
|R/svt_demo.R | *NA*| This demo runs several examples with the two matrices `illc1033` and <br> `bibd_20_10` from the SuiteSparse Matrix Collection.<br>  It is the ${\tt R}$ equivalent to `svt_demo.m`.<br> **Note:** This requires to download the matrices from https://sparse.tamu.edu/|
|R/example43.R| *NA*| Reproduces the Example 4.3 in [[1]](#1) |





##  EXAMPLES:

1. Compute all singular triplets with singular values exceeding `5.1`:

<div align="center">

| Software | Command |
|:--------:|-------------|
|**${\tt R}$:**| `psvd <- svt_irlba(A,sigma=5.1)`|
|**Matlab:** | `[U,S,V,FLAG] = svt_svds(A,'sigma',5.1);`|
|**Matlab:** | `[U,S,V,FLAG] = svt_irlba(A,'sigma',5.1);`|

</div>

2. Compute all singular triplets with singular values exceeding `5.1` given an initial PSVD `psvd0`:

<div align="center">

| Software | Command |
|:--------:|-------------|
|**${\tt R}$:**| `psvd <- svt_irlba(A,sigma=5.1,psvd=psvd0)`|
|**Matlab:** | `[U,S,V,FLAG] = svt_svds(A,'sigma',5.1,'U0',U0,'V0',V0,'S0',S0);`|
|**Matlab:** | `[U,S,V,FLAG] = svt_irlba(A,'sigma',5.1,'U0',U0,'V0',V0,'S0',S0);`|

</div>

3. Compute all singular triplets with singular values exceeding `1.1`, within tolerance of `1e-10` and provide a maximum of `20` singular triplets:

<div align="center">

| Software | Command |
|:--------:|-------------|
|**${\tt R}$:**| `psvd <- svt_irlba(A,sigma=1.1,tol=1e-10,psvdmax=20)`|
|**Matlab:** | `[U,S,V,FLAG] = svt_svds(A,'sigma',5.1,'tol',1d-10,'psvdmax',20);`|
|**Matlab:** | `[U,S,V,FLAG] = svt_irlba(A,'sigma',5.1,'tol',1d-10,'psvdmax',20);`|

</div>

4. Compute the top `6` singular triplets and then continue computing more. <br>Assume that based on the output, the desired threshold is 100<br> and set `psvdmax`, the number of singular triplets to `20`:

<div align="center">

| Software | Commands |
|:--------:|-------------|
|**${\tt R}$:**| `psvd0 <- svt_irlba(A)`  <br> `psvd <- svt_irlba(A,sigma=100,psvd=psvd0,psvdmax=20)`|
|**Matlab:** | `[U,S,V,FLAG] = svt_svds(A);`  <br> `[U,S,V,FLAG] = svt_svds(A,'sigma',100,'U0',U,'V0',V,'S0',S,'psvdmax',20);`|
|**Matlab:** | `[U,S,V,FLAG] = svt_irlba(A);` <br> `[U,S,V,FLAG] = svt_irlba(A,'sigma',100,'U0',U,'V0',V,'S0',S,'psvdmax',20);`|

</div>

**Note:** The user can check if any multiple singular values have been missed by calling the function again with the same parameters and threshold, but with the output from the previous call.

5. Compute the energy percentage 0.9:

<div align="center">

| Software | Command |
|:--------:|-------------|
|**${\tt R}$:**| `psvd <- svt_irlba(A,energy = 0.9)`|
|**Matlab:** | `[U,S,V,FLAG] = svt_svds(A,'energy',0.9);`|
|**Matlab:** | `[U,S,V,FLAG] = svt_irlba(A,'energy',0.9);`|

</div>

    
##  Date Last Modified: 
 5/1/24
  


## References
<a id="1">[1]</a> 
J. Baglama, J.Chavez-Casillas and V. Perovic, "A Hybrid Algorithm for Computing a Partial Singular Value Decomposition Satisfying a Given Threshold", submitted for publication 2024.

[2.]  J. Baglama and V. Perovic, "Explicit Deflation in Golub–Kahan–Lanczos Bidiagonalization Methods, ETNA, Vol. 58, (2023), pp. 164–176.

[3.] B.W. Lewis, J.Baglama, L. Reichel, "The irlba Package", (2021) https://cran.r-project.org/web/packages/irlba/

