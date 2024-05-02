function [U,S,V,FLAG] = svt_svds(A,varargin)
% 
%  DESCRIPTION:
%  Implements Algorithm 1 in reference 1. 
%  A hybrid function for computing all singular triplets above a given threshold. 
%  The function repeatedly calls a Partial Singular Value Decomposition (PSVD)
%  method (svds or irlba) and a block SVD power method to compute a PSVD 
%  of a matrix A. If a user threshold value (sigma) is specified, the function  
%  computes all the singular values above such threshold (sigma). However, if an 
%  energy percentage (<=1) (energy) is specified, this function computes 
%  the energy percentage. i.e.  an r-rank approximation A_r of A where 
%  (||A_r||_F/||A||_F)^2 >= energy. If neither are given then the function 
%  returns a k-PSVD with k=6 by default.   
%
%  INPUT:
%  A - An m x n numeric real matrix A or a function handle.
%
%   PARAMETERS(optional):
%   Comma-separated pairs consisting of the parameter name in single 
%   quotations and its value.
%
%    'sigma' - singular value threshold (>= 0) if missing returns a k-PSVD.
%     (or)
%   'energy' - energy percentage (decimal <= 1) - r-rank approximation A_r of 
%              A where (||A_r||_F/||A||_F)^2 >= energy.  If specified, A cannot
%              be a function handle and cannot be combined with a sigma 
%              value. (default: [])
%     'm'    - # rows of A - required if A is a function handle 
%              (default: []) 
%     'n'    - # cols of A - required if A is a function handle 
%              (default: []) 
%     'U0'   - left singular vectors of an already computed PSVD of A 
%              (default: [])
%     'V0'   - right singular vectors of an already computed PSVD of A
%              (default: [])
%     'S0'   - diagonal matrix of an already computed PSVD of A 
%              (default: [])
%    'tol'   - tolerance used for establishing convergence in the svds 
%              routine (default: sqrt(eps))
%     'k'    - initial number of singular triplets (default: 6)
%  'incre'   - increment added to k (internally doubled each iteration) 
%              (default: 5)
%   'kmax'   - maximum value k can reach (default: min(0.1*min(m,n),100))
%     'p0'   - starting vector for svds (default: p = randn(max(n,m),1))
%  'psvdmax' - maximum dimension of the output PSVD [U,S,V] - output 
%              [U,S,V] will contain the input [U0,V0,S0] if given 
%              (default: max(min(100+size(S0),min(n,m)),k))
%              *NOTE* This function will allocate memory for the full 
%                     matrices U and V. It might run out of memory 
%                     (return system error) on initialization for very 
%                     large A and psvdmax value. 
%  'pwrsvd' -  if set to an integer > 0, each iteration will perform pwrsvd 
%              iterations of a block SVD power method with the output from 
%              svds. If set to 0 then one iteration is performed if 
%              required (e.g. loss of orthogonality of basis  vectors)
%              (default: 0)
% 'display' -  if set to 1 displays some diagnostic information per 
%              iteration (default: 0)
%
%  OUTPUT:
%      U - left singular vectors 
%      S - diagonal matrix of singular values sorted in decreasing order
%      V - right singular vectors
%          *NOTE* - The output U,S,V will include U0,S0,V0 if given
%   FLAG - 0  successful output - either the threshold (sigma) was met or  
%             the energy percentage was satisfied  
%          1  the svds iteration failed to compute any singular triplets - 
%             it outputs last values for U,S,V
%          2  psvdmax was reached before the threshold (sigma) was met or 
%             the energy percentage was satisfied- it outputs the last
%             values for U,S,V
%          3  no singular values above the specified threshold (sigma) 
%             exist - it outputs U=[],V=[],S=[]
%
%  EXAMPLES:
%      1. Compute all singular triplets with singular values exceeding 5.1
%           >> [U,S,V,FLAG] = svt_svds(A,'sigma',5.1);
%
%      2. Compute all singular triplets with singular values exceeding 5.1
%         given an initial PSVD AV0 = U0S0 and A'U0 = V0S0
%           >> [U,S,V,FLAG] = svt_svds(A,'sigma',5.1,'U0',U0,'V0',V0,'S0',S0);
%
%      3. Compute all singular triplets with singular values exceeding 5.1
%         with tolerance of 1d-10 and 20 maximum singular triplets 
%           >> [U,S,V,FLAG] = svt_svds(A,'sigma',5.1,'tol',1d-10,'psvdmax',20);
%
%      4. Compute the top 6 singular triplets and then continue
%           >> [U,S,V,FLAG] = svt_svds(A);  
%         Assume based on output the desired threshold is 100 and set
%         psvdmax, the max. number of singular triplets to 20 
%           >> [U,S,V,FLAG] = svt_svds(A,'sigma',100,'U0',U,'V0',V,'S0',S,'psvdmax',20);
%         the user can check if any multiple singular values have been missed by
%         calling the function again with the same parameters and threshold, but
%         with the output from the previous call
%
%      5. Compute the PSVD with an energy percentage of 0.9
%           >> [U,S,V,FLAG] = svt_svds(A,'energy',0.9);
%    
%  DATE LAST MODIFIED: 
%  5/1/24
%
%  LANGUAGE:
%  MATLAB versions: R2018b ... R2023b (earlier releases will not work)
%  
%  AUTHORS: 
%  James Baglama            email: jbaglama@uri.edu
%  Jonathan Chávez-Casillas email: jchavezc@uri.edu
%  Vasilije Perovic         email: perovic@uri.edu
%
% REFERENCES:
%  1. J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm for 
%     Computing a Partial Singular Value Decomposition Satisfying a Given 
%     Threshold", submitted for publication 2024.
%  2. J. Baglama and V. Perovic, "Explicit Deflation in Golub–Kahan–Lanczos 
%     Bidiagonalization Methods, ETNA, Vol. 58, (2023), pp. 164–176.
% --------------------------------------------------------------------------

% Setting the clock time for display option.
% ------------------------------------------
tStart = tic; 

% Parsing the optional input parameters and set the default values.
% -----------------------------------------------------------------
p = inputParser; % Default setting is not case sensitive.
addParameter(p,'sigma',Inf);
addParameter(p,'energy',[]);
addParameter(p,'U0',[]);
addParameter(p,'V0',[]);
addParameter(p,'S0',[]);
addParameter(p,'m',-1);
addParameter(p,'n',-1);
addParameter(p,'tol',sqrt(eps));
addParameter(p,'psvdmax',[]);
addParameter(p,'p0',[]);
addParameter(p,'kmax',[]);
addParameter(p,'pwrsvd',0);
addParameter(p,'display',0);
addParameter(p,'incre',5);
addParameter(p,'k',6);
parse(p,varargin{:});
p = p.Results;

% Checking the input matrix A
% ----------------------
if isempty(A), error('Missing input matrix A.'); end
if isnumeric(A)
   [m,n] = size(A);
else
   m = p.m; n = p.n;
   if ~isscalar(m) || ~isscalar(n) || m ~=floor(m) || n ~=floor(n)
       error('Incorrect values for m and/or n');
   end
   if (m < 0 || n < 0) 
      error('A is a function handle - missing m and/or n');
   end
end
min_mn = min(m,n); max_mn = max(m,n);

% Checking for errors of input parameters and setting some internal variables.
% ----------------------------------------------------------------------------
kmax = p.kmax; display = p.display; k = p.k; energy = p.energy;
incre = p.incre; sigma = p.sigma; pwrsvd =  p.pwrsvd; p0 = p.p0;  
if isempty(kmax), kmax = max(floor(min(0.1*min_mn,100)),k); end 
if ~isscalar(kmax) || kmax ~= floor(kmax) || kmax < 0 || kmax > min_mn
   error('Incorrect value for kmax.');
end
if ~isscalar(k) || k <= 0 || k >= min_mn || k ~= floor(k)
   error('Incorrect value for k.');
end 
k_org = k;
if ~isscalar(incre) || incre <= 0 || incre ~= floor(incre)
   error('Incorrect value for incre.');
end
if ~isscalar(pwrsvd) || pwrsvd < 0 || pwrsvd~=floor(pwrsvd)
   error('Incorrect value for pwrsvd.');
end
if ~isscalar(display) || (display ~= 0 && display ~= 1)
   error('Incorrect value for display.');
end
if ~isscalar(sigma) || sigma < 0
    error('Incorrect value for sigma.');
end
if ~isempty(energy)
    if ~isscalar(energy) || energy <= 0 || energy > 1
       error('Incorrect value for energy - 0 < energy <= 1.');
    end
    if sigma ~= Inf
       error('Cannot have a value for sigma and energy.');
    end
    if ~isnumeric(A)
       error('A must be a numeric matrix.')
    end
    FronormA = norm(A,'fro')^2;
end
if isempty(p0), p0 = randn(max_mn,1); end
if (m <=n && size(p0,1) ~= n), error('m <= n p0 must be n x 1'); end
if (m > n && size(p0,1) ~= m), error('m > n  p0 must be m x 1'); end

% Setting up A to always be called as a function handle. The numeric 
% matrix A is now included (stored) inside the function handle.
% ------------------------------------------------------------------
if ~isa(A,'function_handle'), A = @(x,t) Afun(A,x,t); end

% Checking the tolerance and setting the internal tolerance values.
% -----------------------------------------------------------------
tol = p.tol; if tol < eps, tol = eps; end 
eps12 = eps^(1/2); Smax = 1; 

% Checking the sizes of the input matrices U0, V0, S0 - this assumes A is an 
% m x n matrix and that the PSVD has the form: A*V0 = U0*S0 and A'*U0=V0*S0.
% --------------------------------------------------------------------------
[mU,kU] = size(p.U0); [nV,kV] = size(p.V0); [k1S,k2S] = size(p.S0); 
Slen = k1S;  
if any([mU kU nV kV k1S k2S]) % checking to see if any input (sizes) are not 0.
   if nV ~= n || mU ~= m || kV ~= kU || k1S ~= k2S || kV~=k1S || kU~=k1S
       error(' Incorrect dimensions of U0, S0, V0.');
   end
   p.S0 = diag(p.S0); 
   if ~issorted(p.S0,'descend')
      [~,I] = sort(p.S0,'descend'); 
      p.U0 = p.U0(:,I); p.V0 = p.V0(:,I); p.S0 = p.S0(I); 
   end
   Smax = p.S0(1); Slen = length(p.S0);
   if p.S0(Slen) < sigma && sigma ~=Inf
       FLAG = 0; V=p.V0; U=p.U0; S=diag(p.S0); return;
   end
end

% Setting the psvdmax value.
% -------------------
psvdmax = p.psvdmax; S0len = Slen;
if isempty(psvdmax)
    psvdmax = min(100,min_mn-Slen); 
else
    psvdmax = min(psvdmax,min_mn-Slen);
end
if k >= min_mn - Slen, error('Incorrect initial k'); end
psvdmax = max(psvdmax,k);
if psvdmax <= 0  
   warning('Input dim of U,S,V >= min(m,n)');
   if ~isempty(p.S0) 
       S = diag(p.S0); V=p.V0; U=p.U0; FLAG = 2;
   else
       S=[]; U=[]; V=[]; FLAG = 3;
   end 
   return;
end
 
% Allocating memory for V, U, and S. 
% NOTE: If psvdmax is too large an out of memory error will occur here.
% ---------------------------------------------------------------------
V(n,psvdmax) = 0; U(m,psvdmax) = 0; S(psvdmax,1) = 0;
if Slen > 0
    V(:,1:Slen) = p.V0;
    U(:,1:Slen) = p.U0;
    S(1:Slen)   = p.S0;
end

% Removing the p. and input varargin from memory.
% -----------------------------------------------
clear varargin p;

%---------------------%
% BEGIN: MAIN PROGRAM %
%---------------------%

% Initializing values for the first call to svds.
% -----------------------------------------------
maxit = 100;        % Initial max # of iterations for svds.
sdim = max(3*k,20); % Initial max dimension of subspace for svds.
iter = 1;           % Overall iteration count.

% Setting up the starting vector and checking if power SVD iteration 
% should be performed with an input.
% -------------------------------------------------------------------------
% Given that A is an m x n matrix with m <= n => the "short vectors" are U.
% Working with the matrix: (I - UU')A
% (I - UU')A(PV_B) = (QU_B)S_B
% A'(I-UU')(QU_B) = (PV_B)S_B + fe'U_B
% Starting vector is p0 => RightStartVector "long vector" (n x 1).
if m <= n
   startvec = 'RightStartVector';
   if Slen > 0 
       if pwrsvd
          [U1,S1,V1] = blocksvdpower(A,U(:,1:Slen),pwrsvd,startvec);
          Slen = size(S1,1); S0len = Slen; 
          U(:,1:Slen) = U1; S(1:Slen) = diag(S1);V(:,1:Slen) = V1; 
       end
      p0 = p0 - V(:,1:Slen)*(V(:,1:Slen)'*p0); 
   end
% Given that A is an m x n matrix with m > n => "short vectors" are V
% Working with the matrix: (I-VV')A'
% (I-VV')A'(QU_B) = (PV_B)S_B
% A'(I-VV')(PV_B) = (QU_B)S_B + fe'V_B
% Starting vector is p0 => LeftStartVector "long vector" (m x 1)
else
    startvec = 'LeftStartVector';
    if Slen > 0
        if pwrsvd
          [U1,S1,V1] = blocksvdpower(A,V(:,1:Slen),pwrsvd,startvec);
          Slen = size(S1,1); S0len = Slen; 
          U(:,1:Slen) = U1; S(1:Slen) = diag(S1);V(:,1:Slen) = V1;
        end
      p0 = p0 - U(:,1:Slen)*(U(:,1:Slen)'*p0); 
   end
end
p0 = p0/norm(p0);

% This is the main while loop.
% ----------------------------
while iter < min_mn

   % Resetting internal values for each iteration.
   % ---------------------------------------------
   svd_fail = 0; FLAG = 1; S1 = []; svd_iter = 1; orth_chk = 0;
   tStart_svds = tic; blksvdpw = 0; blksvdpw_time = 0; 

   % Call to the MATLAB function svds to get the singular triplets.
   % Loop for a maximum of 2 iterations until S1 is not empty. Several 
   % numerical examples do not suggest that more iterations are needed.
   % ------------------------------------------------------------------
   while FLAG && isempty(S1) && svd_iter <= 2
       if Slen == 0
           [U1,S1,V1,FLAG] = svds(@(x,t)A(x,t),[m n],k,'largest',...
           'Tolerance',tol,startvec,p0,'FailureTreatment','drop',...
           'MaxIterations',svd_iter*maxit,'SubspaceDimension',svd_iter*sdim);
        elseif strcmp(startvec,'RightStartVector') 
           [U1,S1,V1,FLAG] = svds(@(x,t)AfunU(@(x,t)A(x,t),x,t,U(:,1:Slen)),...
           [m n],k,'largest','Tolerance',tol,startvec,p0,'FailureTreatment',...
           'drop','MaxIterations',svd_iter*maxit,'SubspaceDimension',svd_iter*sdim);   
        elseif strcmp(startvec,'LeftStartVector')
           [U1,S1,V1,FLAG] = svds(@(x,t)AfunV(@(x,t)A(x,t),x,t,V(:,1:Slen)),...
           [m n],k,'largest','Tolerance',tol,startvec,p0,'FailureTreatment',...
           'drop','MaxIterations',svd_iter*maxit,'SubspaceDimension',svd_iter*sdim);   
       end
       svd_iter = svd_iter+1;
   end
   svds_time = toc(tStart_svds);
   
   % Checking if svds did not compute all k values. svds option 
   % 'FailureTreatment' is set so that the values where it did not 
   % converged are dropped.
   % -----------------------------------------------------------------
   if FLAG || svd_iter > 2
      if isempty(S1) % No computed values from svds - exiting program.
         % resizing U, V, S for output.
         S = S(1:Slen); V = V(:,1:Slen); U = U(:,1:Slen); 
         if ~isempty(S)
            [S,I] = sort(S,'descend'); U = U(:,I); V = V(:,I); S=diag(S);
         end
         if display
            fprintf('MATLAB svds failed to return any singular values.\n');
         end
         return; 
      end 
      svd_fail = 1;  % svds returned some values - continuing.
      maxit = 2*maxit;
   end
   
   % Resetting the k value based on the returned matrices from svds.
   % ---------------------------------------------------------------
   k = size(S1,1); comput_k=k; Slenk = Slen+k; Skindex = (Slen+1):Slenk;

   % Checking orthogonality among the "long" vectors - That is, 
   % the vectors not explicitly orthogonalized against basis vectors.
   % orth_chk is based on Theorem 5 (page 31)
   % Larsen, Rasmus Munk. "Lanczos bidiagonalization with partial 
   % reorthogonalization." DAIMI Report Series 537 (1998).
   % ----------------------------------------------------------------
   if ~isempty(S1) && Slen > 0 && ~svd_fail && pwrsvd == 0
      if strcmp(startvec,'RightStartVector')
         orth_chk = abs(max(V1'*V(:,1:Slen),[],'all','Comparisonmethod',...
                              'abs'));
      else 
         orth_chk = abs(max(U1'*U(:,1:Slen),[],'all','Comparisonmethod',...
                                'abs'));
      end
   end

   % Calling block svd power method. This improves the orthogonality of the 
   % basis vectors and the overall error. However, this can be 
   % computationally expensive, especially if A is large and/or the number 
   % of computed singular triplets is large.
   % If pwrsvd = 0, determination is set by orth_chk.
   % ----------------------------------------------------------------------
   if orth_chk > eps12/sqrt(Slen+k) || pwrsvd || ...
      svd_fail || S1(k,k) < Smax*eps12
       tStart_blk = tic;
       if strcmp(startvec,'RightStartVector') 
          [U1,S1,V1] = blocksvdpower(A,...
          [U(:,1:Slen) U1],max(1,pwrsvd),startvec);
       else
         [U1,S1,V1] = blocksvdpower(A,...
         [V(:,1:Slen) V1],max(1,pwrsvd),startvec);
       end
       Slenk = size(S1,1); Skindex = 1:Slenk;
       blksvdpw = 1; blksvdpw_time = toc(tStart_blk);
   end
  
   % Updating matrices U, V, S, and starting vector p0. 
   % --------------------------------------------------
   U(:,Skindex) = U1; V(:,Skindex) = V1; S(Skindex) = diag(S1);
   if strcmp(startvec,'RightStartVector')
      p0 = p0 - V1*(V1'*p0);
   else
      p0 = p0 - U1*(U1'*p0);
   end
   p0 = p0/norm(p0);

   % Updating the number of computed singular triplets.
   % --------------------------------------------------
   Slen = Slenk;
   
   % Checking convergence.
   % ---------------------
   Smin = min(S(1:Slen));
   Smax = max(max(S(1:Slen)),Smax);
   statement1 = 0; statement2 = 0;
   if isempty(energy)
      statement1 = Smin < sigma;
   else
      statement2 = (norm(S(1:Slen))^2/FronormA) >= energy;
   end
   statement3 = Slen >= (psvdmax+S0len);
   statement4 = Slen >= min_mn;
   if display
      fprintf(' \n');
      fprintf(' SVDS: CPU time = %0.5g secs. \n',svds_time);
      fprintf(' SVDS: computed singvals  : %d   \n',comput_k);
      fprintf(' Current size of PSVD of A: %d \n', Slen);
      fprintf(' Approx. max. singval of A:  %0.5g \n',Smax);
      fprintf(' Approx. min. singval of A: %0.5g \n', Smin);
      if blksvdpw
         fprintf([' Block SVD power method: iterations = %d, ' ... 
                  'CPU time = %0.5g secs. \n'],max(1,pwrsvd),blksvdpw_time);  
      end
   end
   if statement1 || statement2 || statement3 || statement4
     % Convergence or max number of singular triplets reached.
     % Resizing U, V, S for output.
     S = S(1:Slen); V = V(:,1:Slen); U = U(:,1:Slen); FLAG = 0;
     if ~issorted(S,'descend')
        [S,I] = sort(S,'descend'); U = U(:,I); V = V(:,I);
     end
     Slen1 = min([Slen,min_mn,psvdmax]);
     if Slen1 < Slen
        S = S(1:Slen1); V = V(:,1:Slen1); U = U(:,1:Slen1); 
     end   
     % Checking if sigma is set to a threshold value.
     if sigma ~= Inf 
        % Already ordered largest to smallest.  
        I = find(S >= sigma); U = U(:,I); V = V(:,I); S = S(I); 
     end
     if statement2
        for i = Slen:-1:1
            FronormS = norm(S(1:i))^2;
            if FronormS/FronormA < energy
               S = S(1:i+1); V = V(:,1:i+1); U = U(:,1:i+1); break;
            end
        end
     end
     if ~statement1 && ~statement2 && statement3 && ~statement4, FLAG = 2; end
     S = diag(S);
     if isempty(S), FLAG = 3; end
     if display
         if FLAG == 0
           fprintf(' Successful return \n');
         end
        fprintf(' Total time for function svt_svds: %0.5g secs. \n',toc(tStart));
     end
     return;
   end
   
  % Adjusting iteration value of k      ->      on going area of research. 
  % Setting incre = 2*incre is based on the delfation scheme used in the 
  % paper: Li, C., & Zhou, H. "svt: Singular value thresholding in MATLAB", 
  % Journal of statistical software, 81(2), (2017).
  % -----------------------------------------------------------------------
  k = min([k+incre,psvdmax-Slen+S0len,kmax]); incre = 2*incre; iter=iter+1;
  sdim = max(3*k,20); maxit = min(maxit+10,min_mn);

end % end main while loop.

end % end for function svt_svds.

%-------------------%
% END: MAIN PROGRAM %
%-------------------%

% Routine to compute the block svd power method as described in: 
% Bentbib, A.H. and Kanber, A., "Block power method for SVD decomposition", 
% Analele ştiinţifice ale Universităţii Ovidius Constanţa. Seria Matematică, 
% 23(2), (2015) pp.45-58. This improves the orthogonality of the basis vectors 
% and the overall error. Typically, 1 or 2 iterations is sufficient.
% -------------------------------------------------------------------------
function [U,S,V] = blocksvdpower(A,X,maxiter,startvec)
  if strcmp(startvec,'RightStartVector') 
     [U,~]  = qr(X,0);
     for j=1:maxiter
         [V,~]  = qr(A(U,'transp'),0);
         [U,R]  = qr(A(V,'notransp'),0);
     end
     [u_r,S,v_r] = svd(R); 
     V = V*v_r; U = U*u_r;
  else
     [V,~]  = qr(X,0);
     for j=1:maxiter
         [U,~]  = qr(A(V,'notransp'),0); 
         [V,R]  = qr(A(U,'transp'),0);
     end
     % Note: required change in output order for svd of R.
     [v_r,S,u_r] = svd(R); 
     U = U*u_r; V = V*v_r;
  end
end

% Function to use in case a numeric matrix A is given - 
% otherwise function is not used
% -----------------------------------------------------
function y = Afun(A,x,t)
  if strcmp(t,'notransp'), y = A*x; else, y = A'*x; end
end

% Given that A is an m x n matrix with m <= n, then the "short vectors" are 
% U (m x k). Work with the matrix: (I - UU')A. Thus, two-sided
% (I - UU')A(PV_B) = (QU_B)S_B
% A'(I-UU')(QU_B)  = (PV_B)S_B + fe'U_B
% ----------------------------------------------------------------------------
function y = AfunU(A,x,t,U)
if strcmp(t,'notransp') % y = A*x
   % Working with the matrix: (I - UU')A
   y = A(x,t);  y = y - U*(U'*y);
else % y = A'*x
     % Working with the matrix: ((I - UU')A)' = A'(I-UU')
   x = x - U*(U'*x); y = A(x,t);
end
end

% Given A m x n with m > n, then the "short vectors" are V
% Working with the matrix: (I-VV')A'. Thus, two-sided
% (I-VV')A'(QU_B) = (PV_B)S_B
% A'(I-VV')(PV_B) = (QU_B)S_B + fe'V_B
% ----------------------------------------------------------------
function y = AfunV(A,x,t,V)
if strcmp(t,'notransp') % y = A*x
   % Working with the matrix: ((I-VV')A')' = A(I-VV')
   x = x - V*(V'*x); 
   y = A(x,t);
else % y = A'*x
     % Working with the matrix: (I-VV')A'
      y = A(x,t); y = y - V*(V'*y); 
end
end