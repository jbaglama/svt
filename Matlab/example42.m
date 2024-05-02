function example42
% 
%  PAPER:
%  J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm for 
%  Computing a Partial Singular Value Decomposition Satisfying a Given 
%  Threshold", submitted for publication 2024.
% 
%  EXAMPLE 4.2 FROM PAPER:
%  Implementation of the Singular Value Threshold Matrix Completion (SVT-MC)
%  algorithm (Algorithm 1 page 1974) as described in the SVT paper 
%  1. J-F. Cai, E.J. Candès and Z. Shen, "A singular value thresholding algorithm 
%     for matrix completion", SIAM J. Optim Vol. 20, No. 4, 2010 pp. 1956-1982.
%  The MATLAB code presented here follows the implementation of SVT.m written 
%  by E.J. Candès:
%  https://www.convexoptimization.com/wikimization/index.php/Matrix_Completion.m
%  A very similar version by Stephen Becker is posted on 
%  https://candes.su.domains/software/svt/code/.  
%  This example follows the exact SVT-MC algorithm outlined in the proposed
%  paper and is only modified to call PSVD methods: svt_svds, svt_irlba, svt,
%  and svds in place of PROPACK. CPU times are also recorded for comparisson
%  purposes.
%
%  REQUIRED SOFTWARE:  
%  1. svt_irlba.m  and svt_svds.m  
%     https://github.com/jbaglama/svt/
%  2. svt.m
%     https://github.com/Hua-Zhou/svt/%  
% 
%  LANGUAGE:
%  MATLAB versions: R2018b ... R2023b (earlier releases will not work)
%
%  DATE LAST MODIFIED:  
%  5/1/24
%  
%  AUTHORS: 
%  James Baglama            email: jbaglama@uri.edu
%  Jonathan Chávez-Casillas email: jchavezc@uri.edu
%  Vasilije Perovic         email: perovic@uri.edu
% ---------------------------------------------------

% This example is only set-up to run in MATLAB.  
% ---------------------------------------------
uiIsOctave = uioctave; 
if uiIsOctave
    error('Octave environment - requires MATLAB platform\n');
end

% Checking if svt.m., svt_svds.m and svt_irlba.m exists.
% ------------------------------------------------------
if ~exist('svt_svds.m','file'), error('svt_svds.m requried for Example42'); end 
if ~exist('svt_irlba.m','file'), error('svt_irlba.m requried for Example42'); end
if ~exist('svt.m','file'), error('svt.m requried for Example42'); end

% Use this to determine which PSVD method to use:
% method = 1 -> svt_svds
% method = 2 -> svt_irlba
% method = 3 -> svt
% method = 4 -> svds (MATLAB's internal function).
% ------------------------------------------------
method = 1;

% PART 1: Setting up the matrix M, P_Omega, and P_Omega(M).
% ---------------------------------------------------------

% Determining which rank to use.
% 1. rank = 20  = 1%
% 2. rank = 50  = 2.5%
% 3. rank = 100 = 5%
% 4. rank = 200 = 10%
% 5. rank = 300 = 15%
r = 50;

% M - matrix: reserving the size and rank of the matrix. 
%             Also, fixing the random state for the matrix. 
% ---------------------------------------------------------
n1 = 2000; n2 = 20000;
rng(2009); M = randn(n1,r)*randn(r,n2);

% Storing the initial clock time.
% -------------------------------
tStart = tic; 

% NOTE: any modification of n1,n2, and r (or a different type of matrix M) may 
% require changes to the parameters for a successful output of SVT-MC algorithm.
% ------------------------------------------------------------------------------

% P_Omega - matrix: random sampling and permutation - ref [16] in SVT paper 
% --------------------------------------------------------------------------
df = r*(n1+n2-r);                           % Degrees of freedom - ref [16] in SVT paper.
oversampling = 4;                           % ref[16] in SVT sampling slightly more than df.
m = min(oversampling*df,round(.99*n1*n2));  % Number of samplings.   
perob = round(m/(n1*n2)*100);               % Percentage of entries that are observed.
Omega = randperm(n1*n2,m);                  % m random permutation of indices n1*n2.
data = M(Omega);                            % Data is a  1 x m vector of random entries of M.
if ~isempty(find(data == 0))                % If data has zero entries - it needs to add some.
   data = data + eps*randn(size(data));     % Adding noise to data to avoid an error when using
end                                         % MATLAB's sparse with data - indices will not match up later on
norm_data = norm(data);                     % Used for relative stopping criteria
[~,indx] = sort(Omega);                     % Column-major ordering to help with reording speed
[indx_i,indx_j] = ind2sub([n1,n2],Omega);   % Linear index of Omega to multiple index [i,j]
ProjM = sparse(indx_i,indx_j,data,n1,n2,m); % ProjM = P_Omega(M)
normProjM = normest(ProjM,1e-2);            % ProjM can be a very large matrix - can be replace with svds

% Parameters - internal parameters set from SVT paper.
% ----------------------------------------------------
p  = m/(n1*n2);                   % Following Equation (5.2), page 1971. 
tau = round(5*sqrt(n1*n2));       % tau page 1973 (original for square matrices).
delta = 1.2/p;                    % Following Equation (5.1), page 1971. 
maxiter = 3000;                   % Max. # of iteration k_max Alg. 1 page 1974.  
tol = 1e-3;                       % Stopping tol. \epsilon ALg. 1 page 1974. 
incre = 5;                        % Increment size \ell page 1971 - suggested to use 5.
k0 = ceil(tau/(delta*normProjM)); % Initializing k0 from Equation (5.3) page 1972.
r0 = r;                           % Setting up the original rank for output of results.

% Parameters for PSVD methods.
% ----------------------------
psvd_tol = sqrt(eps);
rng(2023); p0 = randn(max(n1,n2),1);
if n1 <= n2
    startvec = 'RightStartVector';
else
    startvec = 'LeftStartVector';
end

% PART 2:  SVT-MC Algorithm. 
% --------------------------

% Step 1: SVT-MC Algorithm. 
% -------------------------
y = k0*delta*data;
Y = k0*delta*ProjM;
[indx_i,indx_j,~] = find(Y); 

% Step 2: SVT-MC Algorithm. 
% -------------------------
r = 0;

% Initializing the timing value.
% ------------------------------
psvdtime = 0; 

% Step 3: SVT-MC Algorithm.
% -------------------------
for k = 1:maxiter 

    % Step 4: SVT-MC Algorithm.
    % -------------------------
    s = min([r + 1, n1, n2]);     
    
    % Steps 5-9: SVT-MC Algorithm.
    % ----------------------------
    if method == 1  % svt_svds method.               
       tpsvd = tic;
       [U,Sigma,V,FLAG] = svt_svds(Y,'sigma',tau,'psvdmax',max(n1,n2),...
                                     'tol', psvd_tol,'k',s,'p0',p0);
       r = size(Sigma,1); sigma = diag(Sigma);  sigma = sigma - tau; 
       Sigma = diag(sigma); psvdtime = psvdtime + toc(tpsvd);
       if FLAG, fprintf('FLAG = %d\n',FLAG); end
       if FLAG && isempty(Sigma), error('svt_svds failed to return any SVs'); end
    elseif method == 2 % svt_irlba method

       tpsvd = tic; 
       [U,Sigma,V,FLAG] = svt_irlba(Y,'sigma',tau,'psvdmax',max(n1,n2),...
                                      'tol',psvd_tol,'k',s,'p0',p0);
       r = size(Sigma,1); sigma = diag(Sigma);  sigma = sigma - tau; 
       Sigma = diag(sigma);psvdtime = psvdtime + toc(tpsvd);
       if FLAG, fprintf('FLAG = %d\n',FLAG); end
       if FLAG && isempty(Sigma), error('svt_irlba failed to return any SVs'); end
    elseif method == 3 % svt method 
       tpsvd = tic;
       [U,Sigma,V,FLAG] = svt(Y,'lambda',tau,'tol',psvd_tol,'k',s);
       r = size(Sigma,1); sigma = diag(Sigma);  sigma = sigma - tau; 
       Sigma = diag(sigma);psvdtime = psvdtime + toc(tpsvd);
       if FLAG, fprintf('FLAG = %d\n',FLAG); end
       if FLAG && isempty(Sigma), error('svt failed to return any SVs'); end
    elseif method == 4 % svds method - original setup scheme used PROPACK 
        [mY,nY] = size(Y); tpsvd = tic; 
        
        % Step 5: SVT-MC Algorithm.
        % -------------------------
        while 1
          
          % Step 6: SVT-MC Algorithm.
          % -------------------------
          % Calling function handle to avoid svds from calling full svd function
          [U,Sigma,V,FLAG] = svds(@(x,t)Yfun(x,t,Y),[mY nY],s,'largest',...
                                  'FailureTreatment','drop','Tolerance',...
                                   psvd_tol,startvec,p0);
          if FLAG && isempty(Sigma), error('svds failed to return any SVs'); end
          if FLAG && ~isempty(Sigma), s = size(Sigma,1); end

          % Step 8: SVT-MC Algorithm.
          % -------------------------
          if (Sigma(s,s) <= tau) || ( s == min(n1,n2) ), break; end

          % Step 7: SVT-MC Algorithm.
          % -------------------------
          s = min(s + incre, min(n1,n2));
         
        end 

        % Step 9: SVT-MC Algorithm.
        % -------------------------
        sigma = diag(Sigma); r = sum(sigma > tau); 
        sigma = sigma(1:r) - tau; Sigma = diag(sigma);
        U = U(:,1:r); V = V(:,1:r); 
        psvdtime = psvdtime + toc(tpsvd); 
    end

    % Step 10: SVT-MC Algorithm
    % This is an expensive routine for large matrices. A better option is  
    % to only compute the entries of X  that are needed. The purpose here  
    % is to compare different PSVD methods, not the overall efficiency.
    % ------------------------------------------------------------------------
    X = (U*Sigma)*V'; x = X(Omega); 
    
    % Step 11: SVT-MC Algorithm.
    % --------------------------
    relRes = norm(x-data)/norm_data;
    if (relRes < tol), break; end
     
    % Step 12: SVT-MC Algorithm.
    % --------------------------
    y = y + delta*(data-x);
    Y = sparse(indx_i,indx_j,y(indx),n1,n2);
     
end % for k loop Step 13: SVT-MC Algorithm

% Step 14: SVT-MC Algorithm  
% X_opt = X = U*Sigma*V' (already computed in step 10).
% -----------------------------------------------------

% Determining which method for output.
% ------------------------------------
fprintf(' \n');
fprintf('******************************************\n');
if method == 1, fprintf('PSVD method: %s\n','svt_svds'); end
if method == 2, fprintf('PSVD method: %s\n','svt_irlba'); end
if method == 3, fprintf('PSVD method: %s\n','svt'); end
if method == 4, fprintf('PSVD method: %s\n','svds'); end 
    
% Output results.
% ---------------
fprintf('total time = %0.5g\n',toc(tStart));  
fprintf('time for PSVD method: %0.5g\n',psvdtime);
fprintf('number of iterations: %d\n',k);
fprintf('number of samplings (m): %d\n',m);
fprintf('percentage of entries that are observed: %0.5g\n',perob);
fprintf('dimension of random M1: %d x %d \n',n1,r0);
fprintf('dimension of random M2: %d x %d \n',r0,n2);
fprintf('rank of M = M1*M2: %d \n',r0);
fprintf('The recovered rank: %d\n',length(diag(Sigma)) );
fprintf('The relative recovery error Frobenius norm: %d\n', norm(M-X,'fro')/norm(M,'fro'))
fprintf('******************************************\n');
fprintf(' \n');

end % example42

% Function to test if Octave is being used.
% MATLAB Central File Exchange.  ioxv.4623 (2024). Is this MATLAB or Octave? 
% https://www.mathworks.com/matlabcentral/fileexchange/23868-is-this-matlab-or-octave
% Retrieved March 14, 2024. 
% ------------------------------------------------------------------------------------
function uiIsOctave = uioctave
  uiIsOctave = false;  
  LIC = license('inuse');
  for elem = 1:numel(LIC)
    envStr = LIC(elem).feature;
    if strcmpi(envStr,'octave')
        uiIsOctave = true;
        break
    end
  end
end

% Matrix-vector product for svds avoid using full svd.
function z = Yfun(x,t,Y)
if strcmp(t,'notransp') % z = Y*x
   z = Y*x;
else % z = Y'*x
      z = Y'*x; 
end
end