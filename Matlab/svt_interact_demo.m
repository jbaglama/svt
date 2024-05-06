function svt_interact_demo
%  
%  INTERACTIVE DEMO FUNCTION FOR THE PAPER:
%  J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm for 
%  Computing a Partial Singular Value Decomposition Satisfying a Given 
%  Threshold", submitted for publication 2024.
% 
%  TO RUN DEMO IN MATLAB OR OCTAVE:
%  >> svt_interact_demo
%
%  REQUIRED SOFTWARE:  
%  svt_irlba.m  and svt_svds.m  
%  https://github.com/jbaglama/svt/
%
%  DATE LAST MODIFIED: 
%  5/1/24
%
%  LANGUAGE:
%  MATLAB versions: R2018b ... R2024a (earlier releases may also work)
%  OCTAVE versions: 8.4.0 (earlier releases may also work)
%  
%  AUTHORS: 
%  James Baglama            email: jbaglama@uri.edu
%  Jonathan Chávez-Casillas email: jchavezc@uri.edu
%  Vasilije Perovic         email: perovic@uri.edu
% ---------------------------------------------------

fprintf([
'*********************************************************************\n' ...
' INTERACTIVE DEMO:                                                   \n' ...
'  [U,S,V,FLAG] = svt_svds(A,PARAMETERS)   Matlab ONLY                \n' ...
'  [U,S,V,FLAG] = svt_irlba(A,PARAMETERS)  Matlab or Octave           \n' ...
'                                                                     \n' ...
'  A hybrid function for computing all singular triplets above a given\n' ...
'  threshold. The function repeatedly calls a Partial Singular Value  \n' ...
'  Decomposition (PSVD) method (svds or irlba) and a block SVD power  \n' ...
'  method to compute a PSVD of a matrix A with all the singular values\n' ...
'  above a threshold (sigma).                                         \n' ...
'                                                                     \n' ...
'  This interactive demo will allow the user to explore svt_svds      \n' ...
'  or svt_irlba with varying thresholds and options for 7 different   \n' ...
'  matrices from the SuiteSparse Matrix Collection.                   \n' ... 
'  This demo requires to download the matrices from                   \n' ...
'  https://sparse.tamu.edu/                                           \n' ...
'  This demo will try to download matrices if not already in path.    \n' ...
'                                                                     \n' ...
' PAPER:                                                              \n' ...
'  J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm  \n' ...
'  for Computing a Partial Singular Value Decomposition Satisfying a  \n' ...
'  Given Threshold", submitted for publication 2024.                  \n' ...
'*********************************************************************\n'])
fprintf(' \n');

% Checking if svt_svds.m and svt_irlba.m exists.
% ----------------------------------------------
if ~exist('svt_svds.m','file') || ~exist('svt_irlba.m','file')
    error('svt_svds.m and svt_irlba.m are required');
end

% Selecting function svt_svds.m or svt_irlba.m to run driver program.
% -------------------------------------------------------------------
fprintf('1. svt_irlba.m  Matlab or Octave  (default - enter) \n');
fprintf('2. svt_svds.m   Matlab ONLY \n');
prompt = 'Select svt_svds.m or svt_irlba.m to run with driver program: ';
svt = input(prompt);
if isempty(svt) || svt ~= 2, svt = 1; end
uiIsOctave = uioctave;
if uiIsOctave && svt == 2
    fprintf('Octave environment - cannot use svt_svds.m\n');
    svt = 1;
end

% This interactive demo will allow the user to explore svt_svds.m or
% svt_irlba.m with varying thresholds and options for 7 different matrices
% taken from the SuiteSparse Matrix Collection at: https://sparse.tamu.edu/ 
% These matrices from SuiteSparse Matrix must be loaded in the MATLAB path 
% or the function will try to download them. Internet connection is required
% for the download. If no internet connection is available, the only example
% that can be run is the diagonal example 1.
%
% 1. diag(1:500)    500 x 500
% 2. illc1033      1033 x 320
% 3. well1850     1,850 x 712
% 4. bibd_20_10     190 x 184,756 -> A = A + LR'; L & R are randn(*,10) matrices
% 5. mhd4800b     4,800 x 4,800   -> A = A + LR'; L & R are randn(*,10) matrices
% 6. Maragal_2      555 x 350
% -----------------------------------------------------------------
fprintf(' \n');
fprintf('1. diag(1:500)    500 x 500 (default - enter) \n');
fprintf('2. illc1033      1033 x 320 \n');
fprintf('3. well1850     1,850 x 712 \n');
fprintf('4. bibd_20_10     190 x 184,756 -> A = A + LR''; L & R are randn(*,10) matrices \n');
fprintf('5. mhd4800b     4,800 x 4,800   -> A = A + LR''; L & R are randn(*,10) matrices \n');
fprintf('6. Maragal_2      555 x 350 \n');
fprintf(' \n');
prompt = 'Select a matrix (1 - 6): ';
example = input(prompt);
if isempty(example), example = 1; end

% Getting the matrix to use for the examples.
% -------------------------------------------

% Setting the value to determine if a matrix cannot be used.
% If it fails to load a matrix, the program defaults to diag(1:500).
% ------------------------------------------------------------------
fail = 0;

% 1. diag(1:500)   500 x 500
% --------------------------
if example == 1
   A=sparse(diag(1:500)); name = 'diag(1:500)';
   leftendpt = 0; rightendpt = 500;
end

% 2. illc1033  1033 x 320
% -----------------------
if example == 2
   name = 'illc1033';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/illc1033.mat';
   leftendpt = 0; rightendpt = 2.1;
end

% 3. well1850  1,850 x 712
% ------------------------
if example == 3
   name = 'well1850';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/well1850.mat';
   leftendpt = 0; rightendpt = 1.7;
end

% 4. bibd_20_10  190 x 184,756 -> A = A + LR''; L & R are randn(*,10) matrices
% ---------------------------------------------------------------------------
if example == 4
   name = 'bibd_20_10';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/JGD_BIBD/bibd_20_10.mat';
   leftendpt = 113.4; rightendpt = 7000;
end

% 5. mhd4800b  4,800 x 4,800   -> A = A + LR''; L & R are randn(*,10) matrices
% -----------------------------------------------------------------------------
if example == 5
   name = 'mhd4800b';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Bai/mhd4800b.mat';
   leftendpt = 0.1; rightendpt = 5000;
end

% 7. Maragal_2  555 x 350
% ------------------------
if example == 6
   name = 'Maragal_2';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/NYPA/Maragal_2.mat';
   leftendpt = 0; rightendpt = 10;
end

% Checking to see if the matrices are in the path - if not, it will try to 
% download the matrix from the web. If it is not in path and the download
% fails, it will switch to use A=sparse(diag(1:500))
% -------------------------------------------------------------------------
if example ~= 1
   name_mat = strcat(name,'.mat');
   if ~exist(name_mat,'file')
      try
         fprintf(strcat(name_mat,'\n'));
         fprintf('matrix not in path \n');
         fprintf('Try getting matrix from web ... \n');
         % urlwrite is required for Octave and in MATLAB (R2024a) it is
         % still usable (but not recommended). MATLAB recommends using
         % the function websave. However websave is not a function in Octave.
         % websave(name_mat,address);
         % ------------------------------------------------------------------
         urlwrite(address,name_mat);
      catch
         fprintf('FAILED getting matrix from web and does not exist in path \n');
         fprintf('Recommend directly downloading the file from the web and placing the file in the path \n');
         fprintf('Changing to diagonal matrix diag(1:500) \n');
         A=sparse(diag(1:500)); name = 'diag(1:500)'; fail = 1;
         leftendpt = 0; rightendpt = 500;
      end
   end
end

% Loading matrix from cell and clear from memory to save space.
% -------------------------------------------------------------
if example ~= 1 && ~fail, load(name_mat); A = Problem.A; clear Problem; end

% Checking again matrix exists in path.
% -------------------------------------
if isempty(A)
   fprintf('FAILED getting matrix from web and does not exist in path \n');
   fprintf('Changing to diagonal matrix diag(1:500) \n');
   A = sparse(diag(1:500)); name = 'diag(1:500)';
   leftendpt = 0; rightendpt = 500;
   example = 1;
end

% Getting size of the matrix A
% ----------------------------
[m,n] = size(A);

% Adding random matrices to matrix A for examples 4 and 5.
% --------------------------------------------------------
if example == 4 || example == 5
   L = randn(m,10); R = randn(n,10); LR = L*R'; A = A + LR;
end

% Getting the threshold value.
% ----------------------------
fprintf(' \n');
prompt = ['Select a threshold value for matrix ',name,' between ', num2str(leftendpt),' and ',num2str(rightendpt),': '];
thres = input(prompt);
if thres > rightendpt || thres < leftendpt
   error_str = ['ERROR: Threshold value not between ', num2str(leftendpt),' and ',num2str(rightendpt),' '];
   error(error_str);
end

% Getting the parameter values.
% -----------------------------
fprintf(' \n');
prompt = 'Select some parameters? Y/N (default N - enter):';
str = input(prompt,'s');
if isempty(str) || ~ischar(str), str = 'N'; end
str = upper(str); if ~strcmp(str,'Y'), str = 'N'; end
if strcmp(str,'Y')
   % Tolerance used for convergence in irlba or svds.
   prompt = 'TOL: Select tolerance (eps <= tol <= 1d-4) (default: sqrt(eps) - enter): ';
   tol = input(prompt); if isempty(tol), tol = sqrt(eps); end
   if tol < eps || tol > 1d-4, error('ERROR: Incorrect TOL value'); end
   % Initial number of singular triplets.
   prompt = 'K: (0 < k <= 20) of initial singular triplets (default: 6 - enter): ';
   k = input(prompt); if isempty(k), k = 6; end
   if k <= 0 || k > 20, error('ERROR: Incorrect K value'); end
   % Intial increment.
   prompt = 'INCRE: (0<incre<=10)intial increment (default: 5 - enter): ';
   incre = input(prompt); if isempty(incre), incre = 5; end
   if incre <= 0 || incre > 10, error('ERROR: Incorrect INCRE value'); end
   % # iters of a block SVD power method.
   prompt = 'pwrsvd: number (0<=pwrsvd <=10) of BLK SVD power iters. (default: 0 - enter): ';
   pwrsvd = input(prompt); if isempty(pwrsvd), pwrsvd = 0; end
   if pwrsvd < 0 || pwrsvd > 10, error('ERROR: Incorrect pwrsvd value'); end
else
 % Setting to default values if there is no user input.
 tol = 1d-8; k = 6; incre = 5; pwrsvd = 0;
end
% Setting the following parameters (not interactive choices) for the demo.
psvdmax = min(n,m);   % Maximum "size" of the output PSVD [U,S,V].
display = 0;          % Displays some iters diagnostic.
kmax = min(m,n);      % Maxiumum value k can reach.

% Initializing input values U0,S0,V0.
% -----------------------------------
U=[]; V=[]; S=[]; compute_more = 0; tot_cputime = 0;

% Call svt_svds.m or svt_irlba.m
% ------------------------------
while 1

      tStart = tic;
      if svt == 1 % Calling svt_irlba.m
         if isempty(S) 
             [U,S,V,FLAG] = svt_irlba(A,'m',m,'n',n,'sigma',thres,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display);
         else
             [U,S,V,FLAG] = svt_irlba(A,'m',m,'n',n,'sigma',thres,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display,...
                                 'U0',U,'V0',V,'S0',S);
         end
      else  % Calling svt_svds.m
         if isempty(S)
            [U,S,V,FLAG] = svt_svds(A,'m',m,'n',n,'sigma',thres,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display);
         else
             [U,S,V,FLAG] = svt_svds(A,'m',m,'n',n,'sigma',thres,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display,...
                                'U0',U,'V0',V,'S0',S);
         end
      end
      outsvt(1) = FLAG;
      outsvt(2) = toc(tStart);  tot_cputime = tot_cputime + outsvt(2);
      if ~isempty(S)
         outsvt(3) = sqrt(norm(A*V - U*S)^2+norm(A'*U - V*S)^2);
         outsvt(4) = sqrt(norm(V'*V-eye(size(V,2)))^2+norm(U'*U-eye(size(U,2)))^2);
      else
        outsvt(3) = -1; outsvt(4)= -1;
      end
      fprintf(' \n');
      fprintf('********************************************************************** \n');
      fprintf('Matrix: %s \n',name);
      if svt == 1, fprintf('SVT_IRLBA \n'); else, fprintf('SVT_SVDS \n'); end
      fprintf('  FLAG = %d\n',outsvt(1));
      fprintf('  threshold = %0.5g\n',thres);
      fprintf('  max. singval =  %0.5g\n',max(diag(S)));
      fprintf('  min. singval (>=%0.5g)  =  %0.5g\n',thres,min(diag(S)));
      fprintf('  # singvals = %0.5g\n',size(S,1));
      fprintf('  cputime for iteration = %0.5g secs\n',outsvt(2));
      if compute_more
         fprintf('  total cputime = %0.5g secs \n',tot_cputime);
      end
      fprintf('  sqrt(||A*V - U*S||^2+||A^T*U - V*S||^2) = %0.5g\n',outsvt(3));
      fprintf('  sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n',outsvt(4));
      fprintf('********************************************************************** \n');
      fprintf(' \n');
      if k >= min(n,m), fprintf('Max. number of singvals computed \n'); return; end

      % Determining if it is needed to get more singular values
      % -------------------------------------------------------
      prompt = 'Continue - change threshold and compute values? Y/N (default N - enter): ';
      str = input(prompt,'s');
      if isempty(str) || ~ischar(str), str = 'N'; end
      str = upper(str); if ~strcmp(str,'Y'), str = 'N'; end
      if strcmp(str,'N'), return; end; compute_more = 1;

      % Getting a new threshold value
      % -----------------------------
      fprintf(' \n');
      rightendpt = thres;
      prompt = ['Select a threshold value for matrix ',name,' between ', num2str(leftendpt),' and ',num2str(rightendpt),': '];
      thres = input(prompt);
      if thres > rightendpt || thres < leftendpt
        error_str = ['ERROR: Threshold value not between ', num2str(leftendpt),' and ',num2str(rightendpt),' '];
        error(error_str);
      end
      if abs(rightendpt - leftendpt) < sqrt(eps)
         fprintf('[leftendpt rightendpt] < sqrt(eps) - exit. \n'); return;
      end
end
end % svt_interact_demo function

% Function to test if Octave is being used.
% MATLAB Central File Exchange.  ioxv.4623 (2024). Is this MATLAB or Octave?
% https://www.mathworks.com/matlabcentral/fileexchange/23868-is-this-matlab-or-octave
% Retrieved March 14, 2024.
% ------------------------------------------------------------------------
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