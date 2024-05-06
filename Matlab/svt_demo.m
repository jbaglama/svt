function svt_demo
%  
%  DEMO ROUTINE FOR THE PAPER:
%  J. Baglama, J.Chávez-Casillas and V. Perovic, "A Hybrid Algorithm for 
%  Computing a Partial Singular Value Decomposition Satisfying a Given 
%  Threshold", submitted for publication 2024.
% 
%  TO RUN THE DEMO IN MATLAB OR OCTAVE:
%  >> svt_demo
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
% ----------------------------------------------------------------------

% Printing to screen the demo information. 
% ----------------------------------------
fprintf([
'*********************************************************************\n' ...
' DEMO:                                                               \n' ...
'  [U,S,V,FLAG] = svt_svds(A,PARAMETERS)   Matlab ONLY                \n' ...
'  [U,S,V,FLAG] = svt_irlba(A,PARAMETERS)  Matlab or Octave           \n' ...
'                                                                     \n' ...
'  A hybrid function for computing all singular triplets above a given\n' ...
'  threshold. The function repeatedly calls a Partial Singular Value  \n' ...
'  Decomposition (PSVD) method (svds or irlba) and a block SVD power  \n' ...
'  method to compute a PSVD of a matrix A with all the singular values\n' ...
'  above a threshold (sigma).                                         \n' ...
'                                                                     \n' ...
'  Demo runs examples with two matrices: illc1033 and bibd_20_10      \n' ...
'  Required matrix download from https://sparse.tamu.edu/             \n' ...
'  Demo will try to download matrices if not already in path          \n' ...
'                                                                     \n' ...
' PAPER:                                                              \n' ...
'  J. Baglama, J.Chavez-Casillas and V. Perovic, "A Hybrid Algorithm  \n' ...
'  for Computing a Partial Singular Value Decomposition Satisfying a  \n' ...
'  Given Threshold", submitted for publication 2024.                  \n' ...
'*********************************************************************\n'])
fprintf(' \n');

% Checking if svt_svds.m and svt_irlba.m exists. 
% ----------------------------------------------
if ~exist('svt_svds.m','file') || ~exist('svt_irlba.m','file') 
    error('svt_svds.m and svt_irlba.m are required'); 
end

% This demo function will run either svt_svds.m (MATLAB only) and
% svt_irlba.m (MATLAB and OCTAVE). It will also try to determine 
% if Octave is being used.
% ----------------------------------------------------------------
uiIsOctave = uioctave; 
if uiIsOctave
    usesvtsvds = 0;
    fprintf('Octave environment - using only svt_irlba.m\n');
else
    usesvtsvds = 1;
    fprintf('Matlab environment - using both svt_svds.m and svt_irlba.m\n');
end

% The following two matrices from the SuiteSparse Matrix Collection 
% https://sparse.tamu.edu/ will be used for this demo. Matrices from 
% SuiteSparse Matrix must be loaded in the path or the function will try 
% to download them. 
% illc1033  (1033 x 320) 
%   - largest singval ~ 2.144 and smallest singval ~ 0.00011353
%   - Large # of repeated singvals ~ 1 starting around singval(106)
% bibd_20_10  (190 x 184,756) -> A = A + LR'; L & R are randn(*,10) matrices
%   - largest singval ~ 7,133  and smallest singval ~ 113.37 
%   - Large # of repeated singvals ~ 113.446 starting around singval(31)\
% --------------------------------------------------------------------------

% Trying to get the matrices for the demo. 
% ----------------------------------------
fail(1:2) = 0;
  
% illc1033  (1033 x 320).
% ----------------------- 
matrix1 = 'illc1033';     
address1 = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/illc1033.mat';
A1 = getmat(matrix1,address1);
if isempty(A1), fail(1)=1; end
[m1,n1] = size(A1);

% bibd_20_10  (190 x 184,756) -> A = A + LR''; L & R are randn(*,10) matrices.
% ----------------------------------------------------------------------------
matrix2  = 'bibd_20_10';  
address2 = 'https://suitesparse-collection-website.herokuapp.com/mat/JGD_BIBD/bibd_20_10.mat';   
A2 = getmat(matrix2,address2);
matrix2 = 'bibd_20_10 + LR''; L & R are randn(*,10) matrices';
if isempty(A2), fail(2)=1; end
[m2,n2] = size(A2);
if ~fail(2)
    L = randn(m2,10); R = randn(n2,10); LR = L*R'; A2 = A2 + LR;
end

% Examples with matrix1: illc1033      
% -------------------------------
if ~fail(1)

   % Initializing the parameters.
   % ----------------------------
   k = 6;                    % Initial number of singular triplets.
   sigma = 1.5;              % Threshold. 
   tol = 1d-8;               % Tolerance used for convergence in irlba or svds.
   incre = 5;                % Initial increment.
   pwrsvd = 0;               % # iters of a block SVD power method.
   psvdmax = min(n1,m1);     % Maximum "size" of the output PSVD [U,S,V].
   display = 0;              % Displays some iters diagnostic.
   kmax = min(m1,n1);        % Maximum value k can reach.
   p0 = randn(max(n1,m1),1); % Starting vector for irlba and svds.

   % Calling the programs svt_irlba and svt_svds.
   % --------------------------------------------
   tStart = tic;
   [U,S,V,FLAG] = svt_irlba(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',...
                                 display,'p0',p0);
   irlba_time = toc(tStart);
   U1=[];S1=[];V1=[];FLAG1=[]; svds_time=[];
   if usesvtsvds
      tStart = tic; 
      [U1,S1,V1,FLAG1] = svt_svds(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',...
                                 display,'p0',p0);
      svds_time = toc(tStart);
   end
   fprintf(' \n');
   fprintf('Example 1: Compute all singvals above sigma = %0.5g \n',sigma);
   fprintf('*******************************************************************\n');
   display_ouput(A1,m1,n1,U,S,V,FLAG,U1,S1,V1,FLAG1,irlba_time,svds_time,...
                 k,sigma,tol,incre,pwrsvd,psvdmax,kmax,usesvtsvds,matrix1);
   fprintf('*******************************************************************\n');

   % Continuing to compute more singular values.
   % -------------------------------------------
   sigma = 0.9; 

   % Calling the programs svt_irlba and svt_svds again.
   % --------------------------------------------------
   tStart = tic;
   [U,S,V,FLAG] = svt_irlba(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display, ...
                                'U0',U,'V0',V,'S0',S,'p0',p0);
   irlba_time = toc(tStart);
   U1=[];S1=[];V1=[];FLAG1=[]; svds_time=[];
   if usesvtsvds
      tStart = tic; 
      [U1,S1,V1,FLAG1] = svt_svds(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display, ...
                                'U0',U1,'V0',V1,'S0',S1,'p0',p0);
      svds_time = toc(tStart);
   end
   fprintf(' \n');
   fprintf('Example 2: Continue Example 1: Compute all singvals above sigma = %0.5g \n',sigma);
   fprintf('*******************************************************************\n');
   display_ouput(A1,m1,n1,U,S,V,FLAG,U1,S1,V1,FLAG1,irlba_time,svds_time,...
                 k,sigma,tol,incre,pwrsvd,psvdmax,kmax,usesvtsvds,matrix1);
   fprintf('*******************************************************************\n');

   % Changing tolerance and sigma.
   % -----------------------------
   tol = 1d-4;
   sigma = 1.2;

   % Calling the programs svt_irlba and svt_svds.
   % --------------------------------------------
   tStart = tic;
   [U,S,V,FLAG] = svt_irlba(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',...
                                 display,'p0',p0);
   irlba_time = toc(tStart);
   U1=[];S1=[];V1=[];FLAG1=[]; svds_time=[];
   if usesvtsvds
      tStart = tic; 
      [U1,S1,V1,FLAG1] = svt_svds(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',...
                                 display,'p0',p0);
      svds_time = toc(tStart);
   end
   fprintf(' \n');
   fprintf('Example 3: Compute all singvals above sigma = %0.5g \n',sigma);
   fprintf('*******************************************************************\n');
   display_ouput(A1,m1,n1,U,S,V,FLAG,U1,S1,V1,FLAG1,irlba_time,svds_time,...
                 k,sigma,tol,incre,pwrsvd,psvdmax,kmax,usesvtsvds,matrix1);
   fprintf('*******************************************************************\n');

   % Continuing to compute more singular values.
   % -------------------------------------------
   sigma = 0.9; 
   pwrsvd = 1;

   % Calling the programs svt_irlba and svt_svds again.
   % --------------------------------------------------
   tStart = tic;
   [U,S,V,FLAG] = svt_irlba(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display, ...
                                'U0',U,'V0',V,'S0',S,'p0',p0);
   irlba_time = toc(tStart);
   U1=[];S1=[];V1=[];FLAG1=[]; svds_time=[];
   if usesvtsvds
      tStart = tic; 
      [U1,S1,V1,FLAG1] = svt_svds(A1,'m',m1,'n',n1,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display, ...
                                'U0',U1,'V0',V1,'S0',S1,'p0',p0);
      svds_time = toc(tStart);
   end
   fprintf(' \n');
   fprintf('Example 4: Continue Example 3 with power svd iteration to compute all singvals above sigma = %0.5g \n',sigma);
   fprintf('*******************************************************************\n');
   display_ouput(A1,m1,n1,U,S,V,FLAG,U1,S1,V1,FLAG1,irlba_time,svds_time,...
                 k,sigma,tol,incre,pwrsvd,psvdmax,kmax,usesvtsvds,matrix1);
   fprintf('*******************************************************************\n');

end % end matrix1: illc1033

% Examples with matrix2: bibd_20_10      
% --------------------------------
if ~fail(2)

   % Initializing the parameters.
   % ----------------------------
   k = 6;                    % Initial number of singular triplets.
   sigma = 4000;             % Threshold.
   tol = 1d-8;               % Tolerance used for convergence in irlba or svds.
   incre = 5;                % Intial increment.
   pwrsvd = 0;               % # iters of a block SVD power method.
   psvdmax = min(n2,m2);     % Maximum "size" of the output PSVD [U,S,V].
   display = 0;              % Displays some iters diagnostic.
   kmax = min(m2,n2);        % Maxiumum value k can reach.
   p0 = randn(max(n2,m2),1); % Starting vector for irlba and svds.

   % Calling the programs svt_irlba and svt_svds.
   % --------------------------------------------
   tStart = tic;
   [U,S,V,FLAG] = svt_irlba(A2,'m',m2,'n',n2,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',...
                                 display,'p0',p0);
   irlba_time = toc(tStart);
   U1=[];S1=[];V1=[];FLAG1=[]; svds_time=[];
   if usesvtsvds
      tStart = tic; 
      [U1,S1,V1,FLAG1] = svt_svds(A2,'m',m2,'n',n2,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',...
                                 display,'p0',p0);
      svds_time = toc(tStart);
   end
   fprintf(' \n');
   fprintf('Example 5: Compute all singvals above sigma = %0.5g \n',sigma);
   fprintf('*******************************************************************\n');
   display_ouput(A2,m2,n2,U,S,V,FLAG,U1,S1,V1,FLAG1,irlba_time,svds_time,...
                 k,sigma,tol,incre,pwrsvd,psvdmax,kmax,usesvtsvds,matrix2);
   fprintf('*******************************************************************\n');

   % Continuing to compute more singular values.
   % --------------------------------------------
   sigma = 200;           

   % Calling the programs svt_irlba and svt_svds again.
   % --------------------------------------------------
   tStart = tic;
   [U,S,V,FLAG] = svt_irlba(A2,'m',m2,'n',n2,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display, ...
                                'U0',U,'V0',V,'S0',S,'p0',p0);
   irlba_time = toc(tStart);
   U1=[];S1=[];V1=[];FLAG1=[]; svds_time=[];
   if usesvtsvds
      tStart = tic; 
      [U1,S1,V1,FLAG1] = svt_svds(A2,'m',m2,'n',n2,'sigma',sigma,'tol',tol,...
                                'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                                'incre',incre,'pwrsvd',pwrsvd,'display',display, ...
                                'U0',U1,'V0',V1,'S0',S1,'p0',p0);
      svds_time = toc(tStart);
   end
   fprintf(' \n');
   fprintf('Example 6: Continue Example 5: Compute all singvals above sigma = %0.5g \n',sigma);
   fprintf('*******************************************************************\n');
   display_ouput(A2,m2,n2,U,S,V,FLAG,U1,S1,V1,FLAG1,irlba_time,svds_time,...
                 k,sigma,tol,incre,pwrsvd,psvdmax,kmax,usesvtsvds,matrix2);
   fprintf('*******************************************************************\n');

end % end matrix2: bibd_20_10 

end % svt_demo
  
% Display the output values.
%---------------------------
function display_ouput(A,m,n,U,S,V,FLAG,U1,S1,V1,FLAG1,irlba_time,svds_time,...
                       k,sigma,tol,incre,pwrsvd,psvdmax,kmax,usesvtsvds,name)

   % Computing the errors for output.
   % --------------------------------
   outirlba(1) = -1; outirlba(2) = -1;
   if ~isempty(S) % Error for svt_irlba
      outirlba(1) = sqrt(norm(A*V - U*S)^2+norm(A'*U - V*S)^2); 
      outirlba(2) = sqrt(norm(V'*V-eye(size(V,2)))^2+norm(U'*U-eye(size(U,2)))^2);
   end
   outsvds(1) = -1; outsvds(2) = -1;
   if ~isempty(S1) % Error for svt_svds
      outsvds(1) = sqrt(norm(A*V1 - U1*S1)^2+norm(A'*U1 - V1*S1)^2); 
      outsvds(2) = sqrt(norm(V1'*V1-eye(size(V1,2)))^2+norm(U1'*U1-eye(size(U1,2)))^2);
   end
   fprintf('Matrix A (%d x %d): %s \n',m,n,name);
   fprintf(' Parameters: \n');
   fprintf('   sigma = %0.5g \n',sigma);
   fprintf('   tol = %0.5g \n',tol);
   fprintf('   psvdmax = %d \n',psvdmax);
   fprintf('   kmax = %d \n',kmax);
   fprintf('   pwrsvd = %d \n',pwrsvd);
   fprintf('   incre = %d \n',incre);
   fprintf('   k = %d \n',k);
   fprintf('-------------------\n');
   fprintf('SVT_IRLBA \n'); 
   fprintf('   FLAG = %d\n',FLAG);
   fprintf('   max. singval =  %0.5g\n',max(diag(S)));
   fprintf('   min. singval (>=%0.5g)  =  %0.5g\n',sigma,min(diag(S)));
   fprintf('   # singvals = %0.5g\n',size(S,1));
   fprintf('   cputime    = %0.5g secs\n',irlba_time);
   if outirlba(1) >= 0
      fprintf('   sqrt(||A*V - U*S||^2+||A^T*U - V*S||^2) = %0.5g\n',outirlba(1));
      fprintf('   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n',outirlba(2));
   else
       fprintf('  svt_irlba output values were empty. \n');
   end
   if usesvtsvds
      fprintf('SVT_SVDS \n'); 
      fprintf('   FLAG = %d\n',FLAG1);
      fprintf('   max. singval =  %0.5g\n',max(diag(S1)));
      fprintf('   min. singval (>=%0.5g)  =  %0.5g\n',sigma,min(diag(S1)));
      fprintf('   # singvals = %0.5g\n',size(S1,1));
      fprintf('   cputime    = %0.5g secs\n',svds_time);
      if outsvds(1) >= 0
         fprintf('   sqrt(||A*V - U*S||^2+||A^T*U - V*S||^2) = %0.5g\n',outsvds(1));
         fprintf('   sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n',outsvds(2));
      else
         fprintf('  svt_svds output values were empty. \n');
      end
   end
end

% Function to try to download matrix from SuiteSparse Matrix Collection.
% ----------------------------------------------------------------------
function A = getmat(name,address)
  name_mat = strcat(name,'.mat'); 
   if ~exist(name_mat,'file')
      try
         fprintf(strcat(name_mat,'\n')); 
         fprintf('matrix not in path \n');
         fprintf('Try getting matrix from web ... \n');
         % urlwrite is required for Octave. In MATLAB (R2024a) it is 
         % still usable (but not recommended). MATLAB recommends using
         % the function websave - however websave is not a function in Octave.
         % websave(name_mat,address);
         % -------------------------------------------------------------------
         urlwrite(address,name_mat); 
         load(name_mat); A = Problem.A; clear Problem;
      catch
         fprintf(strcat(name_mat,'\n'));
         fprintf('FAILED getting matrix from web and does not exist in path \n');
         fprintf('Recommend directly downloading the file from the web and placing the file in the path \n');
         A = []; 
      end       
   else
     load(name_mat); A = Problem.A; clear Problem;
   end
end

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