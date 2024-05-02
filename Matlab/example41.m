function example41
%
%  PAPER:
%  J. Baglama, J.Chavez-Casillas and V. Perovic, "A Hybrid Algorithm for
%  Computing a Partial Singular Value Decomposition Satisfying a Given
%  Threshold", submitted for publication 2024.
%
%  EXAMPLE 4.1 FROM PAPER:
%  Using the codes svt_svds.m, svt_irlba.m, and svt.m, this example demonstrates
%  the algorithm of the paper using 12 matrices with varying sigma values and 
%  the results are reported within the paper. More information about each 
%  matrix and choice of sigma value is given below where the matrix is loaded
%  into the MATLAB workspace. Also, a comparison with the code svt.m from the
%  paper by C. Li and H. Zhou titled "SVT: Singular value thresholding in MATLAB",
%  J. Statist. Softw., 81 (2017), pp. 1–12. All matrices are downloaded
%  from the SuiteSparse Matrix Collection: http://sparse.tamu.edu/
%
%  REQUIRED SOFTWARE:
%  1. svt_irlba.m  and svt_svds.m
%     https://github.com/jbaglama/svt/
%  2. svt.m
%     https://github.com/Hua-Zhou/svt/%
%
%  LANGUAGE:
%  MATLAB versions: R2018b ... R2023b (earlier releases will not work)
%  OCTAVE versions: 8.4.0 (earlier releases may also work)
%
%
%  DATE LAST MODIFIED:
%  5/1/24
%
%  AUTHORS:
%  James Baglama            email: jbaglama@uri.edu
%  Jonathan Chávez-Casillas email: jchavezc@uri.edu
%  Vasilije Perovic         email: perovic@uri.edu
% ---------------------------------------------------

%system('caffeinate -dims &');

% Test if MATLAB or Octave is being used.
% ---------------------------------------
uiIsOctave = uioctave;

% Check if svt.m., svt_svds.m and svt_irlba.m exists
% --------------------------------------------------
if ~uiIsOctave
   if ~exist('svt_svds.m','file'),  error('svt_svds.m requried for Example41'); end
end
if ~exist('svt_irlba.m','file'), error('svt_irlba.m requried for Example41'); end
if ~exist('svt.m','file'),       error('svt.m requried for Example41'); end


% The list of the 12 matrices used is detailed below. Below, (SVT paper) is
% next to each matrix taken from the Li and Zhou paper. All matrices were 
% downloaded from SuiteSparse Matrix Collection: http://sparse.tamu.edu/
% ----------------------------------------------------------------
% Ordered by # of rows
% 1.  bibd_20_10    size: 190 x 184,756       (SVT paper)
% 2.  bibd_22_8     size: 231 x 319,770       (SVT paper)
% 3.  bfwb398       size: 398 x 398           (SVT paper)
% 4.  mhd4800b      size: 4,800 x 4,800       (SVT paper)
% 5.  cryg10000     size: 10,000 x 10,000     (SVT paper)
% 6.  stormG2_1000  size: 528,185 x 1,377,306 (SVT paper)
% 7.  Maragal_2     size: 555 x 350
% 8.  illc1033      size: 1,033 x 320
% 9.  well1850      size: 1,850 x 712
% 10. JP            size: 87,616 x 67,320
% 11. Rel8          size: 345,688 x 12,347
% 12. Rucci1        size: 1,977,885 x 109,900
% -------------------------------------------------------------------------

% Pick a matrix.
% -------------
matrix = 1;

% For Example 4.1 in the paper each routine will be called 'iter' number of
% times with a new random vector. svt_svds.m and svt_irlba.m will use the 
% same random vector on each restart, however svt.m does not have an option
% to specify the starting vector. The output will be the mean overall CPU 
% time and max errors for sqrt(||AV-US||^2+||A^TU-VS||^2) and
% sqrt(||V^TV-I||^2+||U^TU-I||^2).
% -------------------------------------------------------------------------
iter = 1; % For the paper, 'iter' is set to 10

% Set the parameters
% ------------------
k = 6;           % initial number of singular triplets (all routines)
sigma = [];      % threshold - unique to each matrix - see below for value (all routines)
tol = 1d-8;      % tolerance (all routines)
incre = 5;       % initial increment (all routines)
pwrsvd = 0;      % # iters of a block SVD power method
psvdmax = 800;   % max size of the output PSVD
display = 0;     % displays some iters diagnostic
kmax = 100;      % maximum value k can reach
p0 = [];         % starting vector for irlba and svds - reset each iter to randn
fprintf(' \n');

% 1. bibd_20_10  size: 190 x 184,756
%    From SVT paper: In this example, A = A + LR' where L and R are n x 10
%    normally distributed pseudo-random numbers. That is, they are obtained
%    from a call to the function 'randn'. In this example, we use the the
%    same random stream as in SVT paper.
%    Largest SV is on the order of 1,000 while the smallest one is on the 
%    order of 100 with many repeated SVs.
% -------------------------------------------------------------------------
if matrix == 1
   name = 'bibd_20_10';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/JGD_BIBD/bibd_20_10.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('bibd_20_10  size: 190 x 184,756\n');
   [m,n] = size(A);
   if ~uiIsOctave
      s = RandStream('mt19937ar','Seed',2014); RandStream.setGlobalStream(s);
      L = randn(m,10); R = randn(n,10);
   else
      L = randn(m,10); R = randn(n,10);
   end
   A = A + L*R';
   sigma = 466; % set sigma - captures the top 20 singular values.
end

% 2. bibd_22_8  size: 231 x 319,770
%    From SVT paper: In this example, A = A + LR' where L and R are n x 10
%    normally distributed pseudo-random numbers. That is, they are obtained
%    from a call to the function 'randn'. In this example, we use the the
%    same random stream as in SVT paper.
%    Largest SV is on the order of 10,000 while the smallest one is on the 
%    order of 100 with many repeated SVs.
% -------------------------------------------------------------------------
if matrix == 2
   name = 'bibd_22_8';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/JGD_BIBD/bibd_22_8.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('bibd_22_8  size: 231 x 319,770\n');
   [m,n] = size(A);
   if ~uiIsOctave
      s = RandStream('mt19937ar','Seed',2014); RandStream.setGlobalStream(s);
      L = randn(m,10); R = randn(n,10);
   else
      L = randn(m,10); R = randn(n,10);
   end
   A = A + L*R';
   sigma = 435.79; % set sigma - captures the top 20 singular values.
end

% 3. bfwb398 size: 398 x 398
%    From SVT paper: In this example, A = A + LR' where L and R are n x 10
%    normally distributed pseudo-random numbers. That is, they are obtained
%    from a call to the function 'randn'. In this example, we use the the
%    same random stream as in SVT paper.
%    Largest SV is on the order of 1d-5 while the smallest one is on the 
%    order of 1d-6.
% -------------------------------------------------------------------------
if matrix == 3
   name = 'bfwb398';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Bai/bfwb398.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('bfwb398 size: 398 x 398\n');
   [m,n] = size(A);
   if ~uiIsOctave
      s = RandStream('mt19937ar','Seed',2014); RandStream.setGlobalStream(s);
      L = randn(m,10); R = randn(n,10);
   else
      L = randn(m,10); R = randn(n,10);
   end
   A = A + L*R';
   sigma = 2.2d-5; % set sigma - captures the top 50 singular values.
end

% 4. mhd4800b size: 4,800 x 4,800
%    From SVT paper: In this example, A = A + LR' where L and R are n x 10
%    normally distributed pseudo-random numbers. That is, they are obtained
%    from a call to the function 'randn'. In this example, we use the the
%    same random stream as in SVT paper.
%    Largest SV is on the order of 2 while the 50th SV is on the order of 0.09
% ----------------------------------------------------------------------------
if matrix == 4
   name = 'mhd4800b';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Bai/mhd4800b.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('mhd4800b size: 4,800 x 4,800\n');
   [m,n] = size(A);
  if ~uiIsOctave
      s = RandStream('mt19937ar','Seed',2014); RandStream.setGlobalStream(s);
      L = randn(m,10); R = randn(n,10);
   else
      L = randn(m,10); R = randn(n,10);
   end
   A = A + L*R';
   sigma = 0.122; % set sigma - captures the top 50 singular values.
end

% 5. cryg10000 size: 10,000 x 10,000
%    From SVT paper: In this example, A = A + LR' where L and R are n x 10
%    normally distributed pseudo-random numbers. That is, they are obtained
%    from a call to the function 'randn'. In this example, we use the the
%    same random stream as in SVT paper.
%    Largest SV is on the order of 40,000 while the 50th SV is on the order 
%    of 20,000
% -------------------------------------------------------------------------
if matrix == 5
   name = 'cryg10000';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Bai/cryg10000.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('cryg10000 size: 10,000 x 10,000\n');
   [m,n] = size(A);
   if ~uiIsOctave
      s = RandStream('mt19937ar','Seed',2014); RandStream.setGlobalStream(s);
      L = randn(m,10); R = randn(n,10);
   else
      L = randn(m,10); R = randn(n,10);
   end
   A = A + L*R';
   sigma = 21757; % set sigma - captures the top 50 singular values.
end

% 6. stormG2_1000  size: 528,185 x 1,377,306
%    From SVT paper:  Largest SV = 3288 and the 100th SV = 6.7768
% ---------------------------------------------------------------
if matrix == 6
   name = 'stormG2_1000';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Mittelmann/stormG2_1000.mat';
   A = getmat(name,address);
   fprintf('stormG2_1000  size: 528,185 x 1,377,306\n');
   if isempty(A), return; end
   [m,n] = size(A);
   sigma = 632.4603; % set sigma - captures the top 50 singular values.
end

% 7. Maragal_2 size: 555 x 350
%    Largest SV = 10.29508. Its numerical rank is 220
% ---------------------------------------------------
if matrix == 7
   name = 'Maragal_2';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/NYPA/Maragal_2.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('Maragal_2 size: 555 x 350\n');
   [m,n] = size(A);
   sigma = 1d-10; % set sigma - captures the 171 singular values.
end

% 8. illc1033  size: 1,033 x 320
%    Largest SV = 2.144 while the smallest SV = 0.00011353. This matrix has
%    a large number of repeating SV around 1, starting from the 106th SV.
% -------------------------------------------------------------------------
if matrix == 8
   name = 'illc1033';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/illc1033.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('illc1033  size: 1,033 x 320\n');
   [m,n] = size(A);
   sigma = 0.9; % set sigma - captures the top 197 SVs large number repeated ~ 1.
end

% 9. well1850  size: 1,850 x 712
%    Largest SV = 1.793 while the smallest SV = 0.0161. This matrix has a 
%    large number of repeating SV around 1, starting from the 265th SV.
% -------------------------------------------------------------------------
if matrix == 9
   name = 'well1850';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/HB/well1850.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('well1850  size: 1,850 x 712\n');
   [m,n] = size(A);
   sigma = 0.0; % set sigma -captures all SVs large number repeated ~ 1.
end

% 10. JP size: 87,616 x 67,320
%    Largest SV = 4223.096, 50th SV = 1979.547 and 100th SV = 1583.588.
% ---------------------------------------------------------------------
if matrix == 10
   name = 'JP';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Harvard_Seismology/JP.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('JP size: 87,616 x 67,320\n');
   [m,n] = size(A);
   sigma = 1500; % set sigma - captures the top 126 SVs
end

% 11. Rel8 size: 345,688 x 12,347
% Largest SV = 18.2733, 50th SV = 12.4838, 100th SV = 12.1186, and
% 150th SV = 11.892
% ----------------------------------------------------------------
if matrix == 11
   name = 'rel8';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/JGD_Relat/rel8.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('Rel8 size: 345,688 x 12,347\n');
   [m,n] = size(A);
   sigma = 12.5; % set sigma - captures the top 47 singular values.
end

% 12. Rucci1 size: 1,977,885 x 109,900
%    Largest SV = 7.06874,  10th SV =  6.73798 and 25th SV = 6.5622
% -----------------------------------------------------------------
if matrix == 12
   name = 'Rucci1';
   address = 'https://suitesparse-collection-website.herokuapp.com/mat/Rucci/Rucci1.mat';
   A = getmat(name,address);
   if isempty(A), return; end
   fprintf('Rucci1 size: 1,977,885 x 109,900\n');
   [m,n] = size(A);
   sigma = 6.5; % set sigma - captures the top 33 SVs
end

% Common random starting vector for svt_svds and svt_irlba for each iteration.
% ---------------------------------------------------------------------------
p0 = randn(max(n,m),iter);

% Calling svt_svds iter times and output the results for paper.
% -------------------------------------------------------------
if ~uiIsOctave
   fprintf(' \n');
   fprintf('SVT_SVDS\n');
   for i = 1:iter
      tStart = tic;
      [U,S,V,FLAG] = svt_svds(A,'m',m,'n',n,'sigma',sigma,'tol',tol,...
                            'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                            'incre',incre,'pwrsvd',pwrsvd,...
                            'display',display,'p0',p0(:,i));

       if FLAG ~= 0
         fprintf('    iter = %d\n',iter);
         fprintf('    FLAG = %d\n',FLAG);
       end
       cputime(i) = toc(tStart);
       sizeS(i)   = size(S,1);
       if ~isempty(S)
         PSVD_err(i) = sqrt(norm(A*V - U*S)^2+norm(A'*U - V*S)^2);
         VU_err(i)   = sqrt(norm(V'*V-eye(size(V,2)))^2+norm(U'*U-eye(size(U,2)))^2);
       else
         PSVD_err(i) = Inf; VU_err(i) = Inf;
       end
    end
    fprintf('    threshold = %0.5g\n',sigma);
    fprintf('    max # singvals = %0.5g\n',max(sizeS));
    fprintf('    min # singvals = %0.5g\n',min(sizeS));
    fprintf('    avg. total time = %0.5g\n',mean(cputime));
    fprintf('    max sqrt(||A*V - U*S||^2+||A^T*U - V*S||^2) = %0.5g\n',max(PSVD_err));
    fprintf('    max sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n',max(VU_err));
end

% Calling svt_irlba iter times and output the results for paper.
% --------------------------------------------------------------
fprintf(' \n');
fprintf('SVT_IRLBA\n');
for i = 1:iter
    tStart = tic;
    [U,S,V,FLAG] = svt_irlba(A,'m',m,'n',n,'sigma',sigma,'tol',tol,...
                            'psvdmax',psvdmax,'kmax',kmax,'k',k,...
                            'incre',incre,'pwrsvd',pwrsvd,...
                            'display',display,'p0',p0(:,i));

    if FLAG ~= 0
        fprintf('    iter = %d\n',iter);
        fprintf('    FLAG = %d\n',FLAG);
    end
    cputime(i) = toc(tStart);
    sizeS(i)   = size(S,1);
    if ~isempty(S)
       PSVD_err(i) = sqrt(norm(A*V - U*S)^2+norm(A'*U - V*S)^2);
       VU_err(i)   = sqrt(norm(V'*V-eye(size(V,2)))^2+norm(U'*U-eye(size(U,2)))^2);
    else
      PSVD_err(i) = Inf; VU_err(i) = Inf;
    end
end
fprintf('    threshold = %0.5g\n',sigma);
fprintf('    max # singvals = %0.5g\n',max(sizeS));
fprintf('    min # singvals = %0.5g\n',min(sizeS));
fprintf('    avg. total time = %0.5g\n',mean(cputime));
fprintf('    max sqrt(||A*V - U*S||^2+||A^T*U - V*S||^2) = %0.5g\n',max(PSVD_err));
fprintf('    max sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n',max(VU_err));

% Calling svt iter times and output the results for paper.
% --------------------------------------------------------
fprintf(' \n');
fprintf('SVT\n');
for i = 1:iter
    tStart = tic;
    [U,S,V,FLAG] = svt(A,'lambda',sigma,'tol',tol,'k',k,'incre',incre,...
                       'method','deflation');

    if FLAG ~= 0
        fprintf('    iter = %d\n',iter);
        fprintf('    FLAG = %d\n',FLAG);
    end
    cputime(i) = toc(tStart);
    sizeS(i)   = size(S,1);
    if ~isempty(S)
       PSVD_err(i) = sqrt(norm(A*V - U*S)^2+norm(A'*U - V*S)^2);
       VU_err(i)   = sqrt(norm(V'*V-eye(size(V,2)))^2+norm(U'*U-eye(size(U,2)))^2);
    else
      PSVD_err(i) = Inf; VU_err(i) = Inf;
    end
end
fprintf('    threshold = %0.5g\n',sigma);
fprintf('    max # singvals = %0.5g\n',max(sizeS));
fprintf('    min # singvals = %0.5g\n',min(sizeS));
fprintf('    avg. total time = %0.5g\n',mean(cputime));
fprintf('    max sqrt(||A*V - U*S||^2+||A^T*U - V*S||^2) = %0.5g\n',max(PSVD_err));
fprintf('    max sqrt(||V^TV - I||^2+||U^TU-I||^2) = %0.5g\n',max(VU_err));

end % end example41

% Function to download the corresponsing matrix from the SuiteSparse Matrix
% Collection
% -------------------------------------------------------------------------
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
         % urlwrite(address,name_mat); websave(address,name_mat);
         % ------------------------------------------------------------------
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
% -----------------------------------------------------------------------------------
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
