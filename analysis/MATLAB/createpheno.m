% TO DO: Explain here what this script does, and how to use it.
clear

% SCRIPT PARAMETERS
% -----------------
r   = 0.5;  % Prop. variance explained by genetic & population factors.
rho = 0.3;  % Proportion of genetic variance due to QTLs (causal SNPs).
na  = 100;  % Number of QTLs (causal SNPs).

% LOAD GENOTYPE DATA
% ------------------
fprintf('Loading genotype data.\n');
load('../../data/hap550.mat');

% LOAD PCs
% --------
fprintf('Loading PC data.\n');
load('../../data/hap550_pc.mat');

% Initialize the sequence of pseudorandom numbers.
rng(1);

% DATA PRE-PROCESSING STEPS
% -------------------------
% Center the columns of the genotype matrix so that each column has a mean
% of zero.
fprintf('Centering columns of genotype matrix.\n');
p = length(pos);
for i = 1:p
  X(:,i) = X(:,i) - mean(Xsub(:,i));
end

  % Generate (small) polygenic additive effects for the SNPs.
  u = randn(p,1);

  % Generate (large) QTL effects for the SNPs.
  is       = randperm(p);
  is       = is(1:na);
  beta     = zeros(p,1);
  beta(is) = randn(na,1);
  
  % Adjust the additive effects so that we control for (1) the proportion of
  % additive genetic variance that is due to QTL effects (D), and (2) the
  % proportion of variance explained (R). That is, we adjust BETA and U so
  % that
  %
  %   R = A/(A+1)
  %   D = B/A,
  %
  % where I've defined 
  %
  %   A = (U + BETA)' * COV(X) * (U + BETA),
  %   B = BETA' * COV(X) * BETA.
  %
  st   = r/(1-r) * d/var(X*beta,1);
  beta = sqrt(st) * beta;
  switch d
   case 0
    sa = r/(1-r)/var(X*u,1);
   case 1
    sa = 0;
   otherwise
    sa = max(roots2(var(X*u,1),2*dot(X*beta,X*u)/n,var(X*beta,1) - r/(1-r)))^2;
  end
  u = sqrt(sa) * u;
