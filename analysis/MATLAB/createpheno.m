% Simulate a quantitative trait as a linear combination of SNP genotypes and
% population factors (e.g., estimated using principal components) plus
% Gaussian noise. The script parameters r and d control proportion of
% variance in the phenotype explained by genotypes and population factors.
%
% The population factors are stored in matrix Z, and this can be adjusted
% as needed.
clear

% SCRIPT PARAMETERS
% -----------------
r  = 0.5;  % Prop. variance explained by genetic & population factors.
d  = 0.3;  % Proportion of genetic variance due to QTLs (causal SNPs).
na = 100;  % Number of QTLs (causal SNPs).

% LOAD GENOTYPE DATA
% ------------------
fprintf('Loading genotype data.\n');
load('../../data/hap550.mat');

% LOAD PCs
% --------
fprintf('Loading PC data.\n');
load('../../data/hap550_pc.mat');

% This matrix gives the population factors. In this example, we use the
% first principal component to capture population structure.
Z = pc(:,1);

% Initialize the sequence of pseudorandom numbers.
rng(1);

% DATA PRE-PROCESSING STEPS
% -------------------------
% Center the population covariates (columns of Z).
m = size(Z,2);
Z = Z - repmat(mean(Z),m,1);

% Center the columns of the genotype matrix so that each column has a mean
% of zero.
fprintf('Centering columns of genotype matrix.\n');
p = length(pos);
for i = 1:p
  X(:,i) = X(:,i) - mean(X(:,i));
end

% SIMULATE PHENOTYPE ATA
% ----------------------
fprintf('Simulating phenotype data.\n');

% Here we assume the population factors all have the same effect and in
% the same direction.
u = ones(m,1);

% Generate effects for the causal SNPs; the effects for the non-causal SNPs
% are all set to zero.
i    = randperm(p);
i    = i(1:na);
b    = zeros(p,1);
b(i) = randn(na,1);
  
% Adjust the additive effects of the population factors and the SNPs so that
% we control for (1) the proportion of additive genetic variance that is due
% to SNPs ("d"), and (2) the proportion of variance explained by all genetic
% & population factors ("r"). That is, we adjust vectors b and u so that
%
%   r = a/(a+1)
%   d = c/a,
%
% where I've defined 
%
%   a = [u; b]'*cov([Z X])*[u; b],
%   c = b'*cov(X)*b.
%
if d == 0 | d == 1
  error('d cannot be exactly 0 or exactly 1');
end
n  = length(id);
st = r/(1-r) * d/var(X*b,1);
b  = sqrt(st) * b;
zu = Z*u;
xb = X*b;
sa = max(roots2(var(zu,1),...
                2*dot(xb,zu)/n,...
                var(xb,1) - r/(1-r)))^2;
u  = sqrt(sa) * u;
clear zu xb

% Generate the quantitative trait measurements.
y = Z*u + X*b + randn(n,1);

% Center the phenotype.
y = y - mean(y);

% SAVE PHENOTYPE DATA TO FILE
% ---------------------------
fprintf('Saving phenotype data.\n')
simparams = struct('r',r,'d',d,'na',na,'u',u,'b',b);
save('hap550_pheno.mat','y','simparams','-v7.3');
