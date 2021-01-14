function [alpha, misfit]=compute_alpha(Z,N_En,M)

Mzz=1/(N_En)* (Z)*(Z)';  %%stuff I need to compute alpha
misfit=sum(diag(Mzz));
alpha=misfit/M;  %compute alpha

