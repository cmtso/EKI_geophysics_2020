function Un=update_unkown(M,N_En,Un,Un_m,Z,Delta_Z,alpha)

C=1/(N_En-1)* (Delta_Z)*(Delta_Z)'; %covariance me data misfit
E=sqrt(alpha)*randn(M,N_En);  %perturb observations
E=E-mean(E,2);
Z=Z+E;
B=(C+alpha*eye(M))\Z;
for i=1:length(Un)
    Delta_U=Un{1,i}-Un_m{1,i};
    C_u_z=1/(N_En-1)*Delta_U*Delta_Z'; %Cross covariance between data (misfit) and unknown level set
    Un{1,i}=Un{1,i}-C_u_z*B;
end

