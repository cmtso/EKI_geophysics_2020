clear all

%%%% Procedure %%%
% -- set up grid and run forward problem in R2
% -- obs. data format should be in that of protocol.dat for R2

addpath('Tools')
nseed=1;
rng(nseed*129)
%R2_path = pwd ; %'C:\Users\mtso\Downloads\pyr2-master\src\resipy\exe' ;


cond_file='resistivity.dat'; 


%%Write a routine to read grid from Andy's code and change the line below.
R2_Grid=get_R2_Grid('forward_model.dat'); % mustn't be cropped, use f001_res.dat or forward_model.dat
%#R2_Grid.trn = [-4 -7.5] ; % *** make this varargin*** (x,y) coordinate translation so that by subtracting it the center of the inversion domain is at (0,0) , temp, MT 190828

%%define conductivity values (these must be consisten with the ones from
%%the truth)
% note: R2 write files as resistivities
sigma2= 1;
sigma1= 1/100;


L=[16,16]; %dimensions of the 2D domain where we wish to recover conductivity 
%%for this choice of L 16x12 the domain is [-8, 8]x [-6,6]. Change if needed in Set_Grid. 


n=[64,64]; %number of cells for the discretisation of L (for random fields) 

Grid=Set_Grid(n,L);

option=0; %option=1 for variable lengthscale and option=0 for constant lengthscale
%%test code with option=0 first

% %%% Use Grid from Andy's code to define a true conductivity
% sigma_truth=true_cond(R2_Grid,'truth.dat',sigma1,sigma2,Grid,L,Pr);
% 
% %%run Andy's code with this conductivity
% %if (ispc == 1), system('R2.exe'); else system('wine R2.exe'); end
% % system(cstrcat ('wine ',R2_path,'/R2.exe'))  % Octave
load Prior
sigma_prior=sigma;
sigma_truth = dlmread('forward_model.dat'); sigma_truth = 1./ sigma_truth(:,3);
load Results_11.mat
figure(1)
subplot(1,3,1)
scatter(R2_Grid.x,R2_Grid.y,200,sigma_truth,'filled')
axis([-8,8,-6,6])
title('truth')
subplot(1,3,2)
scatter(R2_Grid.x,R2_Grid.y,200,sigma_mean,'filled')
axis([-8,8,-6,6])
title('sigma(E(u))')
subplot(1,3,3)
scatter(R2_Grid.x,R2_Grid.y,200,mean(sigma,2),'filled')
axis([-8,8,-6,6])
title('E(sigma(u))')
figure(2)
subplot(1,2,1)
scatter(R2_Grid.x,R2_Grid.y,200,var(sigma_prior,0,2),'filled');colormap jet
axis([-8,8,-6,6])
caxis([0,0.22])
colorbar
title('variance prior')
subplot(1,2,2)
scatter(R2_Grid.x,R2_Grid.y,200,var(sigma,0,2),'filled');colormap jet
axis([-8,8,-6,6])
caxis([0,0.22])
colorbar
title('variance posterior')

