clear all

%%%% Procedure %%%
% -- set up grid and run forward problem in R2
% -- obs. data format should be in that of protocol.dat for R2

% need this to get reshaping right

addpath('Tools')
nseed=1;
rng(nseed*129)
%R2_path = pwd ; %'C:\Users\mtso\Downloads\pyr2-master\src\resipy\exe' ;


cond_file='resistivity.dat'; 

trn = [24 -10]; %(+ve x to move domain right,-ve y to move domain up,)
%%Write a routine to read grid from Andy's code and change the line below.
R2_Grid=get_R2_Grid('forward_model.dat',trn); % mustn't be cropped, use f001_res.dat or forward_model.dat

%%define conductivity values (these must be consisten with the ones from
%%the truth)
% note: R2 write files as resistivities

%%DON't NEED THIS ANYMORE (SEE Set_prior.m) WARNING!!!
sigma2= 1/250 %10; WARNING!!!
sigma1= 1/2500 %1/100; WARNING!!!


L=[80,20]; %dimensions of the 2D domain where we wish to recover conductivity 
%%for this choice of L 16x12 the domain is [-8, 8]x [-6,6]. Change if needed in Set_Grid. 


%n=[480*2,200]; %number of cells for the discretisation of L (for random fields) 
n=[80,20].*4;

plot(R2_Grid.x,R2_Grid.y)

Grid=Set_Grid(n,L);

option=0; %option=1 for variable lengthscale and option=0 for constant lengthscale
%%test code with option=0 first
n_fields=2;  %number of fields 2 or 3 

Pr=Set_prior(Grid,sigma1,sigma2,option,n_fields);

sigma_truth = dlmread('forward_model.dat'); sigma_truth = 1./ sigma_truth(:,3);

%%generate synthetic data
%cd .. % 
noise=0.02;%0.01; %% percentage of noise added to true data


%%change this routine to read output (i.e. voltages) from Andy's code
data=get_R2_data('R2_forward.dat');

%%get data from e4d for the simulation with the true conductivity
Data.data_noise_free=data;
noise_data1=noise*data;
noise_data2=abs(max(data))*1e-5;%(max(abs(data))-min(abs(data)))*1e-4;
Data.data=data+noise_data1.*randn(length(data),1)+noise_data2.*randn(length(data),1); %%add two components of noise to the data
Data.inv_sqrt_C=diag(1./sqrt(noise_data1.^2+noise_data2.^2)); %inverse of square root of measurement error covariance
save('Data','Data')
%%%%%%%%%%%%%%%%%%We need this for the prior%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%****************************Prior Definition %****************************


N_En=300;  %ensemble size

Un=Get_prior(Grid,Pr,N_En); 
save('Un','Un')
load('Un.mat')

out_file=strcat('Results');
tuning=30;
sigma_mean=Inversion(R2_Grid,Grid,N_En,Pr,Un,Data,out_file,cond_file,tuning);


%**************************** Post processing %****************************
% read vtk file (get mesh info), then add results to the vtk. struct

clear vtk

load('Results.mat') % change iter # if needed

a = 1./unique(sigma_mean);
for i = 1:size(sigma,1)
    zone1_prob(i) = sum(sigma(i,:) < 1e-3) ;
    var_log_res(i) = var(log10(1./sigma(i,:)));

end

vtk = read_vtk() ; 
% ##vtk = add_vtk_scalar(vtk,{"Mean resistivity","Zone 1 probability"}, ...
% ##          [1./sigma_mean (zone1_prob/300)'] ) ; % add scalar variable to the vtk struct

vtk.scalar_data = [vtk.scalar_data log10(1./sigma_mean) zone1_prob'./300 var_log_res'] ;
add_list = {"Log mean resistivity","Zone 1 probability","var(log res)"} ;
vtk.scalar_list(end+1:end+numel(add_list)) = add_list; 

%vtk.polyline = dlmread(fullfile(vtk.folder, 'polyline.txt'));

fieldnames(vtk)

plot_vtk_2D()  % will show drop down menu to let you selct variable

%%% optional
%rectangle('Position',[14 -4 2 3],'edgecolor','r') % 1-4, 14-16
%rectangle('Position',[0 -16 48 16])
%get(gca,'Position')
%0.1191    0.1347    0.7098    0.7903
axis equal
xlim([0 48])
ylim([-10 0])
xlabel('$X$','Interpreter','latex','FontSize',16)
ylabel('$Z$','Interpreter','latex','FontSize',16)
%caxis([-1 2])
hold on, plot([0 24 24 48],[-7 -7 -5 -5],'y-','LineWidth',2), hold off %2 y for mean
title('Mean est. resistivities are 243 and 2362 Ohm m')

colormap(brewermap(21,'Spectral'))
set(gca,'clim',[2.0610    3.8678])
%% generate true fields

val1 = 5; val2 = 1;
pri.len{1}=log(val1)*ones(Grid.N,1); %%this controls the (log) lengthscale in x
pri.len{2}=log(val2)*ones(Grid.N,1); %%this controls the (log) lengthscale in y

%%for this use val1 and val2 proportional to the size of the domain on
%%each direction (e.g. val1=0.1*(length in x))

pri.sigma=1.0; %%control variance
pri.nu=1.0; %%control smoothness (use 1, 2 or 3; value 1 is more realistic for physical properties)
me = [0]; %% mean value of each zone
Field= me +grf2D_fields(Grid, pri,randn(Grid.N,1)); 
temp = reshape(Grid.X,prod(n),1);
temp(:,2) = reshape(Grid.Y,prod(n),1);
imagesc(temp(:,1),temp(:,2),reshape((Field(:,1)),n)')
colorbar

zones = 2*ones( prod(n),1 );

for i=1:size(temp,1)
    if (temp(i,1) < 0)
        if (temp(i,2) < 3), zones(i) = 1; end
    else
        if (temp(i,2) < 5), zones(i) = 1; end
    end
end
scatter(temp(:,1),temp(:,2),[],zones)

v = reshape(reshape(Field,n)',prod(n),1) ;
scatter(temp(:,1),temp(:,2),[],(v))

v(zones == 1) = (1/4)*v(zones == 1) + log10(2500); 
v(zones == 2) = (1/8)*v(zones == 2) + log10(250);
vq = griddata(temp(:,1),temp(:,2),(v),R2_Grid.x,R2_Grid.y);
vq(isnan(vq)) = log10(2500); %BACKGROUND
scatter(R2_Grid.x,R2_Grid.y,[],vq)
%%
clear vtk
vtk = read_vtk() ; 
vtk.scalar_data = [vtk.scalar_data vq] ;
vtk.scalar_list(end+1:end+1) = {"True log resistivity"}; 
%vtk.polyline = dlmread(fullfile(vtk.folder, 'polyline.txt'));
plot_vtk_2D()  % will show drop down menu to let you selct variable

%%% optional
%rectangle('Position',[14 -4 2 3],'edgecolor','r') % 1-4, 14-16
%rectangle('Position',[0 -16 48 16])
axis equal
xlim([0 48])
ylim([-10 0])
xlabel('$X$','Interpreter','latex','FontSize',16)
ylabel('$Z$','Interpreter','latex','FontSize',16)
%caxis([-1 2])
hold on, plot([0 24 24 48],[-7 -7 -5 -5],'r-','LineWidth',3), hold off
title('')




dlmwrite('true_res.dat',[zeros(numel(vq),2) 10.^vq vq])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from Marco 200327

Hi Michael, 
I looked at the folder you mentioned and the subfolder Tools is missing a few files. So I went to the folder "2D" and modified in "Tools" the files "Get_prior" and "physical". Also from "Michael_2D_200325/Tools" I modified "Set_prior" and added "grf2D_fields".
In "Set_prior" you can modify the variance of the prior fields via (for the case of 3 regions):

K(1).per.sigma =0.05;
K(2).per.sigma =0.05;
K(3).per.sigma =0.05;

There were quite a few modifications needed! I hope it works!

I guess you also need a file to generate a truth with variable properties as well? The code below can be used to produce one random field. You can produce two or tree of these for the truth in  your current examples; you just need to play with sigma to make sure is small and the mean (me) so that it has the mean value on each region. Let me know if any of this makes sense?
Have a nice weekend!
M

    pri.len{1}=log(val1)*ones(Grid.N,1); %%this controls the (log) lengthscale in x

    pri.len{2}=log(val2)*ones(Grid.N,1); %%this controls the (log) lengthscale in y

    %%for this use val1 and val2 proportional to the size of the domain on

    %%each direction (e.g. val1=0.1*(length in x))

 

    pri.sigma=0.01; %%control variance

    pri.nu=1.0; %%control smoothness (use 1, 2 or 3; value 1 is more realistic for physical properties)

    Field= me+grf2D_fields(Grid, pri,randn(Grid.N,1));    
