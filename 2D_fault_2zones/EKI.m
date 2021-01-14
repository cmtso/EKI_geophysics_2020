clear all

%%%% Procedure %%%
% -- set up grid and run forward problem in R2
% -- obs. data format should be in that of protocol.dat for R2

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

load('Results_11.mat') % change iter # if needed

a = 1./unique(sigma_mean)
for i = 1:size(sigma,1)
    zone1_prob(i) = sum(sigma(i,:)== 1/a(1)) ;
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
axis equal
xlim([0 48])
ylim([-10 0])
xlabel('$X$','Interpreter','latex','FontSize',16)
ylabel('$Z$','Interpreter','latex','FontSize',16)
%caxis([-1 2])
hold on, plot([0 24 24 48],[-7 -7 -5 -5],'r-','LineWidth',3), hold off
title('Mean est. resistivities are 243 and 2362 Ohm m')


mycolormap = customcolormap(linspace(0,1,4), {'#fbeed7','#ffba5a','#ff7657','#665c84'});
colormap(mycolormap);

%EKI.m

% for i = 1:size(res0,1)
%     if(res0(i,1) > -2.5 && res0(i,1) < 0.5 && res0(i,2) > -3.5 && res0(i,2) < 0.5)
%         res0(i,3) = 1;
%     end
% end
% res0(:,4) = log10(res0(:,3));

%scatter(R2_Grid.x,R2_Grid.y,[],1./sigma_mean)
%%
res = [R2_Grid.x R2_Grid.y];
res(:,3) = 2500;
for i = 1:size(res,1)
    if (res(i,2) > -5)
        res(i,3) = 250;
    end
    if (res(i,2) > -7 && res(i,1) < 24)
        res(i,3) = 250;
    end
end
plot(res(:,3))
res(:,4) = log10(res(:,3));