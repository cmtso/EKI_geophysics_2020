clear all

out_file=strcat('Results');

seed=1243;
rng(seed)


node_file = 'sp20K/sp.1.node';
elem_file = 'sp20K/sp.1.ele';
trn_file = 'sp20K/sp.trn' ; % see any nodal points are translated


%%%%%%%%%%%%%%%%%%%%%%%%%%% define data for inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%
noise=0.02; %% percentage of noise added to true data
data=get_E4D_data('sp20K/baseline.sig.srv');  %%get data from e4d for the simulation with the true conductivity
M=length(data);
noise_data=noise*data; %this is probably not the bests because of the scaling (we need to use a sensible choice here)
data=data+noise_data.*randn(M,1); %%add noise to the data
inv_sqrt_C=diag(1./noise_data); %inverse of square root of measurement error covariance


%%%%%%%%%%%%%%%%%%%%%%%%%%% define initial ensemble %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%read elements' coordinates
ec=read_elements(trn_file,elem_file,node_file);

N_En=100;  %ensemble size
OneN=ones(N_En);
OneN=(1/N_En)*OneN;


%%define conductivity values (these must be consisten with the ones from
%%the truth)
R1=1e-1;
R2=1e-3;

%%%the function below generates an initial ensemble on a regular grid via fft and
%%interpolates it on e4d grid

L=[1e3,1e3,5e2]; %length of the domain
n=[50,50,30];
U=get_ensemble(L,n,ec,N_En); %%ensemble of level set functions (gaussian fields) on e4d elements
U_m=U*OneN; %compute mean

%compute ensemble of conductivities
for en=1:N_En
    sigma(:,en)=R1+(R1-R2).*(U(:,en)<=0);
end


Cond=0;
Z=zeros(M,N_En);
iter=0;
t(1)=0;

%%%EnKI%%%%%%%%%%%%%%%%%%%%%
while (Cond==0)
    iter=iter+1;
    for en=1:N_En
        
        %%% we need a file to save sigma(:,en) in the conductivity file for e4d
        
        
        %%% we need command to run e4d here for this new sigma
        
        
        %%% We now call output from e4d for the run above
        
        Output(:,en)=get_E4D_data('sp20K/baseline.sig.srv');  %%output
        
        
        %%data misfit (weighted by inverse of sqrt of measurement covariance):
        Z(:,en)=inv_sqrt_C*(data-Output(:,en));
    end
    Z_m=Z*OneN;  %compute mean of misfits
    Misfit(iter)=norm(Z_m(:,1)); %this is mean data misfit
    
    %%EKI code to update U
    Delta_Z=Z-Z_m; %compute deviations for data misfit
    C=1/(N_En-1)* Delta_Z*Delta_Z'; %covariance me data misfit
    Mzz=1/(N_En)* (Z)*(Z)';  %%stuff I need to compute alpha
    alpha=sum(diag(Mzz))/M;  %compute alpha
    %%note that the above should be nothing but the average of the squared data
    %%misfits divided by the number of observations
    if (t(iter)+1/alpha>1)
        alpha=1/(1-t(iter));
        disp('EnKI converged') %the sum of 1/alpha's should be 1 for convergence
        Cond=1;
    end
    t(iter+1)=t(iter)+1/alpha;
    Delta_U=U-U_m;
    C_u_z=1/(N_En-1)*Delta_U*Delta_Z'; %Cross covariance between data (misfit) and unknown level set
    
    E=sqrt(alpha)*randn(M,N_En);  %perturb observations
    meanE=E*OneN;
    E=E-meanE;
    Z=Z+E;
    
    rhs=(C+alpha*eye(M))\Z;  %note that I'm working with data misfit (not data).
    U=U-C_u_z*rhs;  %update ensemble of level-sets
    U_m=U*OneN;
    %%compute ensemble of updated conductivities (via level function):
    for en=1:N_En
        sigma(:,en)=R1+(R1-R2).*(U(:,en)<=0);
    end
    %%compute conductivity of the mean of ensemble of level-set functions
    sigma_m1=R1+(R1-R2).*(U_m(:,1)<=0);
    %%compute mean of ensemble of conductivities
    sigma_m=sigma*OneN;
    %%note sigma_m and sigma_m1 are not the same
    save(strcat(out_file, '_', num2str(iter)),  'iter',...'
        'alpha','t', 'Misfit','U_m','U','Z','sigma','sigma_m','sigma_m1')
end










function data=get_E4D_data(srv_file)

fid = fopen(srv_file,'r');
nelec = str2double(fgetl(fid)) ;
fscanf(fid,'%f',[5 nelec])';
fgetl(fid);
fgetl(fid);

nd = str2double(fgetl(fid)) ;
data = fscanf(fid,'%f',[7 nd])';
data = data(:,6);

fclose(fid);
end



% returns ec (ne x 3) for element centre that corresponds ??.sig files

%%
function ec=read_elements(trn_file,elem_file,node_file)
trn = dlmread(trn_file) ;

fid = fopen(elem_file,'r');
ne = fgetl(fid); ne = sscanf(ne,'%d'); ne = ne(1);
elem = fscanf(fid, '%d',[6  ne])'; % col 2-5 are node numbers
fclose(fid);

fid = fopen(node_file,'r');
nn = fgetl(fid); nn = sscanf(nn,'%d'); nn = nn(1);
node = fscanf(fid, '%f',[6  nn])'; % col 2-4 are node coordinates
fclose(fid);
node(:,2:4) = node(:,2:4) + trn ; % add the translation factor back to the node coordinates.

% calculate ec
ec = zeros(ne,3);
for i = 1:ne
    ec(i,:) =  mean( node(elem(i,2:5),2:4) ) ;
end
end

function U=get_ensemble(L,n,ec,N_en)

N=prod(n);
[X,Y,Z]=meshgrid(linspace(-L(1)/2,L(1)/2,n(1)),linspace(-L(2)/2,L(2)/2,n(2)),linspace(-L(3),0,n(3)));
x=reshape(X,N,1);
y=reshape(Y,N,1);
z=reshape(Z,N,1);
for en=1:N_en
    Lambda=L(1)/20; % we need to play with this
    GRF_f= gaussrnd_3D(6,1e1,Lambda,n,L,randn(N,1)); %%generate Gaussian fiel in Fourier space
    GRF=reshape(idctn(reshape(GRF_f,n)),N,1);%% Gaussian fiel in physical space
    U(:,en)= griddata(x,y,z,GRF,ec(:,1),ec(:,2),ec(:,3));%%interpolate to e4d grid
end


    function L = gaussrnd_3D(s,q,l,n,L,xi)
        
        xi = reshape(xi, n);
        [K1,K2,K3] = meshgrid(0:n(2)-1,0:n(1)-1,0:n(3)-1);
        fac=l^2;
        coef = (pi^2*(K1.^2/L(1)^2+K2.^2/L(2)^2+K3.^2/L(3)^2)*fac + 1).^(-s/2);
        L =sqrt(q)*coef.*xi;
        L(1,1,1)=0;
        L = reshape(L,prod(n),1);
        
    end
end
%trn = dlmread('baseline.sig') ;

