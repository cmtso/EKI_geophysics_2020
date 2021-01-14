clear all

addpath('Tools')

cond_file='baseline.sig'; 
node_file = 'sp.1.node';
elem_file = 'sp.1.ele';
trn_file = 'sp.trn' ; % see any nodal points are translated

E4D_Grid=get_E4D_Grid(trn_file,elem_file,node_file);


%%define conductivity values (these must be consisten with the ones from
%%the truth)
sigma2=1e-1;
sigma1=1e-3;
%%define true conductivity anc save it to file
%L=[20,16,15]; %length of the domain
L=[16,12,10]; %length of the domain
n=[30,30,30];%L=[1e3,1e3,500]; %length of the domain


Grid=Set_Grid(n,L);

load Results3
visualise(Grid.n,E4D_Grid.ec,sigma_mean,L)
write_sigma('baseline_estimate_3.sig',sigma_mean)
%visualise(Grid.n,E4D_Grid.ec,sigma(:,200),L)
tempo=Un_m{1,2};
figure
field1=reshape(exp(Un_m{1,3}+tempo(1))/Grid.D(1),Grid.n(1),Grid.n(2),Grid.n(3));
U_temp=zeros(Grid.n(1),Grid.n(2),Grid.n(3));
for iz=1:Grid.n(3)
    U_temp(:,:,iz)=field1(:,:,iz)';
end
visualise_regular_grid(Grid,U_temp)

figure
field1=reshape(exp(Un_m{1,4}+tempo(2))/Grid.D(2),Grid.n(1),Grid.n(2),Grid.n(3));
U_temp=zeros(Grid.n(1),Grid.n(2),Grid.n(3));
for iz=1:Grid.n(3)
    U_temp(:,:,iz)=field1(:,:,iz)';
end
visualise_regular_grid(Grid,U_temp)

figure
field1=reshape(exp(Un_m{1,5}+tempo(3))/Grid.D(3),Grid.n(1),Grid.n(2),Grid.n(3));
U_temp=zeros(Grid.n(1),Grid.n(2),Grid.n(3));
for iz=1:Grid.n(3)
    U_temp(:,:,iz)=field1(:,:,iz)';
end
visualise_regular_grid(Grid,U_temp)



