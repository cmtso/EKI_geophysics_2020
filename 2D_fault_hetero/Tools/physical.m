function [sigma, U]=physical(Grid,R2_Grid,Pr,Un)

RN=Un{1,1};
L_mean=Un{1,2};
L_per.x=Un{1,3};
L_per.y=Un{1,4};
fields_mean=Un{1,5};
fields_per=Un{1,6};
fields_le=Un{1,7};
grf=Get_level_set(Grid,Pr,L_mean,L_per,RN);
U_temp=grf';
U = interp2(Grid.X,Grid.Y,U_temp,R2_Grid.x,R2_Grid.y,'spline',0);
%%compute physical conductivity
%sigma=Pr.sigma1+(Pr.sigma2-Pr.sigma1).*(U<=-0.25);


NF=Pr.n_fields;
cut_off=0.1;

cut=linspace(-cut_off,cut_off, NF-1);

K=Get_Fields(fields_mean,fields_per,NF,Grid,fields_le,Pr,R2_Grid);%,Model,prior,field_le,Inv);


sigma=exp(K{1}).*(U<=cut(1));
for nf=2:NF-1
    sigma=sigma+exp(K{nf}).*((U>cut(nf-1))&(U<=cut(nf)));
end
sigma=sigma+exp(K{end}).*(U>cut(end));
end
function K=Get_Fields(fields_mean,fields_per,NF,Grid,fields_le,prior,R2_Grid)%,Model,prior,field_le,Inv)
K=cell(NF,1);
%vargout=cell(NF,1);

for nf=1:NF
    K{nf}=fields_mean(nf);
    %     vargout{nf}=fields_mean(nf);
    for idim=1:Grid.dim
        pri.len{idim}=fields_le(idim+(nf-1)*NF)*ones(Grid.N,1);
    end
    pri.sigma=prior.K(nf).per.sigma; pri.nu=prior.K(nf).per.nu;
    K_temp=K{nf}+grf2D_fields(Grid, pri,fields_per(1+(nf-1)*Grid.N:Grid.N+(nf-1)*Grid.N,1));    
    K{nf}=interp2(Grid.X,Grid.Y,reshape(K_temp,Grid.n(1),Grid.n(2))',R2_Grid.x,R2_Grid.y,'spline',mean(K_temp)); % 'spline' value was 0 
%    K{nf}=interp2(Grid.X,Grid.Y,reshape(K_temp,Grid.n(1),Grid.n(2))',R2_Grid.x,R2_Grid.y,'nearest');    
end
end