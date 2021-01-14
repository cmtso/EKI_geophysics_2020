function sigma_truth=true_cond(R2_Grid,file,sigma1,sigma2,Grid,L,prior)
x_c=-2.5;
y_c=-2.0;
z_c=-3.0;
r1=1.25;%

x_c2=0;%3;
y_c2=0;
z_c2=-3.0;
r2=0.75;%

% fac1=log(Grid.D(1)*0.65)*ones(Grid.N,1); 
% fac2=log(Grid.D(2)*0.05)*ones(Grid.N,1); 
% fac3=log(Grid.D(3)*0.05)*ones(Grid.N,1); 
% no_fields=2;
% N=Grid.N;
% K=zeros(no_fields*Grid.N,1);
% for ifield=1:no_fields
%     if (ifield==1)
%         pri.len{1}=fac1;
%         pri.len{2}=fac2;
%         pri.len{3}=fac3;
%     else
%         pri.len{1}=fac1;
%         pri.len{2}=fac2;
%         pri.len{3}=fac3;
%     end
%     pri.sigma=prior.K(ifiel   d).per.sigma; pri.nu=prior.K(ifield).per.nu;
%     if ifield==1
%         K(1+(ifield-1)*N:N+(ifield-1)*N,1)=grf3D(Grid, pri,randn(N,1))+log(1e-1)*ones(N,1);
%     else
%         K(1+(ifield-1)*N:N+(ifield-1)*N,1)=grf3D(Grid, pri,randn(N,1))+log(5e-3)*ones(N,1);
%     end
%     
% end
% 
% K_int=K(1:N,1);
% K_mid=K(1+N:2*N,1);
% e_mid=reshape(exp(K_mid),Grid.n(1),Grid.n(2),Grid.n(3));
% e_int=reshape(exp(K_int),Grid.n(1),Grid.n(2),Grid.n(3));
% 
% sigma1 = interp3(Grid.X,Grid.Y,Grid.Z,e_mid,E4D_Grid.x,E4D_Grid.y,E4D_Grid.z,'nearest');    
% sigma2 = interp3(Grid.X,Grid.Y,Grid.Z,e_int,E4D_Grid.x,E4D_Grid.y,E4D_Grid.z);    



%sigma_truth=sigma2+(sigma1-sigma2).*((E4D_Grid.x-x_c2).^2/7.25^2+(E4D_Grid.y-y_c2).^2/0.8^2+(E4D_Grid.z-z_c2).^2/0.3^2<1);
%+(sigma2-sigma1).*((E4D_Grid.x-x_c).^2/1.0^2+(E4D_Grid.y-y_c).^2/3.75^2+(E4D_Grid.z-z_c).^2/0.3^2<1);%...'   

%sigma_truth=sigma1+(sigma2-sigma1).*((  ((R2_Grid.x-x_c2).^2+(R2_Grid.y-y_c2).^2)<1.5^2)&((R2_Grid.z<-1)&(R2_Grid.z>-8))); 
sigma_truth=sigma1+(sigma2-sigma1).*((  ((R2_Grid.x-x_c2).^2+(R2_Grid.y-y_c2).^2)<1.5^2)); 


write_R2_sigma(file,sigma_truth)
%save('Truth','L_mean','RN');