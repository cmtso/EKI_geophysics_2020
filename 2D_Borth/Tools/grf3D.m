function P=grf3D(Grid, prior,xi,fac)


lx=prior.len{1};
ly=prior.len{2};
lz=prior.len{3};
sigma=prior.sigma;
nu=prior.nu;
dim=Grid.dim;

alpha2=sigma^2*2^dim*pi^(dim/2)*gamma(nu+dim/2)/gamma(nu);


Nx=Grid.n(1); Ny=Grid.n(2); Nz=Grid.n(3); N=Nx*Ny*Nz;
hx=Grid.D(1)/Nx; hy=Grid.D(2)/Ny; hz=Grid.D(3)/Nz;

Lambda_x=exp(lx);
Lambda_y=exp(ly);
Lambda_z=exp(lz);

f=alpha2*Lambda_x.*Lambda_y.*Lambda_z*hx*hy*hz;
sqrt_f=sqrt(alpha2*Lambda_x.*Lambda_y.*Lambda_z)*hx*hy*hz;
A=get_Mat3D(Lambda_x,Lambda_y,Lambda_z,Grid,fac);

if (nu==0.5)
    P=A\(sqrt(f).*xi);   
elseif (nu==1.5)
    [R,p]=chol(A);
    y=R\xi;
    P=A\(sqrt_f.*y);    
elseif (nu==2.5)
    y=A\(sqrt(f).*xi);    
    P=A\(hx*hy*hz.*y);    
end

function A=get_Mat3D(Lambda_x,Lambda_y,Lambda_z,Grid,fac)
Nx=Grid.n(1); Ny=Grid.n(2); Nz=Grid.n(3); N=Nx*Ny*Nz;
hx=Grid.D(1)/Nx; hy=Grid.D(2)/Ny; hz=Grid.D(3)/Nz;

lam_x=reshape(Lambda_x,Nx,Ny,Nz);
lam_y=reshape(Lambda_y,Nx,Ny,Nz);
lam_z=reshape(Lambda_z,Nx,Ny,Nz);
Lx=lam_x.^2;
Ly=lam_y.^2;
Lz=lam_z.^2;
%fac=0;
lambda_left=fac*lam_x(1,:,:);
lambda_right=fac*lam_x(Nx,:,:);
lambda_top=fac*lam_y(:,Ny,:);
lambda_bottom=fac*lam_y(:,1,:);
lambda_top2=fac*lam_z(:,:,Nz);
lambda_bottom2=fac*lam_z(:,:,1);


tx=hy*hz/hx; TX=zeros(Nx+1,Ny,Nz);
ty=hx*hz/hy; TY=zeros(Nx,Ny+1,Nz);
tz=hx*hy/hz; TZ=zeros(Nx,Ny,Nz+1);
Average_x=0.5*(Lx(1:Nx-1,:,:)+Lx(2:Nx,:,:));
TX(2:Nx,:,:)=Average_x.*tx;
Average_y=0.5*(Ly(:,1:Ny-1,:)+Ly(:,2:Ny,:));
TY(:,2:Ny,:)=Average_y.*ty;
Average_z=0.5*(Lz(:,:,1:Nz-1)+Lz(:,:,2:Nz));
TZ(:,:,2:Nz)=Average_z.*tz;
TX2=TX;
TY2=TY;
TZ2=TZ;
TX(Nx+1,1:Ny,1:Nz)=(tx./(1/2+lambda_left/hx)).*Lx(1,:,:);
TX(1,1:Ny,1:Nz)=(tx./(1/2+lambda_right/hx)).*Lx(Nx,:,:);
TY(1:Nx,1,1:Nz)=(ty./(1/2+lambda_bottom/hy)).*Ly(:,1,:);
TY(1:Nx,Ny+1,1:Nz)=(ty./(1/2+lambda_top/hy)).*Ly(:,Ny,:);
TZ(1:Nx,1:Ny,1)=(tz./(1/2+lambda_bottom2/hz)).*Lz(:,:,1);
TZ(1:Nx,1:Ny,Nz+1)=(tz./(1/2+lambda_top2/hz)).*Lz(:,:,Nz);


x1=reshape(TX(1:Nx,:,:),N,1); x2=reshape(TX(2:Nx+1,:,:),N,1);
y1=reshape(TY(:,1:Ny,:),N,1); y2=reshape(TY(:,2:Ny+1,:),N,1);
z1=reshape(TZ(:,:,1:Nz),N,1); z2=reshape(TZ(:,:,2:Nz+1),N,1);

x12=reshape(TX2(1:Nx,:,:),N,1); x22=reshape(TX2(2:Nx+1,:,:),N,1);
y12=reshape(TY2(:,1:Ny,:),N,1); y22=reshape(TY2(:,2:Ny+1,:),N,1);
z12=reshape(TZ2(:,:,1:Nz),N,1); z22=reshape(TZ2(:,:,2:Nz+1),N,1);

%z1=reshape(TZ(:,:,1:Nz),N,1); z2=reshape(TZ(:,:,2:Nz+1),N,1);

DiagVecs=[-z22,-y22,-x22,x1+x2+y1+y2+z1+z2,-x12,-y12,-z12];

DiagIndx=[-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
A=spdiags(DiagVecs,DiagIndx,N,N);

DiagVecs2=[ones(N,1)*hx*hy*hz];
DiagIndx2=[0];
C=spdiags(DiagVecs2,DiagIndx2,N,N);




%RHS=sqrt(alpha2*Lambda_x.*Lambda_y.*Lambda_z*hx*hy*hz).*xi;
%RHS=xi*hx*hy;
A=(C+A);



