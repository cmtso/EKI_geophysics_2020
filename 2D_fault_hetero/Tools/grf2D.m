function U=grf2D(Grid, prior,xi)

Nx=Grid.n(1);
Ny=Grid.n(2);
lx=prior.len{1};
ly=prior.len{2};
sigma=prior.sigma;
nu=prior.nu;
dim=Grid.dim;

alpha2=sigma^2*2^dim*pi^(dim/2)*gamma(nu+dim/2)/gamma(nu);
    
Lambda_x=exp(lx);
Lambda_y=exp(ly);

hx=Grid.D(1)/Nx;
hy=Grid.D(2)/Ny;

f=alpha2*Lambda_x.*Lambda_y*hx*hy;
sqrt_f=sqrt(alpha2*Lambda_x.*Lambda_y)*hx*hy;
A=get_Mats(Grid,Lambda_x,Lambda_y);

if (nu==1)
    U=A\(sqrt(f).*xi);   
elseif (nu==2)
    [R,p]=chol(A);
    y=R\xi;
    U=A\(sqrt_f.*y);    
elseif (nu==3)
    y=A\(sqrt(f).*xi);    
    U=A\(hx*hy.*y);    
end
end

function A=get_Mats(Grid,Lambda_x,Lambda_y)

Nx=Grid.n(1); Ny=Grid.n(2); N=Nx*Ny;
hx=Grid.D(1)/Nx; hy=Grid.D(2)/Ny;

lam_x=reshape(Lambda_x,Nx,Ny);
lam_y=reshape(Lambda_y,Nx,Ny);


Lx=lam_x.^2;
Ly=lam_y.^2;

fac=0*1.2;
lambda_left=fac*lam_x(1,:);
lambda_right=fac*lam_x(Nx,:);
lambda_top=fac*lam_y(:,Ny);
lambda_bottom=fac*lam_y(:,1);


tx=hy/hx; TX=zeros(Nx+1,Ny);
ty=hx/hy; TY=zeros(Nx,Ny+1);
Average_x=0.5*(Lx(1:Nx-1,:)+Lx(2:Nx,:));
TX(2:Nx,:)=Average_x.*tx;
Average_y=0.5*(Ly(:,1:Ny-1)+Ly(:,2:Ny));
TY(:,2:Ny)=Average_y.*ty;
TX2=TX;
TY2=TY;
% 
TX(Nx+1,1:Ny)=(tx./(1/2+lambda_left/hx)).*Lx(1,:);
TX(1,1:Ny)=(tx./(1/2+lambda_right/hx)).*Lx(Nx,:);
TY(1:Nx,1)=(ty./(1/2+lambda_bottom/hy)).*Ly(:,1);
TY(1:Nx,Ny+1)=(ty./(1/2+lambda_top/hy)).*Ly(:,Ny);


x1=reshape(TX(1:Nx,:),N,1); x2=reshape(TX(2:Nx+1,:),N,1);
y1=reshape(TY(:,1:Ny),N,1); y2=reshape(TY(:,2:Ny+1),N,1);

x12=reshape(TX2(1:Nx,:),N,1); x22=reshape(TX2(2:Nx+1,:),N,1);
y12=reshape(TY2(:,1:Ny),N,1); y22=reshape(TY2(:,2:Ny+1),N,1);

DiagVecs=[-y22,-x22,x1+x2+y1+y2,-x12,-y12];
DiagIndx=[-Nx,-1,0,1,Nx];
A=spdiags(DiagVecs,DiagIndx,N,N);

DiagVecs=[hx*hy*ones(N,1)];
DiagIndx=[0];
C=spdiags(DiagVecs,DiagIndx,N,N);
A=A+C;
end
%RHS=xi*hx*hy;