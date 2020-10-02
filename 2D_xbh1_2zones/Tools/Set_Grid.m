function Grid=Set_Grid(n,L)
Grid.N=prod(n);
Grid.n=n;
Grid.D=L;
dim=length(n);
Grid.dim=dim;
for idim=1:Grid.dim
    Grid.h(idim)=L(idim)/n(idim);
end
[Grid.X,Grid.Y]=meshgrid( ...
            linspace(-L(1)/2+Grid.h(1)/2,L(1)/2-Grid.h(1)/2,n(1)),...'
            linspace(-L(2)/2+Grid.h(2)/2,L(2)/2-Grid.h(2)/2,n(2))  );
