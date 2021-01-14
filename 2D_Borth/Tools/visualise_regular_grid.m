function visualise_regular_grid(Grid,field)

field=reshape(field,Grid.n(1),Grid.n(2),Grid.n(3));
y_slice=[-2.5,2.5];%[yc(ceil(ny/2)),yc(ny)-dy(1)/2];
x_slice=[-3,3];%[xi(ceil(nx/2)),xi(nx)-dx(1)];
z_slice=[-3.0];%,-5.5];

%x_slice=[-Grid.D(1)/2+Grid.h(1)/2,Grid.D(1)/2-Grid.h(1)/2,Grid.h(1)/2];%[xi(ceil(nx/2)),xi(nx)-dx(1)];
%y_slice=[-Grid.D(2)/2+Grid.h(2)/2,Grid.D(2)/2-Grid.h(2)/2,Grid.h(2)/2];%[xi(ceil(nx/2)),xi(nx)-dx(1)];
%z_slice=-Grid.D(3)/2;%[-Grid.D(3)+Grid.h(3)/2,Grid.D(3)/4-Grid.h(3)/2];%,-Grid.h(3)/2];%[zc(1),zc(ceil(nz/2)),zc(nz)]
slice(Grid.X,Grid.Y,Grid.Z,field,x_slice,y_slice,z_slice,'nearest');
colormap jet

end