function visualise(n,ec,sigma,L)
cor=cell(1,3);
cor{1,1}=linspace(-L(1)/2,L(1)/2,n(1));
cor{1,2}=linspace(-L(2)/2,L(2)/2,n(2));
cor{1,3}=linspace(-L(3),0,n(3));
%for i=1:3
%    cor{1,i}=linspace(-L(i)/2,L(i)/2,n(i)*10);
%end

[X,Y,Z]=meshgrid(cor{1,1},cor{1,2},cor{1,3});
y_slice=[0];%[yc(ceil(ny/2)),yc(ny)-dy(1)/2];
x_slice=[0];%[xi(ceil(nx/2)),xi(nx)-dx(1)];
z_slice=[-5.0];%,-5.5];
field= griddata(ec(:,1),ec(:,2),ec(:,3),sigma,X,Y,Z,'natural');%%interpolate to e4d grid
figure
slice(X,Y,Z,field,x_slice,y_slice,z_slice);%,'cubic');
colormap jet

end