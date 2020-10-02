function plot_resistivity
%% This matlab script plot the _res.dat results with the same unstructured grid as being used in the inversion
% last updated 180322 at Lancaster University
% -- faster IO, better aesthetics, 
% - read both rectangular and tri mesh
% - not saving, allow  interaction

[file folder]=uigetfile('','Select a *res.dat file -change file types to All Types *.*'); %select the file (res.dat file) (both the file name and folder will be saved)

%% open mesh-file file (geometry that will be used to plot the results on)

fid=fopen(fullfile(folder,'mesh.dat'),'r');

%% count # lines
nlines  = 0;
while fgets(fid)~= -1
    nlines = nlines+1;
end
frewind(fid);
%%
%num_triangles = str2double(fgetl(fid));
%%
%% read file-content linewise
data = cell(nlines,1);
for n = 1:nlines
    data{n} = fgetl(fid);
end

empty = []; %% remove empty lines
for i =1:size(data,1)
    wspace = isstrprop(data(i),'wspace');
    if sum (wspace{1})== numel(wspace{1}) %% whether the entire line is just white space
        empty = [empty i]; 
    end
end
data(empty) = [];

line1=sscanf(data{1},'%d %d');
num_triangles=line1(1);
num_points=line1(2);
%%write triangle file and point coordinate file out for further use
%%
lines=1:nlines;
formatSpec='%s\r\n';

fileID=fopen(fullfile(folder,'triangle_points.dat'),'w');
for n= (1:num_triangles)+1
    fprintf(fileID,formatSpec,data{n,:});
end
fclose(fid);

fileID=fopen(fullfile(folder,'point_coo.dat'),'w');
for n= (1:num_points)+1+num_triangles
    fprintf(fileID,formatSpec,data{n,:});
end
fclose(fid)


%
mesh=load(fullfile(folder,'triangle_points.dat'));
coordinates=load(fullfile(folder,'point_coo.dat'));

%% X and Y coordinates of the triangle corners

for i=1:length(mesh)
    x=[coordinates(mesh(i,2),2);coordinates(mesh(i,3),2);coordinates(mesh(i,4),2)]; 
    y=[coordinates(mesh(i,2),3);coordinates(mesh(i,3),3);coordinates(mesh(i,4),3)];
    triangle(i).coo=[x,y];
end
%%%%%%%%%%%%%%%%%%%%%%%


result=load(fullfile(folder,file));

%%%% a mask that leaves out the areas that are unsensitive/non-sensical
%%%% reduces the risk of making overconfident or wrong interpretations
%%%% One can run a second inversion, starting with a different starting
%%%% model and use the method as desribed by Oldenburg and Li (1999) or one
%%%% can use the diagonal of the [Jt Wt W J] matrix as an estimate of
%%%% senstivity over the profile (set RES matrix in line 2 over the R.in
%%%% file to '1'). There are no clear defined boundaries between what
%%%% values constitute sensitive and what values are insensitive areas.
%%%% However, higher values are more sensitive areas. An example of a
%%%% comparison between both methods can be found in the appendix of
%%%% Parsekian et al. (2017). The user is encouraged to assess the depth of
%%%% sensitivity and develop a mask suited for the field conditions

mask=load(fullfile(folder,'f001_sen.dat'));

%%%% The patch command draws triangles with the color of the triangles
%%%% scaled within the caxis values, if the values lie outside the range,
%%%% the color assinged to the triangle will be equal to the minimum or
%%%% maximum value in the caxis.
%%%% The patch command requires: x-coordinates of the corners,
%%%% y-coordinates of the corners, the value assigned (res-value) to the
%%%% triangle, CDataMapping is the command that will scale the values and
%%%% assign it a color to the caxis range on the colormap, EdgeColor is set
%%%% to none, otherwise you get a black edge on each triangle (Very messy
%%%% if you have +10000 triangles, and FaceAlpha is a method to set the
%%%% transparency of the triangle 1.0 is no transparency, 0.00 is fully
%%%% transparent. Set 'Edgealpha' = 0 to remove mesh boundaries.
fig1=figure(1);
clf
for i=1:length(triangle)
%   patch(triangle(i).coo(:,1),triangle(i).coo(:,2),result(i,4),'CDataMapping','scaled','EdgeColor','None','FaceAlpha',1.0,'Edgecolor',[0.75 0.75 0.75],'Edgealpha',0.3) %% 

    if mask(i,4)>-1
        patch(triangle(i).coo(:,1),triangle(i).coo(:,2),result(i,4),'CDataMapping','scaled','EdgeColor','None','FaceAlpha',1.0,'Edgecolor',[0.75 0.75 0.75],'Edgealpha',0.3) %% 
    elseif (mask(i,4)>-2 && mask(i,4)<=-1)
        patch(triangle(i).coo(:,1),triangle(i).coo(:,2),result(i,4),'CDataMapping','scaled','EdgeColor','None','FaceAlpha',0.4,'Edgecolor',[0.75 0.75 0.75],'Edgealpha',0.3)
    else
        patch(triangle(i).coo(:,1),triangle(i).coo(:,2),result(i,4),'CDataMapping','scaled','EdgeColor','None','FaceAlpha',0.2,'Edgecolor',[0.75 0.75 0.75],'Edgealpha',0.3)
    end
end

C = colorbar('so') ; xlabel('X [m]'), ylabel('Y [m]')
set(get(C,'YLabel'),'String','Inverted resistivity [\Omega m]','FontSize',12)
title('R2 ERT inversion','FontSize',12,'FontWeight','bold')
colormap('jet')
caxis([1.999 2.007]); %the range in which you expect the resistivity to be, or set to 'auto'
xl = [-12 60]; yl = [-20 0];
xlim(xl),ylim(yl) % min and max values for the x coordinates and y coordinates

pbaspect([ range(xl)/range(yl) 1 1])  % keep aspect ratio of plot
%%%% saving the figure in a user specified folder
figuredirectory=uigetdir(folder,'Choose the folder to save your figure');
%%%% user specified name
prompt = 'Provide a name for the figure (letters and numbers are accepted) or hit enter';
str = input(prompt,'s');
if isempty(str)
    str='Inverted Resistivity';
end

saveas(gcf,fullfile(figuredirectory,str),'png');

end


% for i=1:length(triangle)
%    patch(triangle(i).coo(:,1),triangle(i).coo(:,2),mask(i,4),'CDataMapping','scaled','EdgeColor','None','FaceAlpha',1.0,'Edgecolor',[0.75 0.75 0.75],'Edgealpha',0.3,'LineStyle','-') %% 
% end

