
function vtk = read_vtk()

% Michael Tso 191010 
% read 2D finite element mesh
% tested on R2 ERT code by Andrew Binley (both rectangular and triangular mesh)

%clear scalar_list scalar_data

% listdlg: https://uk.mathworks.com/help/matlab/ref/listdlg.html#d117e828790

[file folder]=uigetfile('*.vtk','Select a *.vtk file -change file types to All Types *.*'); %select the file (res.dat file) (both the file name and folder will be saved)

  
%% grab scalars line numbers

fid = fopen(fullfile(folder,file)) ;
data = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
idx = strfind(data{1,1},'POINTS') ;
node_line = find(~cellfun(@isempty,idx));
idx = strfind(data{1,1},'CELLS') ;
cell_line = find(~cellfun(@isempty,idx));
idx = strfind(data{1,1},'SCALARS') ;
scalars_line = find(~cellfun(@isempty,idx));

%% read nodes and cells
fid = fopen(fullfile(folder,file)) ;
for i = 1:5
  str = fgetl(fid); 
end
nnode = sscanf(strtrim(str), 'POINTS %d double');
nodes = fscanf(fid,'%f',[3 nnode])';
fgetl(fid);
str = fgetl(fid);
ncell =  sscanf(strtrim(str), 'CELLS %d %d');
cells = fscanf(fid,'%f',[ncell(2)/ncell(1) ncell(1)])';
cells = cells(:,2:end) ; 
cells = cells + 1; % .vtk cells are zero-based
ncell = ncell(1);

%% read scalars
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
str = fgetl(fid);
fgetl(fid);
for i = 1:numel(scalars_line)
    str = data{1,1}{scalars_line(i)} ; 
    scalar_list{i,1} = sscanf(strtrim(str), 'SCALARS %s double 1');
    scalar_data(:,i) = str2num(data{1,1}{scalars_line(i)+2})' ; % str2double didn't work in octave
end
fclose(fid)

vtk.nodes = nodes;
vtk.cells = cells;
vtk.scalar_list = scalar_list;
vtk.scalar_data = scalar_data;
vtk.file = file;
vtk.folder = folder;

end