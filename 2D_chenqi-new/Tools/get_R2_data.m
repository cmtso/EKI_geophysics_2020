
function data=get_R2_data(srv_file)

fid = fopen(srv_file,'r');
%nelec = str2double(fgetl(fid)) ;
%fscanf(fid,'%f',[5 nelec])';
%fgetl(fid);
%fgetl(fid);

nd = str2double(fgetl(fid)) ;
data = fscanf(fid,'%f%f%f%f%f%f %*[^\n]',[6 nd])'; %skip all after the 6th column
data = data(:,6);

if numel(data) < nd
  frewind(fid);
  nd = str2double(fgetl(fid)) ;
  data = fscanf(fid,'%f%f%f%f%f%f',[6 nd])'; % read strictly six columns
  data = data(:,6);
end

if numel(data) < nd
  
  fprintf('WARNING: not all data points read.');
end

fclose(fid);
end