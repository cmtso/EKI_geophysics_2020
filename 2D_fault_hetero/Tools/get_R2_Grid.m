function R2_Grid=get_R2_Grid(filename,trn) 
  
% read R2 grid in the format of R2_forward.dat 
% (obtained from running forward model)
% perform translation given trn, a pair of x,y coordinates

 if nargin < 2
     trn = [0 0] ;
 end


% will translate to make mid point of (x,y) = (0,0)

ec=read_elements(filename);
R2_Grid.x  = ec(:,1) - trn(1);
R2_Grid.y  = ec(:,2) - trn(2);
R2_Grid.ec = [R2_Grid.x R2_Grid.y];
R2_Grid.trn = trn;



end
function ec=read_elements(filename)
    fid = fopen(filename,'r');
    
    %nd = str2double(fgetl(fid)) ;
    ec = fscanf(fid,'%f',[4 Inf])';
    ec = ec(:,1:2);

    fclose(fid); 


end
