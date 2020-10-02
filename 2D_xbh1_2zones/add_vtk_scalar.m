function vtk = add_vtk_scalar(vtk, add_list, add_data)
  
  ### parse error????
  
if (any(size(add_data) == size(vtk.scalar_data,1))
##    if (size(add_data,2)  == size(vtk.scalar_data,1) )
##          add_data = add_data';
##    end
    vtk.scalar_data = [vtk.scalar_data add_data] ;
    vtk.scalar_list(end+1:end+numel(add_list)) = add_list; 
  
  
else
  fprintf("Error: Data dimension does not match those already in the vtk")
end

if ( numel(vtk.scalar_list) > numel(unique(vtk.scalar_list)))
  fprint("Warning: duplicate variable names")
end
          
end