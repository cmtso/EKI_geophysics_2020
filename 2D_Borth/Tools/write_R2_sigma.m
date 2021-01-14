function write_sigma(file,sigma)
    if size(sigma,2) == 1 % making sure sigma is a column vector
        sigma = sigma';
    elseif min(size(sigma)) ~= 1
        fprintf('Quitting! Sigma is not a vector')
        return
    end
    
    % change R2.in (should just change it in template)
    
    %write resistivity.dat (in the format of  _res.dat)
    
    fid = fopen(file,'w');
    fprintf(fid,'\t%8.6e\t%8.6e\t%8.6e\t%8.6e \n', ...
                [zeros(numel(sigma),2)  1./ sigma' log10(1./ sigma')]') ;
    fclose(fid);
end