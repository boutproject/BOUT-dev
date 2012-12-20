function info_file(filename)
% netCDF file viewer
% Input must be string format
% NOTE!! This function only shows variable list in netCDF file directed by user.
% Therefore, there is no return value
% Example)
%       info_file('bout.grd.nc')
%       info_file('data/BOUT.dmp.0.nc')
% 
% Coded by Minwoo Kim(Mar. 2012)

% Check input arguments
if ( nargin < 1 )
    fprintf('\tIt needs input arguments.\n');
    return
end
if ( ~ischar(filename) )    
    fprintf('\tInput must be string format.\n');    
    return
end

% Check netCDF file existence
try
    netcdf.open(filename, 'nc_nowrite');
catch exception
    if ( strcmp(exception.identifier, 'MATLAB:netcdf:open:noSuchFile') )
        fprintf('\tCan not open the file. Please, check file name and path.\n');
        return
    end
end

% Open netCDF file
fid = netcdf.open(filename, 'nc_nowrite');
[~, num_var, ~, ~ ] = netcdf.inq(fid);
% [num_dim, num_var, num_global_atts, unlim_dim_ID ] = netcdf.inq(fid);

fprintf( '\n' );
fprintf( 'Variable list of netCDF file\n\n' );
fprintf( '\t Var_ID\t\t Var_Name\t\tDimension\n' );

for i=1 : num_var

    % Inquire saved variable information
    [var_name, ~, ~, ~] = netcdf.inqVar(fid, i-1);
%     [var_name, xtype, dim_ids, num_atts] = netcdf.inqVar(fid, i-1);
    if ( length(var_name) < 7 )
        fprintf( '\t %d\t\t %s\t\t', i-1, var_name );        
    else
        fprintf( '\t %d\t\t %s\t', i-1, var_name );
    end
        
    var_tmp = netcdf.getVar( fid, i-1 );
    if ( length(size(var_tmp)) == 4 ) % 4 dimension
        var_tmp = permute(var_tmp, [3 2 1 4]); % [Z, Y, X, T] -> [X, Y, Z, T]
    elseif ( length(size(var_tmp)) == 3 ) % 3 dimension
        var_tmp = permute(var_tmp, [3 2 1]); % [Z, Y, X] -> [X, Y, Z]
    elseif ( length(size(var_tmp)) == 2 && sum(size(var_tmp)) > 2 ) % 2 dimension
        var_tmp = permute(var_tmp, [2 1]); % [Y, X] -> [X, Y]   
    end
    var_dim = size(var_tmp);
    fprintf( '\t%d', var_dim(1) );
    for j= 2:(length(var_dim))
        fprintf( ' X %d', var_dim(j) );
    end
    fprintf( '\n' );
    
end

netcdf.close(fid);

end