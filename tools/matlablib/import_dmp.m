function var = import_dmp(path, var_name)
% Import and collect simulation result from netcdf dump files.
% Convert IDL routine(collect) to Matlab
% Example)
%       p = inport_dmp('data', 'P');     
% Two inputs must be string format. If not, function return zero.
% If variable name does not match that in netCDF file, function print error message and return zero.
% NOTE!! This function is case-sensitive. Use 'info_file' to get correct variable names.
% 
% Coded by Minwoo Kim(Mar. 2012)

error("This is currently broken for BOUT++ > v4.0.0. See issue #394")

% Check input arguments
if ( nargin < 2 )
    fprintf('\tBoth dump file path and variable name are requisite input arguments.\n');
    fprintf('\tRetrun value is zero.\n');
    var = 0;
    return
end
if ( ~ischar(path) || ~ischar(var_name) )    
    fprintf('\tInput must be string format.\n\tRetrun value is zero.\n');
    var = 0;
    return
end    

% Check path of BOUT dum files
if ( ~exist(path, 'dir') )
    fprintf('\tThe directory (%s) does not exist.\n', path);
    fprintf('\tRetrun value is zero.\n');
    var = 0;
    return
end

% Check variable name in BOUT dump files
filename = [path '/BOUT.dmp.0.nc'];
fid = netcdf.open(filename, 'nc_nowrite');
try
    netcdf.inqVarID(fid,var_name);
catch exception
    if ( strcmp(exception.identifier, 'MATLAB:netcdf:inqVarID:variableNotFound') )
        fprintf('\tPlease, check the variable name and BOUT dump file path.\n');
        fprintf('\tRetrun value is zero.\n');
        var = 0;
        return
    end    
end

% Import basic parameter for collect certain variable in netcdf dump
% The # of processors in X/Y directions
vid = netcdf.inqVarID( fid, 'NXPE' );   NXPE = netcdf.getVar( fid, vid );
vid = netcdf.inqVarID( fid, 'NYPE' );   NYPE = netcdf.getVar( fid, vid );

% The # of X/Y grid points in each processor
vid = netcdf.inqVarID( fid, 'MXSUB' );  MXSUB = netcdf.getVar( fid, vid );
vid = netcdf.inqVarID( fid, 'MYSUB' );  MYSUB = netcdf.getVar( fid, vid );

% Sizes of the X/Y guard cells
vid = netcdf.inqVarID( fid, 'MXG' );    MXG = netcdf.getVar( fid, vid );
vid = netcdf.inqVarID( fid, 'MYG' );    MYG = netcdf.getVar( fid, vid );

% The # of Z points
% vid = netcdf.inqVarID( fid, 'MZ' );

% nx = NXPE*MXSUB + 2*MXG;
% ny = NYPE*MYSUB;
% nz = netcdf.getVar( fid, vid );

% Investigate dimension of variable to import
vid = netcdf.inqVarID( fid, var_name );
var_tmp = netcdf.getVar(fid, vid);
if ( length(size(var_tmp)) == 4 ) % 4 dimension
    option = 4;
elseif ( length(size(var_tmp)) == 3 ) % 3 dimension
    option = 3;
elseif ( length(size(var_tmp)) == 2 && sum(size(var_tmp)) > 2 ) % 2 dimension
    option = 2;
else % Constant number
    var = var_tmp;
    return
end

netcdf.close(fid);

%  Collect data
% [Z, Y, X, T] / [Z, Y, X] / [Y, X] / Constant number
for i=0 : (NXPE*NYPE-1)
    
    num_ps = num2str(i);
    filename = [path '/BOUT.dmp.' num_ps '.nc'];
    
    fid = netcdf.open(filename, 'nc_nowrite');
    
    vid = netcdf.inqVarID( fid, var_name );
    var_tmp = netcdf.getVar(fid, vid);

    idx_x = mod(i, NXPE);
    idx_y = floor(double(i)/double(NXPE));
    
    switch option
        case 4; % 4D data, [Z, Y, X, T]
            var_tmp = permute(var_tmp, [3 2 1 4]); % [Z, Y, X, T] -> [X, Y, Z, T]
            if ( idx_x == 0 )
                var( (idx_x*MXSUB+1):(idx_x*MXSUB+MXSUB+MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB), :, :)...
                    = var_tmp( 1:end-MXG, (MYG+1):end-MYG, :, :);                
            elseif ( idx_x == (NXPE-1) )
                var( (idx_x*MXSUB+1+MXG):(idx_x*MXSUB+MXSUB+2*MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB), :, :)...
                    = var_tmp( (MXG+1):end, (MYG+1):end-MYG, :, :);                
            else
                var( (idx_x*MXSUB+1+MXG):(idx_x*MXSUB+MXSUB+MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB), :, :)...
                    = var_tmp( (MXG+1):end-MXG, (MYG+1):end-MYG, :, :);                
            end
            
        case 3; % 3D data, [Z, Y, X]
            var_tmp = permute(var_tmp, [3 2 1]); % [Z, Y, X] -> [X, Y, Z]
            if ( idx_x == 0 )
                var( (idx_x*MXSUB+1):(idx_x*MXSUB+MXSUB+MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB), : )...
                    = var_tmp( 1:end-MXG, (MYG+1):end-MYG, : );                
            elseif ( idx_x == (NXPE-1) )
                var( (idx_x*MXSUB+1+MXG):(idx_x*MXSUB+MXSUB+2*MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB), : )...
                    = var_tmp( (MXG+1):end, (MYG+1):end-MYG, : );                
            else
                var( (idx_x*MXSUB+1+MXG):(idx_x*MXSUB+MXSUB+MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB), : )...
                    = var_tmp( (MXG+1):end-MXG, (MYG+1):end-MYG, : );                
            end
            
        case 2; % 2D data, [Y, X]
            var_tmp = permute(var_tmp, [2 1]); % [Y, X] -> [X, Y]
            if ( idx_x == 0 )
                var( (idx_x*MXSUB+1):(idx_x*MXSUB+MXSUB+MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB) )...
                    = var_tmp( 1:end-MXG, (MYG+1):end-MYG );                
            elseif ( idx_x == (NXPE-1) )
                var( (idx_x*MXSUB+1+MXG):(idx_x*MXSUB+MXSUB+2*MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB) )...
                    = var_tmp( (MXG+1):end, (MYG+1):end-MYG );               
            else
                var( (idx_x*MXSUB+1+MXG):(idx_x*MXSUB+MXSUB+MXG), (idx_y*MYSUB+1):(idx_y*MYSUB+MYSUB) )...
                    = var_tmp( (MXG+1):end-MXG, (MYG+1):end-MYG );                
            end
    end
    
    netcdf.close(fid);
    
end

% Convert single format to double format
var = double(var);

end







