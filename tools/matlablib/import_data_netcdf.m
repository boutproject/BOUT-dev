%% matlab file to Import and collect simulation data from netcdf dump files in BOUT++

function varn = import_data_netcdf(path, var_name,nt,ntsp)

% Import and collect simulation result from netcdf dump files.
% Modifications made to file provided by Minwoo Kim(Mar. 2012)
% By Sanat Kumar Tiwari    11Jan2014
% Institute for Plasma Research, India

% Example:
% varn = import_data_netcdf('.','P',20,4)
% loads variable P from current directory 
% varn = import_data_netcdf('.','P',20,1)

% path =  where the .nc data files located, in this example..present
% working directory
% var_name = The variable to be imported
% nt = number of time steps data imported (only for option 4)
% ntsp = time step at which the data to be imported out of total time array
% (only for option 4)
% Last two variables important only for [X,Y,Z,T] format and any number
% will be ok if we wish to plot [X,Y,Z] and [X,Y] type data.

error("This is currently broken for BOUT++ > v4.0.0. See issue #394")

% Check input arguments
if ( nargin < 4 )
    fprintf('\tBoth dump file path and variable name are requisite input arguments.\n');
    fprintf('\tRetrun value is zero.\n');
    varn = 0;
    return
end
if ( ~ischar(path) || ~ischar(var_name) )    
    fprintf('\tInput must be string format.\n\tRetrun value is zero.\n');
    varn = 0;
    return
end    

% Check path of BOUT dum files
if ( ~exist(path, 'dir') )
    fprintf('\tThe directory (%s) does not exist.\n', path);
    fprintf('\tRetrun value is zero.\n');
    varn = 0;
    return
end

% Check variable name in BOUT dump files
filename = [path '/BOUT.dmp.0.nc'];

% Check proper relation between nt and ntsp

% Size of data from each processor
pp = ncinfo(filename); [nnx nny nnz nnt] = pp.Dimensions.Length;
format = pp.Format;

if (nnt < nt*ntsp)
    fprintf('\t Total timesteps of saved data crossed \n');
    varn = 0;
    return
end


try
    ncinfo(filename,var_name);
catch exception
    if ( strcmp(exception.identifier, 'MATLAB:imagesci:netcdf:unknownLocation') )
        fprintf('\tPlease, check the variable name and BOUT dump file path.\n');
        fprintf('\tRetrun value is zero.\n');
        varn = 0;
        return
    end    
end  % End of try

%% This "if loop" is put there as it found that different netcdf versions of library 

if (format=='classic')
    % Import basic parameter for collect certain variable in netcdf dump
    % The # of processors in X/Y directions
    NXPE  = ncread(filename,'NXPE');   NYPE  = ncread(filename,'NYPE');
    % The # of X/Y grid points in each processor
    MXSUB = ncread(filename,'MXSUB');  MYSUB = ncread(filename,'MYSUB');
    % Sizes of the X/Y guard cells
    MXG   = ncread(filename,'MXG');    MYG   = ncread(filename,'MYG');
elseif (format=='netcdf4')
    [iteration MXSUB MYSUB MXG MYG MZ NXPE NYPE] = pp.Attributes(1:end).Value;
end

% Investigate dimension of variable to import
var_tmp = ncinfo(filename,var_name);
if ( length(var_tmp.Size) == 4 ) % 4 dimension
    option = 4;
elseif ( length(var_tmp.Size) == 3 ) % 3 dimension
    option = 3;
elseif ( length(var_tmp.Size) == 2 && sum(var_tmp.Size) > 2 ) % 2 dimension
    option = 2;
else % Constant number
    varn = var_tmp;
    return
end

START1 = [1 1 1 1]; COUNT1 = [nnz nny nnx nt]; STRIDE1= [1 1 1 ntsp];
%  Collect data
% [Z, Y, X, T] / [Z, Y, X] / [Y, X] / Constant number

for ii= 0 : (NXPE*NYPE-1)   % Loop over processor starts
    
    num_ps = num2str(ii);
    filename = [path '/BOUT.dmp.' num_ps '.nc'];
    
    if (option==4)
        var_tmp = ncread(filename, var_name,START1,COUNT1,STRIDE1);
    else
        var_tmp = ncread(filename, var_name);
    end

    idx_x = mod(ii, NXPE);
    idx_y = floor(double(ii)/double(NXPE));
 
% some constants for array operations
    nx1 = idx_x*MXSUB+1; nx2 = idx_x*MXSUB+MXSUB+MXG; 
    ny1 = idx_y*MYSUB+1; ny2 = idx_y*MYSUB+MYSUB;
    
    nx1e = idx_x*MXSUB+1+MXG; nx2e = idx_x*MXSUB+MXSUB+2*MXG;
    
        switch option
        case 4; % 4D data, [Z, Y, X, T]
            var_tmp = permute(var_tmp, [3 2 1 4]); % [Z, Y, X, T] -> [X, Y, Z, T]
            if ( idx_x == 0 )
                varn( nx1:nx2, ny1:ny2, :, :)...
                    = var_tmp( 1:end-MXG, (MYG+1):end-MYG, :, :);                
            elseif ( idx_x == (NXPE-1) )
                varn( nx1e:nx2e, ny1:ny2, :, :)...
                    = var_tmp( (MXG+1):end, (MYG+1):end-MYG, :, :);                
            else
                varn( nx1e:nx2, ny1:ny2, :, :)...
                    = var_tmp( (MXG+1):end-MXG, (MYG+1):end-MYG, :, :);                
            end
            
        case 3; % 3D data, [Z, Y, X]
            var_tmp = permute(var_tmp, [3 2 1]); % [Z, Y, X] -> [X, Y, Z]
            if ( idx_x == 0 )
                varn( nx1:nx2, ny1:ny2, :)...
                    = var_tmp( 1:end-MXG, (MYG+1):end-MYG, :);                
            elseif ( idx_x == (NXPE-1) )
                varn( nx1e:nx2e, ny1:ny2, :)...
                    = var_tmp( (MXG+1):end, (MYG+1):end-MYG, :);                
            else
                varn( nx1e:nx2, ny1:ny2, :)...
                    = var_tmp( (MXG+1):end-MXG, (MYG+1):end-MYG, :);                
            end
            
        case 2; % 2D data, [Y, X]
            var_tmp = permute(var_tmp, [2 1]); % [Y, X] -> [X, Y]
            if ( idx_x == 0 )
                varn( nx1:nx2, ny1:ny2)...
                    = var_tmp( 1:end-MXG, (MYG+1):end-MYG);                
            elseif ( idx_x == (NXPE-1) )
                varn( nx1e:nx2e, ny1:ny2)...
                    = var_tmp( (MXG+1):end, (MYG+1):end-MYG);                
            else
                varn( nx1e:nx2, ny1:ny2)...
                    = var_tmp( (MXG+1):end-MXG, (MYG+1):end-MYG);                
            end
            
        end % End for Switch 
    
end % Loop over processor ends

% Convert single format to double format
varn = double(varn);
end % Function end 
