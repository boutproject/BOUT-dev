function plot_pol_slice(var3d, gfile, period, zangle)
% Project BOUT simulation results to poloidal cross section
% Convert IDL routine(plotpolslice) to Matlab
% There is no output from this function
% Example)
%       plot_pol_slice(var3d, gfile, period, zangle)
%       plot_pol_slice( P(:, :, :, time_step), 'bout.grd.nc', 5 )
%       var3d : 3-dim array, [X Y Z] at specific time step
%       gfile : string format, 'path/grid_file_name'
%       period : default = 1
%       zangle : default = 0
%
% Coded by M. Kim (Mar. 12)

% Check input arguments
if ( nargin < 2)
    fprintf('\tIt needs at least 2 input argumetns.\n');
    return
end
var3d = squeeze(var3d);
var3d_size = size(var3d);
if ( length(var3d_size) ~= 3 )
    fprintf('\tFirst arguments must be 3-dimension array.\n');
    return
end

if ( nargin < 3)
    period = 1;
    zangle = 0;
end
if ( nargin < 4)
    zangle = 0;
end
period = fix( abs(period) );
if ( period < 1 )
    period = 1;
end

% Check grid file path
try
    netcdf.open(gfile, 'nc_nowrite');
catch exception
    if ( strcmp(exception.identifier, 'MATLAB:netcdf:open:noSuchFile') )
        fprintf('\tCan not open the grid file. Please, check file name and path.\n');
        return
    end
end

% Check toroidal shift variable ID
fid = netcdf.open(gfile, 'nc_nowrite');
var_name = 'zShift'; % <------- Modify this line if it is needed
try    
    netcdf.inqVarID(fid, var_name);
catch exception
    if ( strcmp(exception.identifier, 'MATLAB:netcdf:inqVarID:variableNotFound') )
        fprintf('\tCan not find netCDF variable: %s\n', var_name);
        var_name = 'qinty'; %<-------- Modify this line if it is needed
    end    
end
if ( strcmp(var_name, 'qinty') )
    try
        netcdf.inqVarID(fid, var_name);
    catch exception
        if (  strcmp(exception.identifier, 'MATLAB:netcdf:inqVarID:variableNotFound') )
            fprintf('\tCan not find netCDF variable: %s\n', var_name);
            fprintf('\tDouble check the grid file or modify the (var_name) in this function\n');
        end
    end
end
netcdf.close(fid);

% Import variable from grid file
fid = netcdf.open(gfile, 'nc_nowrite');
vid = netcdf.inqVarID(fid, var_name);
zShift = netcdf.getVar(fid, vid);
zShift = double(zShift);
if ( length(size(zShift)) ~= 2 || sum(size(zShift)) == 2 )
    fprintf('\tPlease, check the grid file. zShift is not 2-dimension variable\n');
    return
end
zShift = permute(zShift, [2 1]);

vid = netcdf.inqVarID(fid, 'Rxy');
rxy = netcdf.getVar(fid, vid);
rxy = double(rxy);
rxy = permute(rxy, [2 1]);
vid = netcdf.inqVarID(fid, 'Zxy');
zxy = netcdf.getVar(fid, vid);
zxy = double(zxy);
zxy = permute(zxy, [2 1]);

netcdf.close(fid);

% Main part of this function
nx = var3d_size(1);
ny = var3d_size(2);
nz = var3d_size(3);

dz = 2*pi()/ period/ (nz-1);

nskip = max( abs( zShift(:, 2:end) - zShift(:, 1:end-1) ), [], 1 )/dz - 1;
nskip = round(nskip);
ny2 = ny + sum(nskip);
fprintf('\tNumber of poloidal points in output : %d\n', ny2);

var2d = zeros(nx, ny2);
rxy_n = zeros(nx, ny2);
zxy_n = zeros(nx, ny2);

y_n = 1;
for y = 1:(ny-1)
    % Original point
    for x = 1:nx
        zind = (zangle - zShift(x, y)) / dz + 1; % +1; Matlab array index start from '1'.
        var2d(x, y_n) = zinterp( var3d(x, y, :), zind, nz);
    end
    rxy_n(:, y_n) = rxy(:, y);
    zxy_n(:, y_n) = zxy(:, y);
        
    % Extra point
    for x = 1:nx

        zi0 = (zangle-zShift(x, y)) / dz + 1; % +1; Matlab array index start from '1'.
        zip = (zangle-zShift(x, y+1)) / dz + 1; % +1; Matlab array index start from '1'.
        
        dzi = (zip - zi0) / (nskip(y) + 1);
        range = 1:nskip(y);
        zi = zi0 + range * dzi;
        w = range/(nskip(y) + 1);
        
        var2d(x, y_n+range) = w.*zinterp(var3d(x, y+1, :), zi, nz) + (1-w).*zinterp(var3d(x, y, :), zi, nz);
        rxy_n(x, y_n+range) = w.*rxy(x, y+1) + (1-w).*rxy(x, y);
        zxy_n(x, y_n+range) = w.*zxy(x, y+1) + (1-w).*zxy(x, y);
        
    end
    y_n = y_n + nskip(y) + 1;
    fprintf('\t\ty : %d\t\t%d\n', y, y_n);
end

% Final point
for x= 1:nx
    zind = (zangle-zShift(x, ny)) / dz + 1; % +1; Matlab array index start from '1'.
    var2d(x, y_n) = zinterp(var3d(x, ny, :), zind, nz);    
end
rxy_n(:, end) = rxy(:, ny);
zxy_n(:, end) = zxy(:, ny);

figure;
surf(rxy_n, zxy_n, var2d); view(0,90); shading interp;
xlabel('R [m]', 'FontSize', 11);ylabel('z [m]', 'FontSize', 12);
cmax = max(max(abs(var2d)));
caxis([-cmax cmax]);

% Make color table
red = zeros(1,64);
green = red;
blue = red;

green(1:32) = 255 - 2*abs( (1:32)*4 - 128 );
green(33:64) = 255 - 2*abs( (32:-1:1)*4 - 128 );
blue(1:32) = 255;
blue(33:64) = 255 - 2*abs( (32:-1:1)*4 - 128 );
red(1:32) = 255 - 2*abs( (1:32)*4 - 128 );
red(33:64) = 255;

colormap_red2blue = [red; green; blue];
colormap_red2blue = colormap_red2blue'/255;
colormap(colormap_red2blue);


end


