function result = moment(var, option)
% Obtain moment (rms, dc and ac)
% Convert IDL routine(moment_xyzt) to Matlab
% Example)
%       p_rms = moment(p, 'rms');
% option : it must be string format; 'rms', 'ac' and 'dc' 
% 
% Coded by Minwoo Kim(Mar. 2012)

% Check input arguments
if ( nargin < 1 )
    fprintf('\tIt needs input arguments.\n');
    fprintf('\tReturn value is zero.');
    result = 0;
    return
end
if ( ischar(var) )
    fprintf('\tFirst input must be numeric format.\n');
    fprintf('\tReturn value is zero.');
    result = 0;
    return
end
if ( length(size(var))  ~= 4 )
    fprintf('\tInput variable must be 4-D array.\n');
    fprintf('\tReturn value is zero.');
    result  = 0;
    return
end
if ( nargin < 2 )
    option = 'rms';
    fprintf('\tThe code set the option as "rms".\n');
end

var_size = size(var);

nx = var_size(1);
ny = var_size(2);
nz = var_size(3);
nt = var_size(4);

avg = zeros(nx, ny, nt);
for i = 1:nz
    avg = avg + squeeze(var( :, :, i, :));
end
avg = avg / nz;

switch option
    case {'dc', 'DC'}
        result = avg;
        
    case {'rms', 'RMS'}
        rms = zeros(nx, ny, nt);
        for i = 1:nz
%             rms = rms + squeeze(var( :, :, i, :)) .* squeeze(var( :, :, i, :));
            rms = rms + ( squeeze(var( :, :, i, :)) - avg ).^2;
        end
        rms = sqrt(rms/nz);
        result = rms;
        
    case {'ac', 'AC'}
        ac = zeros(nx, ny, 1, nt);
        tmp = zeros(nx, ny, 1, nt);
        for i = 1:nz
            tmp = tmp + var( :, :, i, :);
        end
        tmp = tmp/ nz; % average for ac
        for i = 1:nz
            ac( :, :, i, :) = var( :, :, i, :) - tmp;
        end
        result = ac;
    otherwise
        fprintf('\tChoose proper option(rms, ac, dc)\n');
        fprintf('\tReturn value is zero');
        result = 0;        
end

end


