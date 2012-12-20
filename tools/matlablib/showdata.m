function showdata(var, scale, delay)
% Visualize BOUT++ result
% Convert IDL routine to Matlab
% Final index is assumed to be time
% Example)
%       showdata(p_rms, 'on', 0.2)
%       Only first argument is requisite
% Coded by Minwoo Kim(Mar. 2012)

delay_def = 0.3;

% Check input arguments
if ( nargin < 1 )
    fprintf('\tIt needs input arguments.\n');    
    return
end

if ( ischar(var) )
    fprintf('\tFirst input must be numeric format.\n');    
end
var = squeeze(var);
var_size = size(var);
if ( length(var_size) == 4 ) % 4 dimension
    fprintf('4-D sturcture visulization is not implemented yet.');
    return
elseif ( length(var_size) == 3 ) % 3 dimension
    option = 3;
elseif ( length(var_size) == 2 && sum(var_size) > 2 ) % 2 dimension
    option = 2;
else % Constant number
    fprintf('Input variable is just constant value.');
    return
end

if ( nargin < 2 )
    scale = 'off';
    delay =delay_def;
end
if ( nargin == 2 )
    if ( ischar(scale) )
        delay = delay_def;
    else
        delay = scale;
        scale = 'off';
    end
end
if ( ischar(delay) )
    delay = delay_def;
end

if ( strcmp(scale, 'on') )
    switch option
        case 2;
            axis_limit_max = max(max(var));
            axis_limit_min = min(min(var));
        case 3;
            axis_limit_max = max(max(max(var)));
            axis_limit_min = min(min(min(var)));
    end
    option = -option;
end

for i = 1:var_size(end)
   
    switch option
        
        case 2;
            plot( var( :, i ), 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 10 );
            set(gca, 'FontSize', 12);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 15);
            pause(delay);
         
        case -2;
            plot( var( :, i ), 'LineWidth', 1.5, 'Marker', 'x', 'MarkerSize', 10 );
            set(gca, 'FontSize', 12);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 15);
            ylim( [axis_limit_min axis_limit_max] );
            pause(delay);
            
        case 3;
            surf( var( :, :, i) );
            set(gca, 'FontSize', 12);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 15);
            pause(delay);
            
        case -3;
            surf( var( :, :, i) );
            set(gca, 'FontSize', 12);
            text_time = ['t = ' num2str(i)];
            xlabel(text_time, 'FontSize', 15);
            zlim( [axis_limit_min axis_limit_max] );
            pause(delay);
    
    end
end

end





