function growth_rate(var)
% Estimate growth rate
% Example)
%       growth_rate(rms(x, y, :))
% --> RMS time trace is plotted
% --> Click two points on plotted graph for selecting fitting range
% --> Funtion takes a linear fitting and print a growth rate
%
% Coded by Minwoo Kim(Mar. 2012)

% Check input arguments
if ( nargin < 1 )
    fprintf('\tIt needs input arguments.\n');
    return
end
if ( ischar(var) )
    fprintf('\tInput argument must be a numeric array.\n');
    return
end
var = squeeze(var);
if ( ~isvector(var) )
    fprintf('\tInput argument must contain time trace of rms value at specific postion. ');
    fprintf('Or, it must be 1-D array\n');
    return
end

timebase = 0 : (length(var)-1);

figure;
plot( var, 'LineWidth', 1.5 );
set(gca, 'Yscale', 'log', 'FontSize', 12);
xlabel('\tau_A', 'FontSize', 15);
ylabel('<P>_{rms}', 'FontSize', 15);
xlim([timebase(1) timebase(end)]);
line_yposition=get(gca, 'Ytick');

time = ginput(2);
ti = find(timebase >= time(1), 1, 'first');
tf = find(timebase <= time(2), 1, 'last');
line( [timebase(ti) timebase(ti)], [line_yposition(1) line_yposition(end)], 'color', 'red' );
line( [timebase(tf) timebase(tf)], [line_yposition(1) line_yposition(end)], 'color', 'red' );

x = timebase(ti:tf);
y = log(var(ti:tf));

if ( ~(size(x) == size(y)) )
    y = y';
end
poly_coef = polyfit(x, y, 1);

figure;
hold on;
plot( timebase, log(var), 'LineWidth', 1.5 );
set(gca, 'FontSize', 12, 'Box', 'on');
xlabel('\tau_A', 'FontSize', 15);
ylabel('ln(<P>_{rms})', 'FontSize', 15);
plot( timebase, poly_coef(1)*timebase+poly_coef(2), '--r', 'LineWidth', 1.5 );
xlim([timebase(1) timebase(end)]);
hold off;

fprintf('\tGrowth rate (slope) : %f (normalized by Alfven frequency)\n', poly_coef(1));

end




