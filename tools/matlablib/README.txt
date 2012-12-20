Simple explanation of Matlab library


var = import_dmp(path, var_name)
	same as IDL routine, collect
	path: string format, 'path_of_BOUT_dmp_file'
	var_name: string format, 'variable_name', case-sensitive


info_file('path/name')
	netCDF file viewer

res = moment(variable, option)
	same as IDL routine, moment_xyzt
	option: string format, 'rms', 'ac' and 'dc'
	return the calculation result

showdata(variable, option, delay)
	same as IDL routine, showdata
	option: string format, for fixed axis scale option; 'on', 'off', default = 'off'
	delay: delay between each frame, default = 0.3

growth_rate(variable)
	variable: 1-dim array
	estimate growth rate of instability

plot_pol_slice(var3d, gfile, period, zangle)
	similar to IDL routine, plot_pol_slice
	var3d: 3-dim array
	gfile: string format, 'path_of_grid_file/name_of_grid_file'
	period: default = 1
	zangle: default = 0