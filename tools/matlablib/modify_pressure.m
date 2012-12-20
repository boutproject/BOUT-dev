% Modify pressure profile
% No radial gradient in pressure at divertor region

fid = netcdf.open('bout.grd_KSTARg006123.03031_nx516ny64_core086sol101.new_p.nc', 'nc_write');

inp_var_list = {'pressure', 'npol', 'nrad', 'psixy', 'psi_axis', 'psi_bndry'};

% [T, Z, Y, X]
% [Y, X]
% vid = netcdf.inqVarID( fid, 'pressure');
% pressure = netcdf.getVar(fid, vid);
% pressure = double(pressure);

for i=1 : length(inp_var_list)
    
    var_name = char(inp_var_list(i));
    vid = netcdf.inqVarID( fid, var_name);
    command = [var_name '=netcdf.getVar(fid, vid);'];
    eval(command);
    command = [var_name '=double(' var_name ');'];
    eval(command);
    
end

psi_norm = (psixy - psi_axis)/(psi_bndry - psi_axis);

p0 = 3100; % top value
ps = 50; % steepness
px0 = 0.96; % pedestal center
pmin = 100; % minimum value

p_new = pmin + p0 * ( 1 - tanh(ps * (psi_norm - px0)) );

p_new(1:npol(1), 1:nrad(1)) = p_new(1, nrad(1));
p_new(end-npol(1)+1:end, 1:nrad(1)) = p_new(1, nrad(1));

% surf(p_new);shading interp;
plot(psi_norm(40,:),p_new(40,:),'LineWidth',2);
% plot(psi_norm(40,:),p_new(40,:),'LineWidth',2, 'Marker', '^', 'MarkerSize', 2);
xlabel('\psi_{norm}','FontSize',14,'FontWeight','Bold');
ylabel('Pressure (Pa)','FontSize',14,'FontWeight','Bold');
title('Pressure profile at mid-plane','FontSize',16,'FontWeight','Bold');
xlim([0.85 1.05]);
set(gcf, 'color', [1 1 1]);
set(gca, 'FontSize', 12, 'FontWeight', 'Bold', 'Box', 'on');

figure;
plot(psi_norm(40,:),p_new(1,:),'LineWidth',2);
xlabel('\psi _{norm}','FontSize',14,'FontWeight','Bold');
ylabel('Pressure (Pa)','FontSize',14,'FontWeight','Bold');
title('Pressure profile at private flux region','FontSize',16,'FontWeight','Bold');
xlim([0.85 1.05]);
set(gcf, 'color', [1 1 1]);
set(gca, 'FontSize', 12, 'FontWeight', 'Bold', 'Box', 'on');

figure;
hold all;
plot(psi_norm(40,:),p_new(40,:),'LineWidth',2);
plot(psi_norm(40,:),p_new(1,:),'LineWidth',2);
legend('Mid-plane', 'PF');
xlabel('\psi _{norm}','FontSize',14,'FontWeight','Bold');
ylabel('Pressure (Pa)','FontSize',14,'FontWeight','Bold');
title('Pressure profile at private flux region','FontSize',16,'FontWeight','Bold');
xlim([0.85 1.05]);
set(gcf, 'color', [1 1 1]);
set(gca, 'FontSize', 12, 'FontWeight', 'Bold', 'Box', 'on');
hold off;

% plot(pressure(40,:));
% figure; plot(pressure(1,:));

% if you comment next line, the new pressure profile is overwritten 
return 
vid = netcdf.inqVarID( fid, 'pressure');
netcdf.putVar(fid, vid, p_new);

netcdf.close(fid);









