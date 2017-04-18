%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% plot linear growth rate and real frequency vs. ktheta*rhoi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ktheta_rhoi
qval=1.4;
rval=0.5;
for n1=nmin:ndel:nmax
    ktheta_rhoi(n1/ndel)=n1*qval/rval*rhos0;
end
R=1.3; a=0.48;  % for cyclone base case
gamfactor=R_Lne*a/R;  % change time unit from a/v_ti to L_ne/v_ti
omgfactor=gamfactor*4;

figure(2)
plot(ktheta_rhoi,gamlin/gamfactor,'-or');
xlabel('k_\theta \rho_i')
ylabel('\gamma  and  \omega/4  (v_{ti}/L_{ne})')
title('Linear growth rate and real frequency')
hold on
plot(ktheta_rhoi,omglin/omgfactor,'-ob');
