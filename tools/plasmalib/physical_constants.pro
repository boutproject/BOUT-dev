; physical_constants.pro
;
; Author:       Ilon Joseph
; Began:        2011/08/22
; Modified:     2011/09/02 
;
; defines a structure to hold important constants of physics 
; sets default plasma parameters: aeff, zeff, meff, loglam 


pro physical_constants

 common constants, c

 amu= 1.66053892e-24; %g
 Me = 9.10938188e-28; % g
 Mp = 1.67262158e-24; % g
 Md = 2.013553212724*amu;
 Mt = 3.0160492*amu;
 Mc = 12.0107*amu;
 Meff = Md;

 zeff=1.0; %Deuterium
 aeff=2.0; %Deuterium
 loglam=14.0; %Default value for Coulomb logarithm

 kb = 1.60217656535e-12; % erg/eV
 ee = 4.8032042510e-10;  % statcoulomb

 ggrav   = 6.6738480e-8;      % erg cm/g^2
 hplanck = 4.135667516e-15; % eV s
 clight  = 2.99792458e10;    % cm/s

 m = {e:me, p:mp, d:md, c:mc, i:mi}
 c = {clight:clight, hplanck:hplanck, ggrav:ggrav, $
      ee:ee, kb:kb, amu:amu, m:m, $
      aeff:aeff, zeff:zeff,  meff:meff, $
      loglam:loglam} 
            
end
