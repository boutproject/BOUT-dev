%======================================================================
% Preprocess for eigenvalue equation solver
% Calculate radial basis functions, arrange mode index, 
% and evaluate matrix elements involving integrals.
% This part was developed based on 
% TRB code [X. Garbet et.al., Phys. Plasmas 8, 2793 (2001)].
%======================================================================

      jmax=256; % number of radial grids
      nmax=15;  % maximum of toroidal mode number
      nmin=5;   % minimum of toroidal mode number
      ndel=5;   % toroidal mode number spacing
      mmax=150; % maximum of poloidal mode number
      mmin=1;   % minimum of poloidal mode number
      mdel=1;   % poloidal mode number spacing
      kyrhoscut=1.3;  % maximal ktheta*rhoi allowed in the calculation
      xmu1=0.1; % viscous damping coefficient
      xmu2=0.0; % hyper-viscous damping coefficient
      r0=0.2;   % inner boundary
      r1=0.9;   % outer boundary

      hermite=1;
      bessel=2;
      basis=bessel;     % choose radial basis function
      if(basis==hermite)
        nwaverad=11;    % number of radial basis functions
        delwavescale=5; % specify width of Hermite function
      else
        nwaverad=50;
      end

status='******  Start initialization  ******'
%-----------------------------------------------------------------------
% Normalized radius and delta_r for integral 
%-----------------------------------------------------------------------
      dr=1/jmax;
      for j=1:jmax-1
         fac(j)=1.0;
         r(j)=j*dr;
      end
      fac(jmax)=0.5;
      r(jmax)=1;
      for j=2:jmax
         drm(j)=r(j)-r(j-1);
      end
      drm(1)=drm(2);
      for j=1:jmax
         rdr(j)=fac(j)*drm(j)*r(j);
      end
      ar=r;

%-----------------------------------------------------------------------
% Set equilibrium profiles : cyclone case used
%-----------------------------------------------------------------------
status='   Setting equilibrium profiles'

      R_Lne=1./0.45;      % R/L_ne, R=major radius, L_ne=1/(dlog(n)/dr) at r=0.5a 
      R_Lte=6.92;         % R/L_Te, L_Te=1/(dlog(Te)/dr) at r=0.5a
      R_Lti=6.92;         % R/L_Ti, L_Ti=1/(dlog(Ti)/dr) at r=0.5a
      Te_T0=1.0;          % Te/Ti ratio
      rmajor=1.3/0.48;    % R/a, a=minor radius
      rhos0=7.086026e-03; % rho_star=rho_i/a
      q0=0.854;           % q(r=0)
      q1=2.184;           % q(r=1)
      for j=1:jmax
        xx=(r(j)-0.5)/0.3;
        enhat(j)=exp(-0.3*R_Lne/rmajor*tanh(xx)); % density profile
        tehat(j)=Te_T0*exp(-0.3*R_Lte/rmajor*tanh(xx)); % Te profile
        that(j)=exp(-0.3*R_Lti/rmajor*tanh(xx)); % Ti profile
        pihat(j)=enhat(j)*that(j); % ion pressure profile
        tau(j)=that(j)/tehat(j);   % tau=Ti/Te
        rln(j)=-R_Lne/rmajor/cosh(xx)^2;  % a*dlog(n)/dr
        rlt(j)=-R_Lti/rmajor/cosh(xx)^2;  % a*dlog(Ti)/dr
        rlte(j)=-R_Lte/rmajor/cosh(xx)^2; % a*dlog(Te)/dr
        rltau(j)=rlt(j)-rlte(j);          % a*dlog(tau)/dr
        dpihat(j)=fac(j)*drm(j)*pihat(j)*(rln(j)+rlt(j)); % dPi/dr * dr
        dthat(j)=fac(j)*drm(j)*that(j)*rlt(j);            % dTi/dr * dr
        denhat(j)=fac(j)*drm(j)*enhat(j)*rln(j);          % dn/dr * dr
        q(j)=q0+q1*r(j)^2; % q-profile
      end
      

%----------------------------------------------------------------------
%  Find resonant surfaces and determine number of basis functions
%-----------------------------------------------------------------------
status='   Finding rational surfaces and determining number of basis functions'

      p0=q0+q1*r0^2;
      p1=q0+q1*r1^2;     

      k=0;
      for n1=nmin:ndel:nmax % toroidal mode number scan
        m1min=p0*n1-mod(p0*n1,1)+1;
        m1min=max(m1min,mmin);
        m1max=p1*n1-mod(p1*n1,1)-1;
        m1max=min(m1max,mmax);

        for m1=m1min:mdel:m1max % poloidal mode number scan
          qtest=m1/n1;
          jqtest=0;
          for j=1:jmax-1
            if(((q(j)-qtest)*(q(j+1)-qtest)) < 0)
               jqtest=j;
            elseif(q(j)==qtest) 
               jqtest=j;
            end
          end

          pkyrhotest=m1/r(jqtest)*rhos0;
         
          if(pkyrhotest<=kyrhoscut)
            for iwaverad=1:nwaverad % radial mode number scan
              k=k+1;
            end
          end
        end
      end
      kmax=k; % total number of basis functions
      
%----------------------------------------------------------------------
%  Assign mode numbers 
%-----------------------------------------------------------------------
status='   Arranging mode index'

      nfb=0*ones(1,kmax);  % toroidal mode number n
      mf=0*ones(1,kmax);   % poloidal mode number m
      jss=0*ones(1,kmax);  % index of rational surface radius
      rss=0*ones(1,kmax);  % radius of rational surface
      qss=0*ones(1,kmax);  % q values on rational surfaces
      kw=0*ones(kmax,nwaverad);  % index of mode radially coupled with k-mode 
      kpw=0*ones(kmax,nwaverad); % index of mode poloidally coupled with k-mode
      kmw=0*ones(kmax,nwaverad); % index of mode poloidally coupled with k-mode    

      k=0;
      for n1=nmin:ndel:nmax
        m1min=p0*n1-mod(p0*n1,1)+1;
        m1min=max(m1min,mmin);
        m1max=p1*n1-mod(p1*n1,1)-1;
        m1max=min(m1max,mmax);

        for m1=m1min:mdel:m1max
          qtest=m1/n1;
          jqtest=0;
          for j=1:jmax-1
            if(((q(j)-qtest)*(q(j+1)-qtest)) < 0)
               jqtest=j;
            elseif(q(j)==qtest) 
               jqtest=j;
            end
          end

          pkyrhotest=m1/r(jqtest)*rhos0;
         
          if(pkyrhotest<=kyrhoscut)
            for iwaverad=1:nwaverad  % iwaverad = p is raidal mode number
              k=k+1;      % k = (n,m,p) 
              nfb(k)=n1;  % toroidal mode number
              mf(k)=m1;   % poloidal mode number
              jss(k)=jqtest;
              rss(k)=r(jqtest);
              qss(k)=qtest;

              for iwaveradp=1:nwaverad
                 kw(k,iwaveradp)=k-iwaverad+iwaveradp;

                 if(m1<m1max) 
                   kpaux=k-iwaverad+nwaverad+iwaveradp;
                   if(kpaux>0) 
                     kpw(k,iwaveradp)=kpaux; % index of m+1 mode
                   else
                     kpw(k,iwaveradp)=0;
                   end
                 else
                   kpw(k,iwaveradp)=0;
                 end

                 if (m1>m1min) 
                   kmaux=k-iwaverad-nwaverad+iwaveradp;
                   if(kmaux>0)
                     kmw(k,iwaveradp)=kmaux; % index of m-1 mode
                   else
                     kmw(k,iwaveradp)=0;
                   end
                 else
                   kmw(k,iwaveradp)=0;
                 end
              end
            end
          end
        end
      end

      for k=1:kmax
         n1=nfb(k);
         for iwaverad=1:nwaverad
           if(kpw(k,iwaverad)~=0)
            if(nfb(kpw(k,iwaverad))~=n1)
               kpw(k,iwaverad)=0;
            end
           end
           if(kmw(k,iwaverad)~=0)
            if(nfb(kmw(k,iwaverad))~=n1)
               kmw(k,iwaverad)=0;
            end
           end
         end
      end
      
%-----------------------------------------------------------------------
%  Plot equilibrium profiles, q profile, and q_mn=m/n for each n-mode 
%-----------------------------------------------------------------------
      figure(1)
      numfig=3+nmax/ndel;
      mfig=numfig/3;
      nfig=numfig/mfig;
      subplot(mfig,nfig,1)
      plot(r,enhat);
      hold on
      plot(r,that,'r');
      hold off
      xlabel('r/a')
      title('T_i:red, n_e:blue')
      subplot(mfig,nfig,2)
      plot(r,-rln);
      hold on
      plot(r,-rlt,'r');
      hold off
      xlabel('r/a')
      title('a*dlnT_i/dr:red, a*dlnn_e/dr:blue')
      subplot(mfig,nfig,3)
      plot(r,rlt./rln);
      xlabel('r/a')
      title('\eta_i = dlnT_i/dlnn_e')
      
      for n1=nmin:ndel:nmax
        subplot(mfig,nfig,3+n1/ndel)
        plot(r,q);
        hold on
        str1=['n=' int2str(n1)];
        text(0.2,q0+0.8*q1,str1);
        xlabel('r/a')
        title('q_{mn}=m/n and q profile') 
        k1=1;
        clear rss1 qss1
        for k=1:kmax
          if(nfb(k)==n1)
            rss1(k1)=rss(k);
            qss1(k1)=qss(k);
            k1=k1+1;
          end
        end
        plot(rss1,qss1,'o','MarkerFaceColor','r','MarkerSize',5);
        hold off
      end
      
%-----------------------------------------------------------------------
%  Radial basis functions 
%-----------------------------------------------------------------------
status='   Calculating radial basis functions:'

      wave=0*ones(jmax,kmax);    % normalized radial basis function
      dwavedr=0*ones(jmax,kmax); % radial derivative of radial basis function

   if(basis==hermite)
      delwave=0*ones(1,kmax);
      %  determine radial width of modes
      for k=1:kmax  
        delwave(k)=rhos0*delwavescale;
        if (rss(k)<(5/7*nwaverad*delwave(k))) 
          delwave(k)=7/5*rss(k)/nwaverad;
        elseif ((1-rss(k))<(4/7*nwaverad*delwave(k))) 
          delwave(k)=7/4*(1.-rss(k))/nwaverad; 
        end
      end
      %  Hermite functions
      status='   ... calculating Hermite functions and their derivatives'
      for k=1:nwaverad:kmax
        delwavek=delwave(k);
        for j=1:jmax
          xaux=(r(j)-rss(k))/delwavek;
          atmp=-(xaux/delwavek+0.5/r(j));

          wave(j,k)=exp(-0.5*xaux^2)/sqrt(r(j));
          wave(j,k+1)=2*xaux*wave(j,k);
          dwavedr(j,k)=atmp*wave(j,k);
    
          for iwaverad=3:nwaverad
            kp=k+iwaverad-1;
            wave(j,kp)=2*xaux*wave(j,kp-1)-2*(iwaverad-2)*wave(j,kp-2);
          end
    
          for iwaverad=2:nwaverad
            kp=k+iwaverad-1;
            dwavedr(j,kp)=atmp*wave(j,kp)+2*(iwaverad-1)/delwavek*wave(j,kp-1);
          end
        end
      end

   else % basis= bessel function 
      %  find zeros of bessel 
      status='   ... finding zeros of Bessel functions'
      alphamp=0*ones(mmax+1,nwaverad);
      zeros_bessel=0*ones(mmax+1,nwaverad);
      for m=1:mmax+1
          if(m==1)
              x0=2;
          else
              x0=alphamp(m-1,1);
          end
        ind=1;
        alphamp(m,1)=fzero(@(x)besselj(m-1,x),x0);
        while (ind < nwaverad)
          y0=besselj(m-1,x0);
          y1=besselj(m-1,x0+0.1);
          if(y0*y1<0) 
            y=fzero(@(x)besselj(m-1,x),[x0 x0+0.1]);
            z=abs(y-alphamp(m,ind));
            if (z>1.0e-13)
              ind=ind+1;  
              alphamp(m,ind)=y;
            end
          end
          x0=x0+0.1;
        end
        m_alphamp=m;
        %here m_alphamp
      end
      %  check zeros of bessel
      for m=1:mmax+1
          for ind=1:nwaverad-1
           if (alphamp(m,ind)>alphamp(m,ind+1))
              %here m
              %here ind
           end
          end
          for ind=1:nwaverad
             zeros_bessel(m,ind)=besselj(m-1,alphamp(m,ind));
          end
      end
      %  Bessel functions
      status='   ... calculating Bessel functions and their derivatives'
      kx=0*ones(1,kmax);
      for k0=1:nwaverad:kmax
          for iwaverad=1:nwaverad
             k=k0+iwaverad-1;
             ktmp_wave=k;
             %here ktmp_wave
             alphaj=alphamp(mf(k)+1,iwaverad);
             kperp(k)=rhos0*alphaj;
             for j=1:jmax
                x=alphaj*r(j);
                wave(j,k)=besselj(mf(k),x);
                dwavedr(j,k)=0.5*alphaj*(besselj(mf(k)-1,x)-besselj(mf(k)+1,x));
             end
          end
      end
   end

   % normalize basis functions
   for k=1:kmax
      anorm=0;
      for j=1:jmax
         anorm= anorm + abs(wave(j,k))^2*rdr(j);
      end
      anorm=sqrt(anorm);
      for j=1:jmax
         wave(j,k)=wave(j,k)/anorm;
         dwavedr(j,k)=dwavedr(j,k)/anorm;
      end
   end

%----------------------------------------------------------------------
%   Parallel wavenumber & Landau damping
%----------------------------------------------------------------------
status='   Calculating parallel wavenumber and Landau damping operator'

      kpar=0*ones(jmax,kmax);      % parallel wavenumber normalized by 1/a 
      akpar=0*ones(kmax,nwaverad); % for Grad_par
      akpar_enhat=0*ones(kmax,nwaverad); % for n * Grad_par
      akpar_tau=0*ones(kmax,nwaverad);   % for Ti/Te * Grad_par
      akpar_that=0*ones(kmax,nwaverad);  % for Ti * Grad_par
      kparglf_sqrtthat=0*ones(kmax,nwaverad); % for Landau damping

      % parallel wavenumber
      for k=1:kmax
         for j=1:jmax
           kpar(j,k)=(mf(k)/q(j)-nfb(k))/rmajor;
         end
      end

      % integrals including parallel wavenumber
      for k=1:kmax
       ktmp_akpar=k;
       %here ktmp_akpar
       for iwaverad=1:nwaverad
        k1=kw(k,iwaverad);
        for j=1:jmax
         atmp=wave(j,k)*wave(j,k1)*rdr(j);
         btmp=kpar(j,k)*atmp;
         akpar(k,iwaverad) = akpar(k,iwaverad) + btmp;
         akpar_enhat(k,iwaverad) = akpar_enhat(k,iwaverad) + enhat(j)*btmp;
         akpar_tau(k,iwaverad) = akpar_tau(k,iwaverad) + tau(j)*btmp;
         akpar_that(k,iwaverad) = akpar_that(k,iwaverad) + that(j)*btmp;
         kparglf_sqrtthat(k,iwaverad) = kparglf_sqrtthat(k,iwaverad) + abs(kpar(j,k))*sqrt(that(j))*atmp;
        end
        kparglf_sqrtthat(k,iwaverad)=sqrt(8./pi)*kparglf_sqrtthat(k,iwaverad);
       end
      end

      clear kpar

%----------------------------------------------------------------------
%    Viscous damping 
%----------------------------------------------------------------------
status='   Calculating viscous damping operator'

   if(basis==hermite)
      kx2=0*ones(1,kmax);
      kx2m=0*ones(1,kmax);
      kx2p=0*ones(1,kmax);
      kkx2m=0*ones(1,kmax);
      kkx2p=0*ones(1,kmax);
      ky=0*ones(1,kmax);
      d0kx2p=0*ones(1,kmax);
      d0kx2m=0*ones(1,kmax);
      xmu=0*ones(1,kmax);
      d0kx2=0*ones(1,kmax);
      
      for k0=1:nwaverad:kmax-nwaverad+1
       for iwaverad=1:nwaverad
         nrad=iwaverad-1;
         k=k0+nrad;
         kx2(k)=(nrad+0.5)*(rhos0/delwave(k))^2;
         if(iwaverad>=3)
            kkx2m(k)=k-2;
            kx2m(k)=-0.5*sqrt(nrad*(nrad-1))*(rhos0/delwave(k))^2;
         end
         if(iwaverad<=(nwaverad-2))
            kkx2p(k)=k+2;
            kx2p(k)=-0.5*sqrt((nrad+1)*(nrad+2))*(rhos0/delwave(k))^2;
         end
       end
      end

      for k=1:kmax
        ky(k)=rhos0*mf(k)/r(jss(k)); 
      end

      % for hyper-viscous damping, radial part of Delp4 is omitted.
      for k=1:kmax
        d0kx2(k) =xmu1*kx2(k); 
        d0kx2p(k)=xmu1*kx2p(k);
        d0kx2m(k)=xmu1*kx2m(k);
        xmu(k)   =xmu1*ky(k)^2 + xmu2*ky(k)^4;
      end
      
   else % basis = bessel function
      d0kx2=0*ones(1,kmax);
      for k=1:kmax
        d0kx2(k) =xmu1*kperp(k)^2 + xmu2*kperp(k)^4; 
      end
   end

%----------------------------------------------------------------------
%  Toroidal coupling via magnetic curvature 
%----------------------------------------------------------------------
status='   Calculating magnetic curvature terms'

      auxatmp=0*ones(jmax,kmax);
      auxbtmp=0*ones(jmax,kmax);
      apwx=0*ones(jmax,kmax);
      amwx=0*ones(jmax,kmax);
      bpwx=0*ones(jmax,kmax);
      bmwx=0*ones(jmax,kmax);
      apwxp=0*ones(1,jmax);
      amwxp=0*ones(1,jmax);
      bpwxp=0*ones(1,jmax);
      bmwxp=0*ones(1,jmax);

      for k=1:kmax
        m=mf(k);
        for j=1:jmax
          auxatmp(j,k)=rhos0/rmajor*fac(j)*drm(j)*wave(j,k);
          auxbtmp(j,k)=ar(j)*auxatmp(j,k);
          apwx(j,k)=enhat(j)*((m+1)*(1+tau(j))+(rln(j)+rltau(j))*tau(j));
          amwx(j,k)=enhat(j)*((m-1)*(1+tau(j))-(rln(j)+rltau(j))*tau(j));
          bpwx(j,k)=enhat(j)*(m+1+rln(j));
          bmwx(j,k)=enhat(j)*(m-1-rln(j));
        end
      end

      for j=1:jmax
          apwxp(j)=enhat(j)*(1+tau(j));
          amwxp(j)=-enhat(j)*(1+tau(j));
          bpwxp(j)=enhat(j);
          bmwxp(j)=-enhat(j);
      end
%----------------------------------------------------------------------
      apw=0*ones(kmax,nwaverad);
      amw=0*ones(kmax,nwaverad);
      bpw=0*ones(kmax,nwaverad);
      bmw=0*ones(kmax,nwaverad);

      for k=1:kmax
        ktmp_apw=k;
        %here ktmp_apw
        for iwaveradp=1:nwaverad
          apw0=0;
          amw0=0;
          bpw0=0;
          bmw0=0;

          kp1=kpw(k,iwaveradp);
          if(kp1~=0)
            for j=1:jmax
               auxp=auxatmp(j,k)*wave(j,kp1);
               auxpp=auxbtmp(j,k)*dwavedr(j,kp1);
               apw0=apw0+apwx(j,k)*auxp+apwxp(j)*auxpp;
               bpw0=bpw0+bpwx(j,k)*auxp+bpwxp(j)*auxpp;
            end
          end

          km1=kmw(k,iwaveradp);
          if(km1~=0)
            for j=1:jmax
               auxm=auxatmp(j,k)*wave(j,kmw(k,iwaveradp));
               auxmm=auxbtmp(j,k)*dwavedr(j,kmw(k,iwaveradp));
               amw0=amw0+amwx(j,k)*auxm+amwxp(j)*auxmm;
               bmw0=bmw0+bmwx(j,k)*auxm+bmwxp(j)*auxmm;
            end
          end

          apw(k,iwaveradp)=apw0;
          amw(k,iwaveradp)=amw0;
          bpw(k,iwaveradp)=bpw0;
          bmw(k,iwaveradp)=bmw0;
         end
      end
      
      clear auxatmp auxbtmp apwx amwx bpwx bmwx 
      clear apwxp amwxp bpwxp bmwxp 

status='****** Initialization done ******'
