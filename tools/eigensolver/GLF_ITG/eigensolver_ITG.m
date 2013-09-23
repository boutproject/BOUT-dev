%======================================================================
% Solve eigenvalue equation for ITG : 
% Ottaviani model used [M. Ottaviani et.al., Phys. Plasmas 6, 3267 (1999)]
%======================================================================
  
  gyroaverage=1; % turn-on(1) or off(0) gyroaverage
  damping_glf=1; % turn-on(1) or off(0) Landau damping term
  damping_c=0;   % turn-on(1) or off(0) viscous damping term

  l_k=0*ones(1,kmax);
  lmaxn=0*ones(1,nmax/ndel);
  
  l=0;
  for k=1:kmax
     if(nfb(k)==nmax)
       l=l+1;
     end
  end
  lmax=l;
  Fphi=0*ones(nmax/ndel,lmax);
  
status='******  Start eigenvalue solver ******'
%----------------------------------------------------------------------
  for n=nmin:ndel:nmax % toroidal mode number scan
%----------------------------------------------------------------------
      nn=n/ndel;
      toroidal_mode_number = n
      
    % Determine matrix size for n mode 
      l=0;
      for k=1:kmax
          if(nfb(k)==n)
            l=l+1;
          end
      end
      lmaxn(nn)=l;
      lmax=lmaxn(nn);
       
    % Match index (k) of basis to  matrix index (l) 

      k_l=0*ones(1,lmax);
      l=0;
      for k=1:kmax
          if(nfb(k)==n)
            l=l+1;
            k_l(l)=k;
            l_k(k)=l;
          end
      end
      
%----------------------------------------------------------------------      
%    ExB, diamagnetic flow terms 
%----------------------------------------------------------------------
status='   Calculating matrix for ExB, diamagnetic flow terms'

      G_pi=0*ones(lmax,lmax);  % -i*G_pi*Lphi = Vpi*Grad(Lphi), Vpi=diamagnetic flow velocity, L=Laplacian
      G_ti=0*ones(lmax,lmax);  % i*G_ti*phi = Vexb*Grad(Ti), Vexb=ExB flow velocity
      G_ni=0*ones(lmax,lmax);  % i*G_ni*phi = Vexb*Grad(n)
      for l=1:lmax
        k=k_l(l);
        btmp=-rhos0*mf(k);
        for l1=1:lmax
         k1=k_l(l1);
         if(mf(k)==mf(k1))
          for j=1:jmax
            atmp=wave(j,k)*wave(j,k1);
            G_pi(l,l1)=G_pi(l,l1)+atmp*dpihat(j);
            G_ti(l,l1)=G_ti(l,l1)+atmp*dthat(j);
            G_ni(l,l1)=G_ni(l,l1)+atmp*denhat(j);
          end
          G_pi(l,l1)=btmp*G_pi(l,l1);
          G_ti(l,l1)=btmp*G_ti(l,l1);
          G_ni(l,l1)=btmp*G_ni(l,l1);
         end
        end
      end
      
%----------------------------------------------------------------------
%    Viscous damping, Laplacian, <J0>=Gamma0^(1/2) 
%----------------------------------------------------------------------
status='   Calculating matrix for Laplacian, gyroaverage, viscous damping operator'

      Dc=0*ones(lmax,lmax); % viscous damping operator
      matL=0*ones(lmax,lmax); % Laplacian = Delp2
      J0=0*ones(lmax,lmax); % gyroaveraging operator

    if(basis==hermite)
      for l=1:lmax
         k=k_l(l);
         j=jss(k);
         kxp=kkx2p(k);
         if(kxp~=0) 
             lxp=l_k(kxp);
             Dc(l,lxp)=-d0kx2p(k);
             matL(l,lxp)=-kx2p(k);
         end
         kxm=kkx2m(k);
         if(kxm~=0)
             lxm=l_k(kxm);
             Dc(l,lxm)=-d0kx2m(k);
             matL(l,lxm)=-kx2m(k);
         end
         Dc(l,l)=-(d0kx2(k)+xmu(k));
         matL(l,l)=-(ky(k)^2+kx2(k));
      end
      J0=inv(1-gyroaverage*0.5*matL);

    else % basis = bessel 
      for l=1:lmax
         k=k_l(l);
         Dc(l,l)=-d0kx2(k);
         matL(l,l)=-kperp(k)^2;
         J0(l,l)=1/(1-gyroaverage*0.5*matL(l,l));
      end
    end
    
      Dc=damping_c*Dc;
      
%----------------------------------------------------------------------
%   Parallel wavenumber & Landau damping
%----------------------------------------------------------------------
status='   Calculating matrix for parallel gradient and Landau damping operator'

      matakpar=0*ones(lmax,lmax); % kpar(parallel wavenumber) = -i*Grad_par
      matakpar_enhat=0*ones(lmax,lmax); % density*kpar
      matakpar_that=0*ones(lmax,lmax); % Ti*kpar
      matakpar_tau=0*ones(lmax,lmax);  % Ti/Te*kpar
      Dglf=0*ones(lmax,lmax); % Landau damping operator = -sqrt(8*Ti/mi/pi)*abs(kpar)
      for l=1:lmax
         k=k_l(l);
         for iwaverad=1:nwaverad 
            l0=l_k(kw(k,iwaverad));
            matakpar(l,l0)=akpar(k,iwaverad);
            matakpar_enhat(l,l0)=akpar_enhat(k,iwaverad);
            matakpar_tau(l,l0)=akpar_tau(k,iwaverad);
            matakpar_that(l,l0)=akpar_that(k,iwaverad);
            Dglf(l,l0)=-kparglf_sqrtthat(k,iwaverad);
         end
      end

      Dglf=damping_glf*Dglf;
      
%----------------------------------------------------------------------
%   Vorticity = M * phi
%----------------------------------------------------------------------
status='   Calculating matrix M for voriticity(= M * potential)'

      matM=0*ones(lmax,lmax); % Vorticity=(ne/Te-ne*Delp2)phi = M*phi
      matMinv=0*ones(lmax,lmax); % inverse of M

  if(basis==hermite)
      for l=1:lmax
         k=k_l(l);
         j=jss(k);
         kxp=kkx2p(k);
         if(kxp~=0) 
             lxp=l_k(kxp);
             matM(l,lxp)=enhat(j)*kx2p(k);
         end
         kxm=kkx2m(k);
         if(kxm~=0)
             lxm=l_k(kxm);
             matM(l,lxm)=enhat(j)*kx2m(k);
         end
         matM(l,l)=enhat(j)/tehat(j)*(1+tehat(j)*(ky(k)^2+kx2(k)));
      end
      matMinv=inv(matM);    

  else % basis = bessel
      for l=1:lmax
         k=k_l(l);
         j=jss(k);
         matM(l,l)=enhat(j)/tehat(j)*(1+tehat(j)*kperp(k)^2);
         matMinv(l,l)=1/matM(l,l);
      end
  end
      
%----------------------------------------------------------------------
%  Coupling between m and m+1,m-1 modes via magnetic curvature 
%----------------------------------------------------------------------
status='   Calculating matrix for toroidal coupling via magnetic curvature'

      aw=0*ones(lmax,lmax);
      bw=0*ones(lmax,lmax);
      for l=1:lmax
         k=k_l(l);
         j=jss(k);
         for iwaverad=1:nwaverad 
            kp1=kpw(k,iwaverad);
            if(kp1~=0) 
               lp=l_k(kp1);
               aw(l,lp)=apw(k,iwaverad);
               bw(l,lp)=bpw(k,iwaverad);
            end

            km1=kmw(k,iwaverad);
            if(km1~=0)
               lm=l_k(km1);
               aw(l,lm)=amw(k,iwaverad);
               bw(l,lm)=bmw(k,iwaverad);
            end
         end
      end
      
      
%----------------------------------------------------------------------      
%   Matrix for eigenvalue equation
%----------------------------------------------------------------------
status='   Setting up matrix for eigenvalue equation'

%--------------------------------------------------
%               |A11 A12 A13|    |phi |
%  |AA|*|F| =   |A21 A22 A23|  * |Ti  | = omega*|F|
%               |A31 A32 A33|    |Vpar|
%---------------------------------------------------

      A11=0*ones(lmax,lmax);
      A12=0*ones(lmax,lmax);
      A13=0*ones(lmax,lmax);
      A21=0*ones(lmax,lmax);
      A22=0*ones(lmax,lmax);
      A23=0*ones(lmax,lmax);
      A31=0*ones(lmax,lmax);
      A32=0*ones(lmax,lmax);
      A33=0*ones(lmax,lmax);

      A11= matMinv* ( G_pi*matL + G_ni*J0 + i*Dc*matM - aw );
      A12= -matMinv* bw;
      A13= matMinv* matakpar_enhat;

      A21= G_ti*J0;
      A22= i*Dc + i*Dglf;
      A23=2/3*matakpar_that;       

      A31= matakpar + matakpar_tau;
      A32= matakpar;
      A33= i*Dc;

      clear G_pi G_ti G_ni 
      clear Dc Dglf matM matMinv matL J0
      clear matakpar matakpar_enhat matakpar_that matakpar_tau
      clear aw bw 

      AA=0*ones(3*lmax,3*lmax);

      AA(1:lmax,1:lmax)=A11(:,:);
      AA(lmax+1:2*lmax,1:lmax)=A21(:,:);
      AA(2*lmax+1:3*lmax,1:lmax)=A31(:,:);

      AA(1:lmax,lmax+1:2*lmax)=A12(:,:);
      AA(lmax+1:2*lmax,lmax+1:2*lmax)=A22(:,:);
      AA(2*lmax+1:3*lmax,lmax+1:2*lmax)=A32(:,:);

      AA(1:lmax,2*lmax+1:3*lmax)=A13(:,:);
      AA(lmax+1:2*lmax,2*lmax+1:3*lmax)=A23(:,:);
      AA(2*lmax+1:3*lmax,2*lmax+1:3*lmax)=A33(:,:); 

      clear A11 A12 A13 A21 A22 A23 A31 A32 A33      

%----------------------------------------------------------------------      
%   Obtain eigenvalues & eigenvectors
%----------------------------------------------------------------------
status='   Calculating eigenvalues and eigenvectors'

      F=0*ones(3*lmax,3*lmax);  % eigenvector
      omg=0*ones(3*lmax,3*lmax); 
 
      [F,omg]=eig(AA);  % solve eigenvalue equation
      
      omgvec=0*ones(3*lmax,1);  % real part of eigenvalue
      gamvec=0*ones(3*lmax,1);  % imaginary part of eigenvalue
      
      for l=1:3*lmax
        omgvec(l)=real(omg(l,l));
        gamvec(l)=imag(omg(l,l));
      end
      [gam,index]=max(gamvec);  % find maximal growth rate
      gamlin(nn)=gam;
      omglin(nn)=omgvec(index); 

      % eigenvector for potential
      for l=1:lmax
        Fphi(nn,l)=F(l,index);
      end
%----------------------------------------------------------------------   
     clear AA F omg omevec gamvec
%----------------------------------------------------------------------   
  end
  
status='****** Eigenvalue solver finished ******'
  


