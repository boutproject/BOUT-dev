%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour plot of eigenmode structure for potential  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=-pi:0.01:pi;
mft=mf';
matexp=exp(i*mft(1:kmax)*theta);
wavemap=wave(:,1:kmax);
rcore=r';
x=rcore*cos(theta);
y=rcore*sin(theta);

%for n1=nmin:ndel:nmax

n1=10; % choose toroidal mode number

  nn=n1/ndel;
  figure(2+nn)
  ph=0*ones(kmax,1);
  l=0;
  for k=1:kmax
    if(nfb(k)==n1)
       l=l+1;
       ph(k)=Fphi(nn,l);
     end
  end
  matph=ph(1:kmax,1)*ones(1,length(theta));
  matexp1=matph.*matexp;
  phi=wavemap*matexp1;
  phi=real(phi);
  pcolor(x,y,phi);
  shading('flat');
  colormap('jet');
  colorbar('vert');
  xlabel('r/a')
  str1=['Eigenmode stucture for n=' int2str(n1)];
  title(str1);
  
%end
