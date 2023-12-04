% layers  calculate admittance and internal field in layered tranducer
%         structure
%
% Usage:  [y,ua,ub]=layers(f,ml,mt,mr)
% 
% Input:    f   frequency 
%           ml  materials on "left" side, ml=[m0;m1;...;mN1]
%           mt  transducer structure, mt=[mN1+1;mN1+2;...;mN2-1]
%           mr  materials on "right" side, mr=[mN2;mN2+1;...;mN-1] 
%               each row mn, n in {0,1,...,N-1} of the material specification
%               has the format
%               mn=[layer thickness, velocity, impedance, loss, ...
%                   relative permittivity, conductivity, ... 
%                   electromechanical coupling] 
% 
% Output:   y   admittance 
%           ua  rightmoving wave amplitudes per unit voltage 
%           ub  leftmoving wave amplitudes per unit voltage
%
% See also: dfield 
%
% Written by Einar Halvorsen, Vestfold University College, 2006
%
% Changelog
% Version Date         Author   Comment
% 1.00    Sep 27, 2006 EH       Initial  version
% 1.01    Sep 25, 2010 EH       Misprint in header
% 1.05    Feb  8, 2011 EH       Changed sign error in matrix elements for
%                               piezoelectric coupling, affected multilayered mt
%                               Cleaned up comments at the end 
%                               Expanden documentation
function [y,ua,ub]=layers(f,ml,mt,mr)
dofield=nargout>1;
sizef=size(f);
f=f(:);
y=zeros(sizef);
eps0=8.854187816e-12;
NF=length(f);
N1=size(ml,1)-1;
N2=N1+size(mt,1)+1;
N=N2+size(mr,1);
NL=size(ml,1);
NR=size(mr,1);
NT=size(mt,1);
ND=2*NL+2*NT+2*NR+1;
if dofield
   ua=zeros(NF,NL+NT+NR);
   ub=zeros(NF,NL+NT+NR);
end
mall=[ml; mt; mr];
for nf=1:length(f)
   w=2*pi*f(nf);
   %%% build system matrix
   sysm=zeros(ND,ND);
   % first only the mechanical dofs
   kd=mall(1,4)*w*w;
   k=w/mall(1,2)-i*kd;
   if (abs(imag(k))>0.1*abs(real(k)))
      fprintf(1,'warning: extreme losses reconsider model');
   end
   if isinf(mall(1,1))
      pf=0;
   else
      pf=exp(-i*k*mall(1,1));
   end
   sysm(1,1)=1;
   sysm(1,2)=-pf;
   pzm=ml(1,3);
   ppf=pf;
   for n=1:N-1
      zm=mall(1+n,3);
      kd=mall(1+n,4)*w*w;
      k=w/mall(1+n,2)-i*kd;
      if isinf(mall(1+n,1))
         pf=0;
      else
         pf=exp(-i*k*mall(1+n,1));
      end
      n1=2*n;
      n2=2*n+1;
      sysm(n1,n2-2)=ppf;
      sysm(n1,n2-1)=1;
      sysm(n1+1,n2-2)=-pzm*ppf;
      sysm(n1+1,n2-1)=pzm;
      sysm(n1,n2)=-1;
      sysm(n1,n2+1)=-pf;
      sysm(n1+1,n2)=zm;
      sysm(n1+1,n2+1)=-zm*pf;
      pzm=zm;
      ppf=pf;
   end
   n1=n1+2;
   sysm(n1,n2)=pf;
   sysm(n1,n2+1)=-1;
   % then fill in elements associated with electrical dofs
   % from transducer left boundary
   n1=2*N2;
   ke=mall(N2,7);
   km=sqrt( mall(N2,3)*mall(N2,2)/(eps0*mall(N2,5)) )*ke;
   sysm(n1+1,end)=-i*km/w;
   n1=2*(N1+1);
   ke=mall(N1+2,7);
   km=sqrt( mall(N1+2,3)*mall(N1+2,2)/(eps0*mall(N1+2,5))  )*ke;
   sysm(n1+1,end)=i*km/w;
   pkm=km;
   for n=(N1+2):(N2-1)
      n1=2*n;
      ke=mall(n+1,7);
      km=sqrt( mall(n+1,3)*mall(n+1,2) / (eps0*mall(n+1,5)) )*ke;
      sysm(n1+1,end)=i*(km-pkm)/w; % corrected in v1.05
      pkm=km;
   end
   iepsacc=0;
   kmmax=-inf;
   for n=(N1+1):(N2-1)
      iepsacc=iepsacc+mall(n+1,1)/(eps0*mall(n+1,5));
      n2=2*n+1;         
      kd=mall(n+1,4)*w*w;
      k=w/mall(n+1,2)-i*kd;
      pf=exp(-i*k*mall(1+n,1));
      ke=mall(n+1,7);
      km=sqrt(  mall(n+1,3)*mall(n+1,2)/(eps0*mall(n+1,5)) )*ke;
      kmmax=max(kmmax,km);
      sysm(end,n2)=km*(pf-1);
      sysm(end,n2+1)=km*(1-pf);
   end
   sysm(end,end)=iepsacc;
   %% scale 
   zscale=max(mall(:,3));
   sysm(3:2:(end-1),:)=sysm(3:2:(end-1),:)/zscale;
   sysm(end,1:(end))=sysm(end,1:(end))/kmmax;
   chfact=zscale*w/kmmax;
   sysm(:,end)=chfact*sysm(:,end);
   %% solve
   x=zeros(ND,1);
   x(end)=1/kmmax;
   sol=sysm\x;
   y(nf)=i*w*chfact*sol(end);
   ua(nf,:)=sol(1:2:(end-1));
   ub(nf,:)=sol(2:2:(end-1));   
end


