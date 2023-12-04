% dfield  calculate displacement field in multilayer structure
function u=dfield(m,x,ua,ub,f)
x=x(:); f=f(:);
u=NaN*ones(length(x),length(f));
lw=m(:,1);
nl=length(lw);
l1=zeros(size(lw)); l2=l1;
if isinf(lw(1))
   l1(1)=-Inf;
   l2(1)=0;
else
   l1(1)=0;
   l2(1)=lw(1);
end
for n=2:nl
   l1(n)=l2(n-1);
   l2(n)=l1(n)+lw(n);
end
for nf=1:length(f)
   w=2*pi*f(nf);
   for n=1:nl
      indx=find(x>=l1(n) & x<= l2(n));
      xn=x(indx);
      k=w/m(n,2);
      alph=m(n,4)*w*w;
      el=exp(-(i*k+alph)*(xn-l1(n)));
      er=exp((i*k+alph)*(xn-l2(n)));
      ux=ua(nf,n)*el+ub(nf,n)*er;
      u(indx,nf)=ux;
   end
end
