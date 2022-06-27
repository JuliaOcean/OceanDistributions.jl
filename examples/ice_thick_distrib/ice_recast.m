function [V,C,H]=ice_recast(v,c,h);
%object: recast ice categories after thermodynamical ice growth/melt
%	(v increments) may result in ice categories outgrowing
%	their limits ( v/c>h(n) or v/c<h(n-1) ). Here the 
%	recasting to smaller categories is done synchroneously for 
%	all cats., and then the recasting to larger cats. similarly.
%inputs:	v is the ice volume in each cat
%		c is the ice concentration in each cat
%		h is the cat thickness limits
%outputs:	V,C,H are the recasted vol.,conc., and thick.
%note:	below, for growing ice, the fraction of ice cover getting
%	recasted (0<beta<=1) can be changed (default is beta=1).

V=v; C=c; H=V./c;
nh=size(V,1);
nx=size(V,2);

%1) recast to thinner categories

Vmin=c(3:end,:).*h(2:end-1,:);
V0=v(3:end,:);
C0=c(3:end,:);

dV0=V0.*(V0<Vmin); 
V(2:end-1,:)=V(2:end-1,:)+dV0;
V(3:end,:)=V(3:end,:)-dV0;

dC0=C0.*(V0<Vmin); 
C(2:end-1,:)=C(2:end-1,:)+dC0;
C(3:end,:)=C(3:end,:)-dC0;

%2) recast to thicker categories

Vmax=c(2:end-1,:).*h(2:end-1,:);
V1=v(2:end-1,:);
C1=c(2:end-1,:);

beta=1*ones(size(V1));
% beta=0.5*ones(size(V1));
% beta(1:10,:)=0.25+0.75*([1:10]'-1)/(10)*ones(1,nx);
V2=V1+(1./beta-1).*(V1-Vmax);
alpha=beta.*V2./V1;
beta(V1==0|V1<=Vmax)=0;
alpha(V1==0|V1<=Vmax)=0;

dV1=alpha.*V1;
V(2:end-1,:)=V(2:end-1,:)-dV1;
V(3:end,:)=V(3:end,:)+dV1;

dC1=beta.*C1; 
C(2:end-1,:)=C(2:end-1,:)-dC1;
C(3:end,:)=C(3:end,:)+dC1;

