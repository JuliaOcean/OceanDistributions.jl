function [V,C,H]=ice_remap(v,c,h);
%(Lipscomb's approach using linear distributions; simplified.)
%
%object: remap ice categories after thermodynamical ice growth/melt
%	(v increments) may result in ice categories outgrowing
%	their limits ( v/c>h(n) or v/c<h(n-1) ). Here we assume a
%   homogeneous distribution (linear more generally; see Lipscomb 01)
%   around the updated h1=v/c and compute the overlapping integral 
%   within the prescribed neighbor category limits. The result (dC&dV) 
%   is then moved to that neighboring category (thicker or thinner).
%inputs:	v is the ice volume in each cat
%		c is the ice concentration in each cat
%		h is the cat thickness limits
%outputs:	V,C,H are the recasted vol.,conc., and thick.

V=v; C=c; H=V./c;
nh=size(V,1);
nx=size(V,2);

%0) invert formula for new distribution
%======================================

dh=0.25;

c1=c(2:end-1,:);%category concentration
v1=v(2:end-1,:);%category volume
h1=v1./c1;%moved category center
HL1=h1-dh/2;%moved lower limit
HR1=h1+dh/2;%moved upper limit

HL0=h(1:end-2,:);%fixed lower limit
HR0=h(2:end-1,:);%fixed upper limit
h0=(HL0+HR0)/2;%fixed category center

etaR=HR1-HL1;
etaN=h1-HL1;
g1=12*c1./(etaR.^3).*(etaN-etaR/2);
g0=6*c1./(etaR.^2).*(2/3*etaR-etaN);
%new distribution is 0 outside [HL1 HR1] and otherwise is
%g(h)=g0+g1*(h-HL1);
%
%note: as compared with Lipscomb 01, I simplified the 
%computation of lower/upper boundaries, which implies that
%the distributions are homogeneous -- so g(h)=g0 in fact.
%after simplification:
g1=0*ones(size(etaR)); g0=c1/dh;

%1) recast to thinner categories
%===============================

%1.1) compute integrals between HR0 and HR1, which
%   will be transfered to the next category
dC0=(g0-g1.*HL1).*(HL0-HL1)+1/2*g1.*(HL0.^2-HL1.^2);
dV0=1/2*(g0-g1.*HL1).*(HL0.^2-HL1.^2)+1/3*g1.*(HL0.^3-HL1.^3);

%avoid leaving exponentially small numbers : full recasting
tmp1=V(2:end-1,:); tmp2=C(2:end-1,:);
ii=find(tmp1<1e-4|tmp2<1e-4);
dV0(ii)=tmp1(ii); dC0(ii)=tmp2(ii);

%avoid category jump (HL1<HL0-dh) that leads to 
%dC0>C or dV0>V and ultimately to negative values
tmp1=V(2:end-1,:); tmp2=C(2:end-1,:);
ii=find(HL1<HL0-dh);
dV0(ii)=tmp1(ii); dC0(ii)=tmp2(ii);

%1.2) mask out points where no remapping
ii=find(HL0==0|HL1>=HL0|c1==0|v1==0); dC0(ii)=0; dV0(ii)=0;

%1.3) time step V and C:
V(2:end-1,:)=V(2:end-1,:)-dV0;
V(1:end-2,:)=V(1:end-2,:)+dV0;

C(2:end-1,:)=C(2:end-1,:)-dC0;
C(1:end-2,:)=C(1:end-2,:)+dC0;

%2) recast to thicker categories
%===============================

%2.1) compute integrals between HR0 and HR1, which
%   will be transfered to the next category
dC1=(g0-g1.*HL1).*(HR1-HR0)+1/2*g1.*(HR1.^2-HR0.^2);
dV1=1/2*(g0-g1.*HL1).*(HR1.^2-HR0.^2)+1/3*g1.*(HR1.^3-HR0.^3);

%avoid leaving exponentially small numbers : full recasting
tmp1=V(2:end-1,:); tmp2=C(2:end-1,:);
ii=find(tmp1<1e-4|tmp2<1e-4);
dV1(ii)=tmp1(ii); dC1(ii)=tmp2(ii);

%avoid category jump (HR1>HR0-dh) that leads to 
%dC0>C or dV0>V and ultimately to negative values
tmp1=V(2:end-1,:); tmp2=C(2:end-1,:);
ii=find(HR1>HR0+dh);
dV0(ii)=tmp1(ii); dC0(ii)=tmp2(ii);

%2.2) mask out points where no remapping
ii=find(HR1<=HR0|c1==0|v1==0); dC1(ii)=0; dV1(ii)=0;

%2.3) time step V and C:
V(2:end-1,:)=V(2:end-1,:)-dV1;
V(3:end,:)=V(3:end,:)+dV1;

C(2:end-1,:)=C(2:end-1,:)-dC1;
C(3:end,:)=C(3:end,:)+dC1;
    
