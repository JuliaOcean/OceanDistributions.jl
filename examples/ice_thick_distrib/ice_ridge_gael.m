function [g2,v2]=ice_ridge_bitz(g1,v1);
%object : remove excess area (sum(g1,1)-1>0) in a case where
%           we cant further open leads (g1(1,.)=0) so
%           we need to raise ridges (while conserving volume)
%input: g1 is the thickness distribution (pdf.*dh)
%       h is the corresponding thickness
%note:  the first direction is thickness space
%       additional directions are for geographic space

% fprintf('this should to be re-coded to support sharp distributions\n');
% fprintf('e.g. thermo+advdiff run has failed when all area is in layer 2');

% j=20;
nh=size(g1,1);
dG=sum(g1)-1;
%dG=dGdtforRidges*dt;

g2=g1; v2=v1;

%initialize conservation tests:
test1=sum(v2,1)'./sum(v1,1)';
test2=sum(g2,1)';
test3=NaN*test1; test3(abs(dG)<=1e-8)=0;

la=0*dG; mu=la;
kk=find(la==0&abs(dG)>1e-8);
j=1;
while ~isempty(kk)&j<nh;
    j=j+1;

    a=sum(g1(1:j,kk),1);
    b=sum(g1(j+1:end,kk),1);
    A=sum(v1(1:j,kk),1);
    B=sum(v1(j+1:end,kk),1);
    D=dG(kk);
    
    %mask out points where no ice for i<=j
    m=(b.*A-a.*B); m=1*(abs(m)>1e-8); m(find(~m))=NaN;
    a=a.*m; b=b.*m; A=A.*m; B=B.*m;
    
    la=B./(b.*A-a.*B).*D;
    mu=-A./(b.*A-a.*B).*D;
    %treat case when e.g. one thickness only
    ii=find(~isfinite(la)); la(ii)=0; mu(ii)=0;
    %treat the case of too little ice until j
    ii=find(la<-1); la(ii)=0; mu(ii)=0;    
    
    g2(1:j,kk)=g1(1:j,kk).*(1+ones(j,1)*la);
    g2(j+1:end,kk)=g1(j+1:end,kk).*(1+ones(nh-j,1)*mu);

    v2(1:j,kk)=v1(1:j,kk).*(1+ones(j,1)*la);
    v2(j+1:end,kk)=v1(j+1:end,kk).*(1+ones(nh-j,1)*mu);
    
    %update conservation tests
    test1=sum(v2,1)'./sum(v1,1)';
    test2=sum(g2,1)';
    %keep track : completed ridging
    tmp1=find(la<0&abs(D)>1e-8);
    test3(kk(tmp1))=j;
    %iterate : left for ridging
    tmp1=find(la==0&abs(D)>1e-8);
    kk=kk(tmp1);
    
end;

%if g is a very sharp distribution, then the 
%above mail fail to find a working j; in this case
%we simply reduce conc (rasing h) homogeneously; and
%will complement that by moving ice to appropriate h.
if ~isempty(kk);
    fac=1./(1+dG(kk));
    g2(:,kk)=g1(:,kk).*(ones(j,1)*fac);
    %update conservation tests
    test1=sum(v2,1)'./sum(v1,1)';
    test2=sum(g2,1)';
    %keep track : completed ridging
    tmp1=find(la<0&abs(D)>1e-8);
    test3(kk(tmp1))=-j;
    kk=[];
end;



err=0;

if ~isempty(kk);
    fprintf('failed to ridge\n');
    err=1;
end;

if max(abs(test1-1))>1e-8&~isnan(test3);
    fprintf('failed to conserve volume\n');
    err=2;
end;

if max(abs(test2-1))>1e-8&~isnan(test3);
    fprintf('failed to conserve area\n');
    err=3;
end;

if err>0; 
    error('ice_ridge_gael : abort %i \n',err);
end;
