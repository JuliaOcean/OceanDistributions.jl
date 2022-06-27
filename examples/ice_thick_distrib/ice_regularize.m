function [V,C]=ice_regularize(V,C,it);
%object : get rid of negatives and small positives
%   that will be moved to other existing ice classes
%   if possible or to open water otherwise

global nh nt nx;
global reg2ice reg2ocn;
if it==1; reg2ice=zeros(nt+1,2); reg2ocn=zeros(nt+1,2); end;

tmp1=V(2:end,:);
tmp2=C(2:end,:);
%identify small pos:
tmp3=1*(tmp1>=1e-6&tmp2>=1e-6);
%compute the small pos and large pos reservoirs:
tmp1_large=sum(tmp1.*tmp3,1);
tmp1_small=sum(tmp1.*(1-tmp3),1);
tmp2_large=sum(tmp2.*tmp3,1);
tmp2_small=sum(tmp2.*(1-tmp3),1);
%if possible then merge with other cats:
ii=find(tmp1_large>100*abs(tmp1_small)&tmp2_large>100*abs(tmp2_small));
tmp11=(tmp1_large+tmp1_small)./tmp1_large;
tmp11=ones(size(tmp1,1),1)*tmp11;
V(2:end,ii)=(tmp1(:,ii).*tmp3(:,ii)).*tmp11(:,ii);
tmp22=(tmp2_large+tmp2_small)./tmp2_large;
tmp22=ones(size(tmp2,1),1)*tmp22;
C(2:end,ii)=(tmp2(:,ii).*tmp3(:,ii)).*tmp22(:,ii);
%otherwise dump in open water:
ii=find(tmp1_large<=100*abs(tmp1_small)|tmp2_large<=100*abs(tmp2_small));
V(2:end,ii)=tmp1(:,ii).*tmp3(:,ii);
C(2:end,ii)=tmp2(:,ii).*tmp3(:,ii);
V(1,ii)=V(1,ii)+tmp1_small(ii);%dump in open water
C(1,ii)=C(1,ii)+tmp2_small(ii);%dump in open water

%keep track of what has been merged with other cats:
ii=find(tmp1_large>100*abs(tmp1_small)&tmp2_large>100*abs(tmp2_small));
tmp1=nansum(tmp1_small(ii));
tmp2=nansum(tmp2_small(ii));
if tmp1~=0|tmp2~=0;
    %         fprintf('%d pos merged %g %g\n',it,tmp1,tmp2);
end;
reg2ice(it,:)=reg2ice(it,:)+[tmp1 tmp2];
%keep track of what has been dumped to open water:
ii=find(tmp1_large<=100*abs(tmp1_small)|tmp2_large<=100*abs(tmp2_small));
tmp1=nansum(tmp1_small(ii));
tmp2=nansum(tmp2_small(ii));
if tmp1~=0|tmp2~=0;
    %         fprintf('%d pos removed %g %g\n',it,tmp1,tmp2);
end;
reg2ocn(it,:)=reg2ocn(it,:)+[tmp1 tmp2];

