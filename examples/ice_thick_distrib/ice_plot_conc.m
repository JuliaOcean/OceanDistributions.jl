function []=ice_plot_conc(v,c,h);

aa=8; bb=4;
nh=size(v,1)-1;
nx=size(v,2)-1;
xplot=ones(nh+2,1)*[1:nx+2];
yplot=[0:nh+1]'*ones(1,nx+2);
% zplot=NaN*xplot; zplot(1:nh+1,1:nx+1)=log10(abs(c));
zplot=NaN*xplot; zplot(1:nh+1,1:nx+1)=c; zplot(zplot==0)=NaN;

figure; set(gcf,'Units','Normalized','Position',[0.1 0.3 0.4 0.6]);
for ii=1:4;
%     if ii>2; AA=[1:aa]; else; AA=[nh+1-aa:nh+1]; end;
    if ii>2; AA=[1:aa]; else; AA=[64-aa+1:64]-38; end;
    if mod(ii,2)==1; BB=[1:bb]; else; BB=[nx+1-bb+1:nx+1]; end;
    subplot(2,2,ii);
    %     imagescnan(zplot(AA,BB));
    pcolor(xplot(AA,BB),yplot(AA,BB),zplot(AA,BB));
%     caxis([-2 -1.5]);
    caxis([0 0.04]);
    caxis([-1 1]*1e-8);
end;
