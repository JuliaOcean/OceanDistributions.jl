
%main model parameters:
%======================

doThorndike=0; doBitz=0; doLipscomb=1; suff='lips';
% doThorndike=0; doBitz=1; doLipscomb=0; suff='bitz';
% doThorndike=1; doBitz=0; doLipscomb=0; suff='thor';

choiceIC=2;
choiceAdv=1;
choiceDiff=0;
doThermo=1;
doHeatCons=1;
displayChecks=1;

%test for parameter inconsistency:
if ~doLipscomb&~doBitz&~doThorndike; error('choose Lipscomb OR Bitz OR Thorndike'); end;
if doLipscomb&doBitz; error('choose Lipscomb OR Bitz OR Thorndike'); end;
if doThorndike&doBitz; error('choose Lipscomb OR Bitz OR Thorndike'); end;
if doThorndike&doLipscomb; error('choose Lipscomb OR Bitz OR Thorndike'); end;
if doThorndike&doHeatCons; error('doHeatCons AND doThorndike : un-supported'); doHeatCons=0; end;

%print parameters to screen:
fprintf('==== model parameters: ====\n')

if ~doThermo;
    fprintf('thermal forcing :  NO\n');
else;
    fprintf('thermal forcing :  YES\n');
    if doLipscomb;
        fprintf(['    model :        Lipscomb \n']);
    elseif doBitz;
        fprintf(['    model :        Bitz \n']);
    elseif doThorndike;
        fprintf(['    model :        Thorndike \n']);
    else;
        fprintf(['    model :        (none) \n']);
    end;
end;

if doHeatCons; tmp='YES'; else; tmp='NO'; end;
fprintf(['    conserve :     ' tmp '\n']);

if choiceAdv==0; tmp='(none)';
elseif choiceAdv==1; tmp='to the right';
elseif choiceAdv==2; tmp='randomly';
end;
fprintf(['advect  :          ' tmp '\n']);
if choiceAdv>0; tmp='YES'; else; tmp='NO'; end;
fprintf(['diffuse :          ' tmp '\n']);

fprintf('====== model run ... ======\n')


%plotting parameters:
%====================
freqPlot=NaN;
freqPlot=15;
pausePlot=1e-4;
posPlot=[0.5 0.4 0.5 0.5];

%grid parameters:
%================
global nh nt nx;

dt=0.1; nt=2000;
% dh=0.25; nh=160;
dh=0.25; nh=120;
% dh=0.05; nh=600;
dx=0.2; nx=50;

h=[0 0.5:nh-0.5]*dh;%center of each bin (0, dh/2, 3dh/2, etc)
H=[0:nh]*dh;%edge of each bin (0, dh, 2dh, etc)
x=[0:nx]*dx;
t=[0:nt];

[xx,HH,tt]=meshgrid(x,H,t);
[xx,hh,tt]=meshgrid(x,h,t);

%the forcing f and and cumulative distribution G are at bin edges, as H is
ff=1*cos(tt/2000*12*pi);%seasonal cycle
% ff=1*cos(tt/nt*12*pi);%seasonal cycle
% ff=1*cos(tt/nt*4*pi);%seasonal cycle
% ff=ff.*( (1-hh/nh/dh).^2 );%growth decreases with h
%
%make the forcing assymetric to have net ice growth
% ii=find(ff<0); ff(ii)=ff(ii)/1.2;

%advective velocity
if choiceAdv==1;
    %flow from left to right, converging gradually
    u=0.05*( (1-xx/nx/dx) );%m/s
    % u=0.001+0*xx;%m/s
    % u=0.05*( (1-xx/nx/dx).*2 );%m/s
elseif choiceAdv==2;
    %random velocity field, smoothed
    u=0.2*randn(nx+1,nt+1);
    u=runmean(u,5,1);
    u=runmean(u,25,2);
    u=repmat(u,[1 1 nh+1]);
    u=permute(u,[3 1 2]);
elseif choiceAdv==0;
    u=0*xx;
end;

%diffusivity (if advection is included, do diffusion too)
if choiceAdv==0|choiceDiff==0;
    kappa=0*xx;
else;
    kappa=0.02+0*xx;
end;

%initial condition:
if choiceIC==0;%all open water
    c=ones(nh+1,nx+1,nt+1); c(2:end,:,:)=0; %ice concentration
    v=zeros(nh+1,nx+1,nt+1); %ice volume
elseif choiceIC==1;%half open water, half ice bloc
    c=zeros(nh+1,nx+1,nt+1); c(1,:,:)=1/2; c(60*nh/120,:,:)=1/2; v=c.*hh;
%     c=zeros(nh+1,nx+1,nt+1); c(1,:,:)=1/2; c(90*nh/120,:,:)=1/2; v=c.*hh;
elseif choiceIC==2;%homogeneous ice distribution over 60 points
    c=zeros(nh+1,nx+1,nt+1); nn=60; c(2:1+nn,:,:)=1/nn; v=c.*hh;
elseif choiceIC==4;%fully homogeneous ditrib (original adv only pb)
    c=zeros(nh+1,nx+1,nt+1); nn=nh; c(2:1+nn,:,:)=1/nn; v=c.*hh;
end;

%switch off thermo forcing:
if ~doThermo; ff=ff*0; end;

dc_adv=zeros(nh+1,nx+1,nt+1);
dv_adv=zeros(nh+1,nx+1,nt+1);
dc_leads=zeros(nh+1,nx+1,nt+1);

dc_ridges=zeros(nh+1,nx+1,nt+1);
dv_ridges=zeros(nh+1,nx+1,nt+1);

dc_dif=zeros(nh+1,nx+1,nt+1);
dv_dif=zeros(nh+1,nx+1,nt+1);

dGthermo=zeros(nh+1,nx+1,nt+1);

if ~isnan(freqPlot);
    myfig=figure; set(gcf,'Units','Normalized','Position',posPlot);
end;
for it=1:nt;
    
    %place holder for advection/leads/ridges/diff
    %and the potential recasting that ensues
    c(:,:,it+1)=c(:,:,it); cstar=c(:,:,it);
    v(:,:,it+1)=v(:,:,it); vstar=v(:,:,it);
    
    if doThorndike;
        v(:,:,it+1)=c(:,:,it+1).*hh(:,:,it+1);
        vstar=v(:,:,it+1);
    end;
    
    %====== PART 1 ====== diffusion =====
    tmp1=(kappa(:,1:nx,it)+kappa(:,2:nx+1,it))/2;
    tmp2=-(c(:,2:nx+1,it)-c(:,1:nx,it))/dx;
    dc_dif(:,1,it)=-tmp1(:,1).*tmp2(:,1)*dt/dx;
    dc_dif(:,2:nx,it)=-diff(tmp1.*tmp2,1,2)*dt/dx;
    dc_dif(:,nx+1,it)=tmp1(:,nx).*tmp2(:,nx)*dt/dx;
    cstar=cstar+dc_dif(:,:,it);%step forward
    %
    tmp1=(kappa(:,1:nx,it)+kappa(:,2:nx+1,it))/2;
    tmp2=-(v(:,2:nx+1,it)-v(:,1:nx,it))/dx;
    dv_dif(:,1,it)=-tmp1(:,1).*tmp2(:,1)*dt/dx;
    dv_dif(:,2:nx,it)=-diff(tmp1.*tmp2,1,2)*dt/dx;
    dv_dif(:,nx+1,it)=tmp1(:,nx).*tmp2(:,nx)*dt/dx;
    vstar=vstar+dv_dif(:,:,it);%step forward
    
    %====== PART 2 ====== advection ======
    tmp1=u(:,1:nx,it);
    tmp2a=c(:,1:nx,it);
    tmp2b=c(:,2:nx+1,it);
    %upwind scheme
    tmp2=tmp2a; ii=find(tmp1<0); tmp2(ii)=tmp2b(ii);
    dc_adv(:,1,it)=-tmp1(:,1).*tmp2(:,1)*dt/dx;
    dc_adv(:,2:nx,it)=-diff(tmp1.*tmp2,1,2)*dt/dx;
    dc_adv(:,nx+1,it)=tmp1(:,nx).*tmp2(:,nx)*dt/dx;
    cstar=cstar+dc_adv(:,:,it);%step forward
    %
    tmp1=u(:,1:nx,it);
    tmp2a=v(:,1:nx,it);
    tmp2b=v(:,2:nx+1,it);
    %upwind scheme
    tmp2=tmp2a; ii=find(tmp1<0); tmp2(ii)=tmp2b(ii);
    dv_adv(:,1,it)=-tmp1(:,1).*tmp2(:,1)*dt/dx;
    dv_adv(:,2:nx,it)=-diff(tmp1.*tmp2,1,2)*dt/dx;
    dv_adv(:,nx+1,it)=tmp1(:,nx).*tmp2(:,nx)*dt/dx;
    vstar=vstar+dv_adv(:,:,it);%step forward
    
    %====== PART 3 ===== leads =====
    tmp1=u(:,1:nx,it);
    tmp3=[-tmp1(1,1)/dx -diff(tmp1(1,:),1,2)/dx tmp1(1,nx)/dx]*dt;%divergence
    tmp4=min([tmp3;cstar(1,:)]);%the part for leads opening/closing
    dc_leads(1,:,it)=-tmp4(1,:);
    cstar=cstar+dc_leads(:,:,it);%step forward
    %====== PART 3.5 ====== ridges =====
    [ctmp,vtmp]=ice_ridge_gael(cstar,vstar);
    dc_ridges(:,:,it)=ctmp-cstar;
    dv_ridges(:,:,it)=vtmp-vstar;
    cstar=cstar+dc_ridges(:,:,it);%step forward
    vstar=vstar+dv_ridges(:,:,it);%step forward
    
    c(:,:,it+1)=c(:,:,it+1)+dc_dif(:,:,it);
    v(:,:,it+1)=v(:,:,it+1)+dv_dif(:,:,it);
    c(:,:,it+1)=c(:,:,it+1)+dc_adv(:,:,it);
    v(:,:,it+1)=v(:,:,it+1)+dv_adv(:,:,it);
    c(:,:,it+1)=c(:,:,it+1)+dc_leads(:,:,it);
    %dv_leads=0
    c(:,:,it+1)=c(:,:,it+1)+dc_ridges(:,:,it);
    v(:,:,it+1)=v(:,:,it+1)+dv_ridges(:,:,it);
    
    if doThorndike; cPreTh=c(:,:,it+1); end;
    
    %====== PART 0 ===== THERMO =====
    %1) time step ice volume due to thermo dynamics
    
    %1.1) compute volume increment
    %1.1.1) raw growth rate
    dv=ff(:,:,it)*dt.*c(:,:,it+1);
    %1.1.2) preclude excessive melting
    dv2=dv; dv2(2:end,:)=max(-v(2:end,:,it+1),dv(2:end,:));
    %1.1.3) move excess ice melting to open water heat
    dv2(1,:)=dv2(1,:)+sum(dv(2:end,:)-dv2(2:end,:),1);
    %1.1.35) check conservation then overwrite dv
    if abs(sum(dv(:))-sum(dv2(:)))>1e-6;
        error('loosing conservation');
    else;
        dv=dv2;
    end;
    %1.1.3) move excess open water warming to ice melt
    tmp1=sum(v(2:end,:,it+1)+dv(2:end,:),1);%is >=0 (potential melting left)
    tmp2=min(dv(1,:),0);%now <=0 (required warming or melting)
    tmp3=max(tmp1+tmp2,0);%residual ice after using open ocean warming to melt
    dvMelt=tmp3-tmp1;%part that can be used to melt ice
    dvWarm=dv(1,:)-dvMelt;%part that can be used to warm open ocean
    tmp1=v(2:end,:,it+1)+dv(2:end,:);%distribution of available ice
    tmp2=ones(nh,1)*(dvMelt./sum(tmp1,1));%fraction of melting
    %now actually move warming to melting, if needed, and in proportion
    ii=find(sum(tmp1,1)>0);
    dv2=dv; dv2(1,ii)=dvWarm(ii);
    dv2(2:end,ii)=dv(2:end,ii)+tmp1(:,ii).*tmp2(:,ii);
    %1.1.4) check conservation then overwrite dv
    if abs(sum(dv(:))-sum(dv2(:)))>1e-6;
        error('loosing conservation');
    else;%to more closely reproduce ice_thorndike results : set doHeatCons to 0
        dv=dv+doHeatCons*(dv2-dv);
    end;
    %1.1.5) apply dv(1) as warming/cooling of open water before melt/freeze
    %   dv<0 (warming) => add to v (<0)
    %   dv>0 (cooling) => add to v until v=0
    dv2=doHeatCons*min(-v(1,:,it+1),dv(1,:));
    v(1,:,it+1)=v(1,:,it+1)+dv2;
    dv(1,:)=dv(1,:)-dv2;
    
    %1.2) apply increment, except for open water
    v(2:end,:,it+1)=v(2:end,:,it+1)+dv(2:end,:);
    %1.2.1) make sure that we dont leave conc without vol
    tmp1=1*(v(2:end,:,it+1)==0); %mask empty classes
    tmp1=tmp1.*c(2:end,:,it+1); %collect conc
    dc=[sum(tmp1,1); -tmp1]; %move conc to open water
    c(:,:,it+1)=c(:,:,it+1)+dc;
    %1.3) if open water 'grows' : expand ice cover at local h (or median)
    %   - add dv(1) into vv(2)
    %   - convert open water left to h ice (dc(1)=-dc(2))
    ii=find(dv(1,:)>0);
    if ~isempty(ii);
        1;
    end;
    %use small h:
    hmin=HH(1,:,it)+(HH(2,:,it)-HH(1,:,it))/10;
    %use mean h of cat 1
    %     hmin=hh(2,:,it);
    %use current h:
    %     hmin=v(2,:,it+1)./c(2,:,it+1);
    %     hmin2=hh(2,:,it);
    %     hmin(find(isnan(hmin)))=hmin2(find(isnan(hmin)));
    cmax=c(1,:,it+1);
    dc=min(cmax,dv(1,ii)/hmin(ii));
    %
    v(2,ii,it+1)=v(2,ii,it+1)+dv(1,ii);
    c(1,ii,it+1)=c(1,ii,it+1)-dc(ii);
    c(2,ii,it+1)=c(2,ii,it+1)+dc(ii);
    %1.4) if thin ice melts : shrink ice cover at h
    %   (this also mean : lateral ice melt for thin ice)
    %   - dv(2) has already been added to v(2)
    %   - compute dc(2) to keep h constant, move dc(2) to open water
    ii=find(dv(2,:)<0);
    if ~isempty(ii);
        1;
    end;
    cstar=c(2,ii,it+1).*v(2,ii,it+1)./(v(2,ii,it+1)-dv(2,ii));
    dc=cstar-c(2,ii,it+1);
    c(1,ii,it+1)=c(1,ii,it+1)-dc;
    c(2,ii,it+1)=c(2,ii,it+1)+dc;
    
    [V,C]=ice_regularize(v(:,:,it+1),c(:,:,it+1),it);
    v(:,:,it+1)=V;
    c(:,:,it+1)=C;
    
    if doLipscomb;
        [V,C,H]=ice_remap(v(:,:,it+1),c(:,:,it+1),HH(:,:,it+1));
    elseif doBitz;
        [V,C,H]=ice_recast(v(:,:,it+1),c(:,:,it+1),HH(:,:,it+1));
    else;
        %none -> thorndike overwrite the above c and v thermo increments
    end;
    v(:,:,it+1)=V;
    c(:,:,it+1)=C;
    
    if doThorndike;
        %0) switch to cumulative distribution
        G=cumsum(cPreTh);
        %1) a la Thorndike : step H to Hstar (new characteristics), then interpolate G back to H
        % 1.1) compute hstar:
        Hstar=HH(:,:,it)+ff(:,:,it)*dt;
        % 1.2) interpolate G (now at hstar) back to h:
        for ix=1:nx+1;
            tmp1=interp1(Hstar(:,ix),G(:,ix),HH(2:nh,ix,it));%too diffusive
%             tmp1=interp1(Hstar(:,ix),G(:,ix),HH(2:nh,ix,it),'cubic');%ok
            dGthermo(2:nh,ix,it)=tmp1-G(2:nh,ix);
        end;
        %2) open water growth (dirac at h=0);
        % >> new ice cover, of dh/2 thickness, from open water
        % 2.1) move dirac at 0 : d[0] >>> d[fdt] >>> new volume G(1)fdt
        % 2.2) reshape to d[dh/2] thickness >>> add 2/dh*G(1)fdt to first box g
        % 2.2) equivalently : reshape to homogeneous distribution
        % 2.3) that is : remove it from G(1)
        mskGrow=1*(ff(1,:,it)>0);
        dGthermo(1,:,it)=dGthermo(1,:,it)-2/dh*G(1,:).*ff(1,:,it)*dt.*mskGrow;
        %3) open water lead opening (ice melt)
        % >> assume homogeneous distribution g(0+)=g(dh/2)=g(h-)=dG/dh
        mskMelt=1*(ff(1,:,it)<=0);
        dGdh=(G(2,:)-G(1,:))/dh;
        dGthermo(1,:,it)=dGthermo(1,:,it)-dGdh.*ff(1,:,it)*dt.*mskMelt;
        %4) step G forward and then c & v
        Gstar=G+dGthermo(:,:,it);
        cstar=Gstar; cstar(2:end,:)=Gstar(2:end,:)-Gstar(1:end-1,:);
        c(:,:,it+1)=cstar;
        v(:,:,it+1)=cstar.*hh(:,:,it+1);
    end;
    
    if ~isempty(find(isnan(v(:,:,it+1))|isnan(c(:,:,it+1))));
        error('model blew up');
    end;
    
    %display:
    if ~isnan(freqPlot)&mod(it,freqPlot)==0;
        figure(myfig);
        subplot(2,1,1); imagesc(c(:,:,it+1)); caxis([-0.01 0.1]); title(num2str(it+1)); axis([0.5 51.5 0 nh]); axis xy;
        subplot(2,1,2); imagesc(cumsum(c(:,:,it+1),1)); caxis([-0.1 1.1]); title(num2str(it+1)); axis([0.5 51.5 0 nh]); axis xy;
        if isfinite(pausePlot); pause(pausePlot); else; pause; end;
    end;
end;
%display:
if ~isnan(freqPlot);
    figure(myfig);
    subplot(2,1,1); 
    imagesc(c(:,:,it+1)); caxis([-0.01 0.1]); 
    axis([0.5 51.5 0 nh]); axis xy; 
    ylabel('ice thickness'); xlabel('position');
    title(['final distribution (at t='  num2str(it*dt) ')']);     
    subplot(2,1,2); 
    imagesc(cumsum(c(:,:,it+1),1)); caxis([-0.1 1.1]); 
    axis([0.5 51.5 0 nh]); axis xy;
    ylabel('ice thickness'); xlabel('position');
    title(['final concentration (at t='  num2str(it*dt) ')']);     
end;

%store result:
eval(['c_' suff '=c;']); eval(['v_' suff '=v;']);
% figure; plot(c_thor(:,1,end),'g'); hold on; plot(c_bitz(:,1,end),'r'); plot(c_lips(:,1,end),'b');
% figure; plot(c_thor(:,end,end),'g'); hold on; plot(c_bitz(:,end,end),'r'); plot(c_lips(:,end,end),'b');

%plot ice totals and ice distrib in select points
%================================================
figure; set(gcf,'Units','Normalized','Position',posPlot);
%
subplot(2,1,1);
plot(squeeze(mean(sum(c(2:end,:,:),1),2)),'r'); hold on;
plot(squeeze(mean(sum(0.1*v(2:end,:,:),1),2))); 
legend('concentation','volume/10'); 
aa=axis; aa(2)=nt+1; axis(aa); grid on;
title('total ice'); xlabel('time');
%
subplot(2,1,2); 
G=cumsum(c,1); 
cc=[0.01 0.2:0.2:0.8 0.99];
% cc=[0.1:0.2:0.9];
contour(squeeze(G(:,1,:)),cc,'b'); hold on;
% contour(squeeze(G(:,round(nx*1/4),:)),cc,'b'); hold on;
%
% contour(squeeze(G(:,round(nx*1/3),:)),cc,'g');
% contour(squeeze(G(:,round(nx*2/3),:)),cc,'m');
contour(squeeze(G(:,round(nx*1/2),:)),cc,'g');
%
contour(squeeze(G(:,nx+1,:)),cc,'r');
% contour(squeeze(G(:,round(nx*3/4),:)),cc,'r');
%
legend('left boundary','interior','right boundary','Orientation','horizontal'); 
title('ice distribution (20% between curves)');
ylabel('ice thickness'); xlabel('time');


if displayChecks;
    %display volume time series: conservation check
    %===========================
    tmp1=[zeros(nx+1,1) cumsum(squeeze(ff(1,:,1:it-1)*dt),2)];
    tmp1=tmp1+squeeze(sum(v(:,:,1),1))'*ones(1,nt);
    tmp2=squeeze(sum(v(:,:,1:it),1));
    tmp3=squeeze(v(1,:,1:it));
    %
    tmp1=mean(tmp1,1); tmp2=mean(tmp2,1); tmp3=mean(tmp3,1);
    %
    figure; set(gcf,'Units','Normalized','Position',posPlot);
    subplot(2,1,1); plot(tmp1); hold on; plot(tmp3,'ro'); plot(tmp2-tmp3,'gx');
    legend('forcing','ocn heat','ice vol');
    subplot(2,1,2); plot(tmp2-tmp1,'k'); title('res');
    
    %display concentration time series: conservation check
    %==================================
    tmp1=squeeze(c(1,:,2:nt+1));
    tmp2=squeeze(sum(c(2:end,:,2:nt+1),1));
    %
    tmp1=mean(tmp1,1); tmp2=mean(tmp2,1);
    %
    figure; set(gcf,'Units','Normalized','Position',posPlot);
    subplot(2,1,1); plot(tmp1,'.-'); hold on; plot(tmp2,'rx-');
    axis([0 nt+1 -0.1 1.1]); legend('open ocean conc','ice conc');
    subplot(2,1,2); plot(1-tmp1-tmp2,'k'); title('res');
end;

