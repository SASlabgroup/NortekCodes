% OLD version (neglects cross-spectra) 

clear all 

load('TTM_APLUWxIMU_Jun2012.mat'),
load('TTM_APLUWvector_Jun2012_processed_300s.mat')

m = min( [ size(xIMU.Saxf,2) size(ADV.Surawf,2) ]);


%% change acc units from g to m/s^2
xIMU.Saxf = xIMU.Saxf * 9.8^2;
xIMU.Sayf = xIMU.Sayf * 9.8^2;
xIMU.Sazf = xIMU.Sazf * 9.8^2;

%% bin spectra by mean velocity
du = 0.2;
ubar = [0.1:du:max(ADV.ubar)];

for ui = 1:length(ubar),  
    uinds = find( abs(ADV.ubar(1:m)-ubar(ui)) < du./2 );
    ADV.binnedSurawf(:,ui) = nanmean( ADV.Surawf(:,uinds), 2 );
    xIMU.binnedSaxf(:,ui) = nanmean( xIMU.Saxf(:, uinds), 2 );
end

%% binned spectra motion correction, with frequency filter

twopi = 2*pi;
clear filter
% fcutoff = 7e-2;
% 
% filter( ADV.f < fcutoff ) = 0; %filter(  ADV.f >= 7e-2 ) = 5e-3;
% %filter = logspace(1e-11,1e-10,length(ADV.f));
% filter( find( ADV.f >= fcutoff) ) = linspace(0,20,length( find( ADV.f >= fcutoff) ));
% %filter = 1e-3*filter' * ( nansum(ADV.binnedSurawf) ./ nansum(xIMU.binnedSaxf)  ).^-1;
% filter = (filter' * ones(1,length(ubar)));

h = 8; % mooring length at instrument
fcutoff = ubar./ h;
scale = 14 * ones(1,length(ubar));
%scale = 10*nansum(ADV.binnedSurawf);
%scale = 2*nansum(xIMU.binnedSaxf);
for ui = 1:length(ubar),
    filter( ADV.f < fcutoff(ui) , ui ) = 0;  % no motion correction a frequencies below lowest mooring mode
    filter( ADV.f >= fcutoff(ui) , ui ) = linspace(0,scale(ui), sum(ADV.f >= fcutoff(ui)) );
end

ADV.binnedtrueSurawf = ADV.binnedSurawf - ( filter .* xIMU.binnedSaxf .* ( ( twopi .* ADV.f) * ones(1,length(ubar)) ).^(-2) ); 


%% plot binned results


ax = [1e-2 10 1e-5 2];

figure(1),clf
cmap = colormap;
for i=1:length(ubar),
    cmapindex = round(ubar(i) ./ max(ubar) * 64);    if cmapindex==0, cmapindex=1;, else end
    loglog(ADV.f,ADV.binnedSurawf(:,i),'color',cmap(cmapindex,:),'linewidth',2),
    hold on
end
%loglog(ADV.f,nanmean(ADV.Surawf(:,good),2),'k-','linewidth',3),
axis(ax)
set(gca,'Fontsize',18,'fontweight','demi')
title('Raw ADV')
xlabel('frequency [Hz]'),ylabel('S_{uu} [m^2 s^{-2}/HZ]'),
plot([2e-1 2e0],1e-2*[2e-1 2e0].^(-5/3),'k--','linewidth',1.5),
text(4e-1,1e-1,'f^{-5/3}','fontsize',16,'fontweight','demi')
cb = colorbar; axes(cb), set(gca,'YTickLabel',['0.3'; '0.6'; '0.9'; '1.2'; '1.5'; '1.8';])
set(gca,'Fontsize',18,'fontweight','demi')
ylabel('<u> [m/s]')

figure(2),clf
cmap = colormap;
for i=1:length(ubar),
    cmapindex = round(ubar(i) ./ max(ubar) * 64);    if cmapindex==0, cmapindex=1;, else end
    loglog(ADV.f,xIMU.binnedSaxf(:,i),'color',cmap(cmapindex,:),'linewidth',2),
    hold on
end
%loglog(ADV.f,nanmean(xIMU.Saxf(:,good),2),'k-','linewidth',3),
axis(ax)%[ax(1:2) 1e-7 1e-1])
set(gca,'Fontsize',18,'fontweight','demi')
title('Acceleration')
xlabel('frequency [Hz]'),ylabel('S_{aa} [m^2 s^{-4}/HZ]')
cb = colorbar; axes(cb), set(gca,'YTickLabel',['0.3'; '0.6'; '0.9'; '1.2'; '1.5'; '1.8';])
set(gca,'Fontsize',18,'fontweight','demi')
ylabel('<u> [m/s]')

    
figure(3), clf
for i=1:length(ubar),
    cmapindex = round(ubar(i) ./ max(ubar) * 64);    if cmapindex==0, cmapindex=1;, else end
    loglog(ADV.f,ADV.binnedtrueSurawf(:,i),'color',cmap(cmapindex,:),'linewidth',2),
    hold on
end
axis(ax)
%loglog(ADV.f,nanmean(ADV.trueSuuf(:,good),2),'k-','linewidth',3),
%loglog(ADV.f,nanmean(ADV.Surawf(:,good),2),'k-','linewidth',3),
set(gca,'Fontsize',18,'fontweight','demi')
title('True TKE')
xlabel('frequency [Hz]'),ylabel('TKE [m^2 s^{-2}/HZ]')
plot([2e-1 2e0],1e-2*[2e-1 2e0].^(-5/3),'k--','linewidth',1.5),
text(4e-1,1e-1,'f^{-5/3}','fontsize',16,'fontweight','demi')
cb = colorbar; axes(cb), set(gca,'YTickLabel',['0.3'; '0.6'; '0.9'; '1.2'; '1.5'; '1.8';])
set(gca,'Fontsize',18,'fontweight','demi')
ylabel('<u> [m/s]')


% %% all spectra motion correction (x dir)
% 
% filter( ADV.f < 1e-1 ) = 0;
% filter(  ADV.f >= 1e-1 ) = 1;
% 
% ADV.trueSuuf(:,1:m) = ADV.Surawf(:,1:m) - ( (filter'*ones(1,m)) .* xIMU.Saxf(:,1:m) .* ( ( 2 .* pi .* ADV.f) * ones(1,m) ).^(-2) ); 
% %ADV.trueSuuf(:,1:m) = ADV.Surawf(:,1:m) - ( (filter'*ones(1,m)) .* xIMU.Saxf(:,1:m) ) ; 
% ADV.trueSuuf( ADV.trueSuuf <= 0 ) = NaN;
% 
% %% plots (x dir)
% 
% ax = [1e-2 10 1e-5 2];
% cutin = 0;
% good = find(ADV.ubar(1:m) > cutin);
% 
% figure(1),clf
% cmap = colormap;
% for i=good,
%     cmapindex = floor(ADV.ubar(i) ./ max(ADV.ubar) * 64);
%     loglog(ADV.f,ADV.Surawf(:,i),'color',cmap(cmapindex,:)),
%     hold on
% end
% loglog(ADV.f,nanmean(ADV.Surawf(:,good),2),'k-','linewidth',3),
% axis(ax)
% set(gca,'Fontsize',18,'fontweight','demi')
% title('Raw TKE')
% 
% figure(2),clf
% cmap = colormap;
% for i=good,
%     cmapindex = floor(ADV.ubar(i) ./ max(ADV.ubar) * 64);
%     loglog(ADV.f,xIMU.Saxf(:,i),'color',cmap(cmapindex,:)),
%     hold on
% end
% loglog(ADV.f,nanmean(xIMU.Saxf(:,good),2),'k-','linewidth',3),
% axis(ax)%[ax(1:2) 1e-7 1e-1])
% set(gca,'Fontsize',18,'fontweight','demi')
% title('Acceleration')
% 
%     
% figure(3), clf
% for i=good,
%     cmapindex = floor(ADV.ubar(i) ./ max(ADV.ubar) * 64);
%     loglog(ADV.f,ADV.trueSuuf(:,i),'color',cmap(cmapindex,:)),
%     hold on
% end
% axis(ax)
% loglog(ADV.f,nanmean(ADV.trueSuuf(:,good),2),'k-','linewidth',3),
% %loglog(ADV.f,nanmean(ADV.Surawf(:,good),2),'k-','linewidth',3),
% set(gca,'Fontsize',18,'fontweight','demi')
% title('True TKE')
% 
% 
% %% spectra motion correction (y dir)
% 
% m = min( [ size(xIMU.Sayf,2) size(ADV.Svrawf,2) ]);
% 
% filter( ADV.f < 1e-1 ) = 0;
% filter(  ADV.f >= 1e-1 ) = 1;
% twopi = 2*pi;
% 
% ADV.trueSvvf(:,1:m) = ADV.Svrawf(:,1:m) - ( (filter'*ones(1,m)) .* xIMU.Sayf(:,1:m) .* ( ( twopi .* ADV.f) * ones(1,m) ).^(-2) ); 
% 
% 
% %% plots (y dir)
% 
% ax = [1e-2 10 1e-5 2];
% cutin = 0;
% good = find(ADV.ubar(1:m) > cutin);
% 
% figure(1),clf
% cmap = colormap;
% for i=good,
%     cmapindex = floor(ADV.ubar(i) ./ max(ADV.ubar) * 64);
%     loglog(ADV.f,ADV.Svrawf(:,i),'color',cmap(cmapindex,:)),
%     hold on
% end
% loglog(ADV.f,nanmean(ADV.Svrawf(:,good),2),'k-','linewidth',3),
% axis(ax)
% set(gca,'Fontsize',18,'fontweight','demi')
% title('Raw TKE')
% 
% figure(2),clf
% cmap = colormap;
% for i=good,
%     cmapindex = floor(ADV.ubar(i) ./ max(ADV.ubar) * 64);
%     loglog(ADV.f,xIMU.Sayf(:,i),'color',cmap(cmapindex,:)),
%     hold on
% end
% loglog(ADV.f,nanmean(xIMU.Sayf(:,good),2),'k-','linewidth',3),
% axis(ax)
% set(gca,'Fontsize',18,'fontweight','demi')
% title('Acceleration')
% 
%     
% figure(3), clf
% for i=good,
%     cmapindex = floor(ADV.ubar(i) ./ max(ADV.ubar) * 64);
%     loglog(ADV.f,ADV.trueSvvf(:,i),'color',cmap(cmapindex,:)),
%     hold on
% end
% axis(ax)
% loglog(ADV.f,nanmean(ADV.trueSvvf(:,good),2),'k-','linewidth',3),
% set(gca,'Fontsize',18,'fontweight','demi')
% title('True TKE')
