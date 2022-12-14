% motion correct ADV on TTM

clear all
close all
clc

load('/Volumes/Data/ADMIRALTY/ADMIRALTY_Jun2012/TTM_xIMU_Jun2012/TTM_APLUWxIMU_Jun2012_burst.mat')
load('/Volumes/Data/ADMIRALTY/ADMIRALTY_Jun2012/TTM_Vectors/TTM_APLUWvector_Jun2012_despikedbursts.mat')

m = min( [ size(burstaccx,2) size(burstu,2) ] );

%% fix acc units
burstaccx = burstaccx * 9.8;
burstaccy = burstaccy * 9.8;
burstaccz = burstaccz * 9.8;

%% make raw crossspectra

for i=1:m,
    
    n = length(burstaccx);
    w = window(@hamming,n);
    
    axfft(:,i) = fft( detrend( burstaccx(:,i) ) .* w );
    ayfft(:,i) = fft( detrend( burstaccy(:,i) ) .* w);
    azfft(:,i) = fft( detrend( burstaccz(:,i) ) .* w);
    
    ufft(:,i) = fft( detrend( burstu(:,i) ) .* w);
    vfft(:,i) = fft( detrend( burstv(:,i) ) .* w);
    wfft(:,i) = fft( detrend( burstw(:,i) ) .* w);

end

% make normalized spectra

Saxax = axfft .* conj(axfft) ./ n;
Suu = ufft .* conj(ufft) ./ n;
Xuax = ufft .* conj(axfft) ./ n;



%% ensembles by mean speed
du = 0.2;
ubar = [0:du:max(abs(burstu(:)))];

 clear binned*
for ui = 1:length(ubar), 
    uinds = find( abs( nanmean(abs(burstu(:,1:m))) - ubar(ui) ) < du./2 );
    binnedSaxax(:,ui) = nanmean( Saxax(:,uinds), 2 );
    binnedSuu(:,ui) = nanmean( Suu(:,uinds), 2 );
    binnedXuax(:,ui) = nanmean( Xuax(:,uinds), 2 );
end

%% remove neg frequencies 
fs = 32;
f = [(fs/n):(fs/n):fs/2];

binnedSaxax([(n/2+1):n],:) = [];
binnedSuu([(n/2+1):n],:) = [];
binnedXuax([(n/2+1):n],:) = [];


%% merge frequencies
merge = 5;

for i = merge:merge:(n/2) 
	finalSaxax(i/merge,:) = nanmean( binnedSaxax((i-merge+1):i , : ) );
    finalSuu(i/merge,:) = nanmean( binnedSuu((i-merge+1):i , : ) );
    finalXuax(i/merge,:) = nanmean( binnedXuax((i-merge+1):i , : ) );
end

finalf = f( [(merge/2-.5):merge:(n/2)] );

figure(1),clf
loglog(finalf,finalSuu)

figure(2), clf
loglog(finalf,finalSaxax)

figure(3), clf
loglog(finalf,abs(finalXuax))

%% motion correction
tp = 2 * pi;

for ui = 1:length(ubar),
    trueSuu(:,ui) = finalSuu(:,ui)    -   2 * real( finalXuax(:,ui) ) .* (tp * finalf(:,ui))^(-1)   +    finalSaxax(:,ui).* (tp * finalf(:,ui)).^(-2);
end

figure(4), clf
loglog(finalf,trueSuu)

