
function [ freq, SSE, k, PP, PU, PV, UU, VV, UV, dir1, spread1, dir2, spread2, dir2_swell, spread2_swell, Hsig_swell, Hsig_ig, centroid_swell, centroid_ig, in, out, PsqP, PsqU, Depth ]  =  PUVspectra( file , window_length, merge, plot );

%
% Matlab function to calculate and plot spectra from a sinlge PUV sensor
%    using time series stored in a "mmddhhhh.ID.mat" file
%    (where all measurements are in one file... see RAW-processing codes)
%
% Input is a file string 'mmddhhhh.ID.mat', the window size (in seconds), 
%    an integer with the number of freq bands to merge, 
%    and a toggle to plotting option (1=ON, 0=OFF)
% 
%      >>  PUVspectra (  file,  window,  merge,  plot? (1 or 0)   );
%
%    ** note that window size should be 2^N for FFT efficiency **
%          and the number of bands to merge should be odd
%
% Output is freq (Hz) and power spectral density of sea surface elevation (cm^2/Hz).
%    ALSO wavenumber from linear theory, autospectra of bottom pressure,
%    cross-spectra of pressure and veloctiy PU and PV, 
%    autospectra of U and V velocities (x and y), cross UV,
%    direction1, spread1, direction2, spread2,
%    weighted mean swell direction2, weighted mean swell spread2, 
%    significant swell height, significant infragravity height, 
%    centroid swell frequency, centriod infragravity frequency,
%    energy fluxes 'in' and 'out' of the beach (x component),
%    cross spectra of P^2,P and of P^2,U  (sum to get bispectra),
%    and mean water depth (pres) for the record
%
% function also makes 2 PLOTS of pressure spectra, velocity spectra 
%    (including coherence & phase), wave direction, equivalent pressure,
%    and energy flux.  (Velocities are converted to "equivalent pressure" 
%    to graphically check quality of data)
%
% Use function "crossPUVspectra.m" to compare data from two different PUVs
%
% J. Thomson, September 2003, 
%               *** modified to included higher order spectra (Aug 2005) **

warning off MATLAB:divideByZero; close all;

% use these lines to run this as script instead of a function:
	% file = input('enter filename (including path): ', 's');
	% merge = input('enter number of f-bands to merge: ');
	% window_length = input('enter window size (seconds): ');

% LOAD DATA -------------------------------------------------------------------
% binary .mat (with appropriate dir path) with entire time series
% convention is to keep variables lowercase in time domain,
% and then use capitals in frequency domain
load( ['/home/jthomson/NCEX/DATA/Triton/' file(1:4) '/' file ] );
pts = length(p);                 % record length in data points
secs = pts / sample_rate;        % record length in seconds
pcleaned, vcleaned               % print quality strings to screen
clear k l m;                     % clear the beam correlation variables 
Depth = mean(p) + 55;		 % depth for the entire record (triton p-gage dof~60cm)
%------------------------------------------------------------------------------



% WINDOW THE DATA (use 75 percent overlap)--------------------------------------
% WINDOW LENGTH SHOULD BE 2^n FOR FFT EFFICIENCY
w = sample_rate * window_length;    % window length in data points
windows = 4*(pts/w -1)+1;   % number of windows, the 4 comes from a 75% overlap
% loop to create a matrix of time series, where COLUMN = WINDOW 
for q=1:windows, 
	p_window(:,q) = p(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
	u_window(:,q) = u(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
	v_window(:,q) = v(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
end
%-------------------------------------------------------------------------------



% WATER DEPTH (use later for depth correction)----------------------------------
% offsets from central transducer for SonTek Tritons 
% (see triton.pdf, which is in inches!!!)
doff = 78;			  % distance of cental tranducer of bottom
				  % (measured by diver during deployment)
doffu = doff + 10;                % distance (cm) of sample volume off bottom 
doffp = doff - 23.25;             % distance (cm) of p gage of bottom
depth = mean(p_window) + doffp;   % avg p (cm) depth for that column (window)
%-------------------------------------------------------------------------------



% DETREND THE WINDOWED DATA-----------------------------------------------------
% remove the mean of each window
p_window_nomean = p_window - ( ones(w,1) * mean(p_window)  ) ; 
u_window_nomean = u_window - ( ones(w,1) * mean(u_window)  ) ; 
v_window_nomean = v_window - ( ones(w,1) * mean(v_window)  ) ; 
% loop to remove quadratic trend (attempt to reduce tidal leakage)
secs = [sample_rate^-1:sample_rate^-1:(w/sample_rate)]';
for q=1:windows
quadfitp = polyfit(secs,p_window_nomean(:,q),2);
p_window_detrend(:,q) = p_window_nomean(:,q) - (quadfitp(3) + secs.*quadfitp(2) + (secs.^2).* quadfitp(1) );
quadfitu = polyfit(secs,u_window_nomean(:,q),2);
u_window_detrend(:,q) = u_window_nomean(:,q) - (quadfitu(3) + secs.*quadfitu(2) + (secs.^2).* quadfitu(1) );
quadfitv = polyfit(secs,v_window_nomean(:,q),2);
v_window_detrend(:,q) = v_window_nomean(:,q) - (quadfitv(3) + secs.*quadfitv(2) + (secs.^2).* quadfitv(1) );
end
%------------------------------------------------------------------------------




% TAPER THE DATA (use a Hanning type window)-----------------------------------
% form taper matrix (columns of taper coef)
taper = sin ([1:w] * pi/w )' * ones(1,windows); 
% taper each window
p_window_taper = p_window_detrend .* taper;
u_window_taper = u_window_detrend .* taper;
v_window_taper = v_window_detrend .* taper;
% now find the correction factor (comparing old/new variance)
factp = sqrt( var(p_window_detrend) ./ var(p_window_taper) );
factu = sqrt( var(u_window_detrend) ./ var(u_window_taper) );
factv = sqrt( var(v_window_detrend) ./ var(v_window_taper) );
% and correct for the change in variance
% (mult each window by it's variance ratio factor)
p_window_ready = (ones(w,1)*factp).* p_window_taper;
u_window_ready = (ones(w,1)*factu).* u_window_taper;
v_window_ready = (ones(w,1)*factv).* v_window_taper;
% check & report
if abs(  var(p_window_ready) - var(p_window_detrend)  ) > 0.1,
  disp('******************************')
  disp('Problem preserving variance variance');
  disp('******************************')
  else end
%------------------------------------------------------------------------------



% SPECTRA (FFT)-----------------------------------------------------------------
% calculate Fourier coefs
P_window = fft(p_window_ready);
Psq_window = fft(p_window_ready.^2);  % include for skewness and higher order spectra
U_window = fft(u_window_ready);
V_window = fft(v_window_ready);
% second half of fft is redundant, so throw it out
P_window( (w/2+1):w, : ) = [];
Psq_window( (w/2+1):w, : ) = [];
U_window( (w/2+1):w, : ) = [];
V_window( (w/2+1):w, : ) = [];
% throw out the mean (first coef) and add a zero (to make it the right length)  
P_window(1,:)=[];  Psq_window(1,:)=[]; U_window(1,:)=[]; V_window(1,:)=[]; 
P_window(w/2,:)=0; Psq_window(w/2,:)=0; U_window(w/2,:)=0; V_window(w/2,:)=0; 
% POWER SPECTRA (auto-spectra)
PP_window = ( P_window .* conj(P_window) );
UU_window = ( U_window .* conj(U_window) );
VV_window = ( V_window .* conj(V_window) );
% CROSS-SPECTRA 
PU_window = ( P_window .* conj(U_window) );
PV_window = ( P_window .* conj(V_window) );
UV_window = ( U_window .* conj(V_window) );
% THIRD ORDER SPECTRA (equals sum of bispectra, see Henderson write-up, 6/30/2005)
PsqP_window = ( Psq_window .* conj(P_window) );
PsqU_window = ( Psq_window .* conj(U_window) );
% -----------------------------------------------------------------------------



% MERGE FREQUENCY BANDS -------------------------------------------------------
% raw fft has w/2 frequency bands before merging... merge to improve stastics
% number of bands to merge is an input to function
for i = merge:merge:(w/2) 
	PP_window_merged(i/merge,:) = mean( PP_window((i-merge+1):i , : ) );
	UU_window_merged(i/merge,:) = mean( UU_window((i-merge+1):i , : ) );
	VV_window_merged(i/merge,:) = mean( VV_window((i-merge+1):i , : ) );
	PU_window_merged(i/merge,:) = mean( PU_window((i-merge+1):i , : ) );
	PV_window_merged(i/merge,:) = mean( PV_window((i-merge+1):i , : ) );
	UV_window_merged(i/merge,:) = mean( UV_window((i-merge+1):i , : ) );
	PsqP_window_merged(i/merge,:) = mean( PsqP_window((i-merge+1):i , : ) );
	PsqU_window_merged(i/merge,:) = mean( PsqU_window((i-merge+1):i , : ) );
end
% freq range and bandwidth
n = (w/2) / merge;                         % number of f bands
Nyquist = .5 * sample_rate;                % highest spectral frequency 
bandwidth = Nyquist/n ;                    % freq (Hz) bandwitdh
% find middle of each freq band, ONLY WORKS WHEN MERGING ODD NUMBER OF BANDS!
freq= 1/(window_length) + bandwidth/2 + bandwidth.*[0:(n-1)] ; 
% -----------------------------------------------------------------------------



% ENSEMBLE AVERAGE THE WINDOWS -------------------------------------------------
% take the average of all windows at each freq-band
% and divide by N*sample_rate to get power spectral density
% the two is b/c Matlab's fft output is the symmetric FFT, and we did not use the redundant half (so need to multiply the psd by 2)
PP = mean( PP_window_merged.' ) / (w/2 * sample_rate );    
UU = mean( UU_window_merged.' ) / (w/2 * sample_rate  );
VV = mean( VV_window_merged.' ) / (w/2 * sample_rate  );
PU = mean( PU_window_merged.' ) / (w/2 * sample_rate  ); 
PV = mean( PV_window_merged.' ) / (w/2 * sample_rate  ); 
UV = mean( UV_window_merged.' ) / (w/2 * sample_rate  ); 
PsqP = mean( PsqP_window_merged.' ) / (w/2 * sample_rate  ); 
PsqU = mean( PsqU_window_merged.' ) / (w/2 * sample_rate  ); 
%--------------------------------------------------------------------------



% DEPTH CORRECTION -------------------------------------------------------------
% find correction/conversion rcoefs at each f-band 
% to calc sea surface elevation and convert velocities to pressure units 
g = 981;    % gravity, cm/s^2
for index = 1:n
   f = freq(index);
   fcutoff = 0.25 ;  % cutoff frequency 
		     % (beyond which correction will be too large to trust)
   if f < fcutoff,
   % find k for each f with function wavenumber.m (a Newton-Ralphson iteration)
   k(index) = wavenumber(f,mean(depth));
   correction(index) =  cosh(k(index)*mean(depth)) / cosh(k(index)*doffp);
   convert(index) = (2*pi*f/(g*k(index))) * cosh(k(index)*doffp) / cosh(k(index)*doffu) ;
   % record conversion coefs to check them
   else convert(index)=0; correction(index)=0; k(index)=0; end % beyond fcutoff
end   
% correct pressure for attentuation to get sea-surface elevation
   SSE = PP .* (correction.^2) ; 
% convert velocities into "equivalent pressure" (to compare with PP)
   UUpres = UU .* (convert.^2) ; 
   VVpres = VV .* (convert.^2) ; 
   UVpres = UV .* (convert.^2) ; 
% convert x-spectra as well (but conversion is not squared here)
   PUpres = PU .* convert ; 
   PVpres = PV .* convert ; 
% should check that:   PP^2 = UUpres^2 + VVpres^2,  which is true for linear propagating waves
%-------------------------------------------------------------------------------



% COHERENCE & PHASE, etc -------------------------------------------------------
% Cospectrum & Quadrature:
coPU = real(PU);   quPU = imag(PU);
coPV = real(PV);   quPV = imag(PV);
coUV = real(UV);   quUV = imag(UV);
% Coherence & Phase at each freq-band
% *** note that it's important to calc this AFTER all merging and ensemble avg.
cohPU = sqrt( (coPU.^2 + quPU.^2) ./ (PP.* UU) );
phPU  = 180/pi .* atan2( quPU , coPU );  
cohPV = sqrt((coPV.^2 + quPV.^2)./ (PP.* VV));
phPV  = 180/pi .* atan2( quPV , coPV );  
cohUV = sqrt((coUV.^2 + quUV.^2)./(UU .* VV));
phUV  = 180/pi .* atan2( quUV , coUV );  
% -----------------------------------------------------------------------------



%  WAVE DIRECTION & SPREAD------------------------------------------------------
%  before rotation 0 deg is for waves headed towards positive x (usually EAST)
%  NOTE: alt spread values are can be used to form cos^(2s)(theta/2) spreading function,
%        not as uncertainty in dir, see Donelan et al 1985, p. 44, or Longuet-HIggins et al 1963
%        and, of course, Kuik et al, JPO, 1988
a1 = coPU ./ sqrt( PP .* ( UU + VV ) );
b1 = coPV ./ sqrt( PP .* ( UU + VV ) );
dir1 = rad2deg ( atan2(b1,a1) );          
spread1 = ( sqrt( 2 .* ( 1-sqrt(a1.^2 + b1.^2) ) ) );
% turn direction into 360deg compass heading...
% using direction FROM which waves are coming (like wind)
	%dir1 =  270 + dir1;
	%t = find(dir1 > 360);
	%dir1(t) = dir1(t) - 360;
% other method
a2 = (UU - VV) ./ (UU + VV);
b2 = 2 .* coUV ./ ( UU + VV );
dir2 = rad2deg ( atan2(b2,a2)/2 ); 
spread2 = sqrt(abs( 0.5 - 0.5 .* ( a2.*cos(2.*deg2rad(dir2)) + b2.*cos(2.*deg2rad(dir2)) )  ));
% turn direction into 360deg compass heading...
% using direction FROM which waves are coming (like wind)
	% dir2 = 270 + dir2;
	% t = find(dir2 > 360);
	% dir2(t) = dir2(t) - 360;
% -----------------------------------------------------------------------------



% ESTIMATE ENERGY FLUX (cm^3 s^-1 Hz^-1 ) ------------------------------------------------
% spectral estimator from T.H.C. Herbers, using LFDT group velocity **OR** shallow water:
% robust to partial standing waves, but assumes shore normal propagation
%Cg = 0.5 * (2*pi*freq).* k .* ( 1 + (2*k*mean(depth))./sinh(2*k*mean(depth)) ); 
rho = 0.001028 ;  g = 981;   % density (Kg cm^-3) and gravity (cm s^-2)
Cg = sqrt(g*Depth);   % shallow water version, see Sheremet et al, 2002
const = rho * g * Cg .* ((cosh(k*mean(depth))).^2) ./ ((cosh(k*doffp)).^2); % correct for attenuation 
% energy flux (by freq) in cartesian coordinates (assumes X propagation dominates... valid near beach)
in  = 1/4 .* Cg .* ( abs(PP) + ( ( 2*pi*freq./(g*k) ).^2 ).*abs(UU) + 2.*(2*pi*freq./(g*k)).*real(PU) );
out = 1/4 .* Cg .* ( abs(PP) + ( ( 2*pi*freq./(g*k) ).^2 ).*abs(UU) - 2.*(2*pi*freq./(g*k)).*real(PU) );
%in = const .* (0.5 .* ( abs(PP) + ( ( 2*pi*freq./(981*k) ).^2 ).*( abs(UU) - abs(VV) ) ) + (2*pi*freq./(981*k)).*real(PU) );
%out = const .* (0.5 .* ( abs(PP) + ( ( 2*pi*freq./(981*k) ).^2 ).*( abs(UU) - abs(VV) ) ) - (2*pi*freq./(981*k)).*real(PU) );
% ----------------------------------------------------------------------------



% SPECTRAL WEIGHTED AVERAGES & STATS -----------------------------------------
% find indices of freq bands
ig = find(freq>0.005 & freq<0.03);
swell = find(freq>0.05 & freq<0.2);
% significant wave height (cm):
Hsig_swell = 4 * sqrt( sum( SSE(swell) * bandwidth ) );      
Hsig_ig = 4 * sqrt( sum( SSE(ig) * bandwidth ) );            
% swell direction:
dir1_swell = sum( dir1(swell).*PP(swell) ) / sum(PP(swell)) ;
spread1_swell = sum( spread1(swell).* PP(swell) ) / sum(PP(swell)) ;
dir2_swell = sum( dir2(swell).* PP(swell) ) / sum(PP(swell)) ;
spread2_swell = sum( spread2(swell).* PP(swell) ) / sum(PP(swell)) ;
dir_swell = mean( [ dir1_swell dir2_swell ] );
spread_swell = mean( [ spread1_swell spread2_swell ] );
% centriod frequency
centroid_ig = sum ( freq(ig).* PP(ig) ) / sum ( PP(ig) ) ;
centroid_swell = sum ( freq(swell).* PP(swell) ) / sum ( PP(swell) ) ;
% ----------------------------------------------------------------------------



% DEGREES OF FREEDOM and level of no significant coherence --------------------
% DOF = 2 * (# independent windows) * (# bands merged)
DOF = 2 * pts/w * merge; 
chi2 = chi2pdf([1:100],DOF);  within95 = find(chi2 > 0.05*max(chi2));
low = within95(1)/find(chi2==max(chi2));
high = within95(length(within95))/find(chi2==max(chi2)); 
% 95% significance level for zero coherence
SIG = sqrt(6/DOF);
phPUsig = phPU( find(cohPU > SIG) );
freqPUsig = freq( find(cohPU > SIG) );
phPVsig = phPV( find(cohPV > SIG) ); 
freqPVsig = freq( find(cohPV > SIG) );
% ------------------------------------------------------------------------------



% PLOTTING ---------------------------------------------------------------------
if plot==1,

figure(1), % power spectral density of sea surface, and wave direction
H(1) = subplot(2,1,1);
loglog ( freq, SSE, 'linewidth',4 ), hold on, 
xlabel('frequency (Hz)'), ylabel('spectral energy density (cm^2/Hz)'),
axis( [ 0.002 0.3 5 10^4] ), legend('sea-surface elevation',2),
title([ file ',  bandwidth = ' num2str(bandwidth,2) ',  DOF = ' num2str(DOF)]),
text(0.007,8,['H_{sig}^{IG} = ' num2str(Hsig_ig,2) ' cm']),
text(0.08,8,['H_{sig}^{SW} = ' num2str(Hsig_swell,3) ' cm']),
loglog([0.03 0.03],[low*10^3 high*10^3],'k','linewidth',5), text(0.032,10^3,'95%'),
H(2) = subplot(2,2,3);
polar( deg2rad(dir1(ig)),freq(ig),'gx'),hold on, polar( deg2rad(dir2(ig)),freq(ig),'r+' )
,
title('INFRAGRAVITY'),
% legend('dir1','dir2',-1),
ylabel('Direction, relative to (+) x'),
H(3) = subplot(2,2,4);
polar( deg2rad(dir1(swell)),freq(swell),'gx'),hold on, polar( deg2rad(dir2(swell)),freq(swell),'r+' ),
title('SWELL')
text(-.08,-.27,['DIR2 = ' num2str(dir2_swell,3) ' +/- ' num2str(spread2_swell/2,2)])

figure(2), % pressure & velocity
K(1) = subplot(2,1,1);
loglog ( freq, PP, 'b',  freq, UUpres, 'r', freq, VVpres, 'g', freq, UUpres+VVpres, 'y--'), hold on
ylabel('energy density (cm^2/Hz)') 
axis( [ 0.002 0.7 0.2 10^4] ),
legend('PP','(\omega/gk) UU','(\omega/gk) VV','(\omega/gk) (UU+VV)',1),
title([ file ',  bandwidth = ' num2str(bandwidth,2) ',  DOF = ' num2str(DOF)])
loglog([0.03 0.03],[low.*10^3  high.*10^3],'k','linewidth',5), text(0.032,10^3,'95%')
K(2) = subplot(4,1,3);
semilogx( freq,cohPU,'r',  freq,cohPV, 'g', [0.002 0.7], [SIG SIG],'k:')
ylabel('Coherence') 
axis( [ 0.002 0.7 0 1] ),
legend('PU','PV',1),text(0.0025,SIG+0.05,'95% level')
K(3) = subplot(4,1,4);
semilogx ( freqPUsig,phPUsig,'ro',  freqPVsig,phPVsig, 'gs', [0.002 0.7], [0 0],'k--')
ylabel('Phase (deg)'), xlabel('frequency (Hz)') 
axis( [ 0.002 0.7 -180 180] ),
legend('PU','PV',1),
set(K(3),'YTick',[-180 -90 0 90 180]);

% DIR plot...
% figure(4), % directions and spreads 
% L = semilogx( freq,dir1,'r', freq,dir2,'b', freq,(dir1+spread1./2),'r:', freq,(dir1-spread1./2),'r:', freq,(dir2+spread2./2),'b:', freq,(dir2-spread2./2),'b:');
% ylabel('wave direction, zero is along (+) x'), xlabel('frequency (Hz)'),
% legend('dir1','dir2',4), axis([0 inf -180 180])
% title([ file ',  bandwidth = ' num2str(bandwidth,2) ',  DOF = ' num2str(DOF)])
% text(0.07,100,['DIR2 = ' num2str(dir2_swell,3) ' +/- ' num2str(spread2_swell/2,2)])

% "Z" test plot...
% figure(5), % quality control (including higher freq)
% J = loglog ( freq, PP,'b', freq, (UUpres + VVpres),'m' ); 
% xlabel('frequency (Hz)'), ylabel('spectral energy density (cm^2/Hz)') 
% legend('PP (bottom)','UU + VV (converted to p)',2),
% title([ file ',  bandwidth = ' num2str(bandwidth,2) ',  DOF = ' num2str(DOF)])

% figure(3), % energy fluxes
% loglog(freq,posX,'c',freq,negX,'c--',freq,posY,'y',freq,negY,'y--','linewidth',3),axis([10^-3 0.5 0 inf])
% legend('(+) x','(-) x','(+) y','(-) y',2)
% ylabel('energy flux (cm^3 s^{-1} Hz^{-1})') 
% title([ file ',  bandwidth = ' num2str(bandwidth,2) ',  DOF = ' num2str(DOF)])

% little detail (re-order active figures)
figure(1)

else end
%-----------------------------------------------------------------------------
