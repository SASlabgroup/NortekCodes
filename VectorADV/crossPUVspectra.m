
function [ freq, cohP1P2, phP1P2, SSE1, SSE2, P1P2, Hsig_swell_1, Hsig_swell_2, Hsig_ig_1, Hsig_ig_2 ] = crossPUVspectra( mmddhhhh, id1 , id2 , window_length, merge, plottoggle);

%
% Matlab function to calculate and plot spectra from two PUV sensors
%    using time series stored in a "mmddhhhh.ID.mat" files
%    (where all measurements are in one file... see RAW-processing dir)
%
% Inputs are strings 'mmddhhhh' and instrument ids, the window size (in seconds), and 
%    the integer number of f-bands to merge, and a toggle switch for plotting 
% 
%      >> crossPUVspectra( mmddhhhh, id1, id2, window, merge, plottoggle );
%
%    ** note that the window size should be 2^N for FFT efficiency **
%
% Output is frequency (Hz), coherence, phase (deg), 
%    power spectral density of sea surface elevation 1, 
%    power spectral density of sea surface elevation 2, cross-power spectral 
%    density of the two pressure records, significant wave height 1, 
%    significant wave height 2, infragravity significant wave height 1, 
%    infragravity significant wave height 2, 
%    
% Function also makes a PLOT of the two SSE spectra,
%    and the coherence and phase (from cross-spectra) 
% 
% Currently this function only deals with pressure data
%
% Use similar script "PUVspectra.m" to look at data from a single PUV
%
% J. Thomson, September 2003

warning off MATLAB:divideByZero; close all;

% use these two lines to run this as script instead of a function:
	% file = input('enter filename (including path): ', 's');
	% merge = input('enter number of f-bands to merge: ');

% LOAD DATA -------------------------------------------------------------------
% .mat files (with appropriate dir path) with entire time series
% convention is to keep variables lowercase in time domain,
% and then use capitals in frequency domain
load(['~/NCEX/DATA/Triton/' mmddhhhh(1:4) '/' mmddhhhh '.' id1 '.mat' ] ); p1=p;
pts = length(p1);                 % record length in data points
load(['~/NCEX/DATA/Triton/' mmddhhhh(1:4) '/' mmddhhhh '.' id2 '.mat' ] ); p2=p;
if pts ~= length(p2), disp('PROBLEM: records are different lengths'), else end
secs = pts / sample_rate;        % record length in seconds
%------------------------------------------------------------------------------



% WINDOW THE DATA (use 75 percent overlap)--------------------------------------
% WINDOW LENGTH SHOULD BE 2^n FOR FFT EFFICIENCY
% window_length = ;                  % window length in secs (now an input to function)
w = sample_rate * window_length;     % window length in data points
windows = 4*(pts/w -1)+1;   % number of windows, the 4 comes from a 75% overlap
% loop to create a matrix of time series, where COLUMN = WINDOW 
for q=1:windows, 
	p1_window(:,q) = p1(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
	p2_window(:,q) = p2(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
end
%-------------------------------------------------------------------------------



% WATER DEPTH (use later for depth correction)----------------------------------
doff = 78;                        % distance of central transducer off bottom
                                  % (measured by diver during deployment)
doffu = doff + 10;                % distance (cm) of sample volume off bottom 
doffp = doff - 23.25;             % distance (cm) of p-gage off bottom
depth1 = mean(p1_window) + doffp;   % avg p (cm) depth for that column (window)
depth2 = mean(p2_window) + doffp;   % avg p (cm) depth for that column (window)
%-------------------------------------------------------------------------------



% DETREND THE WINDOWED DATA-----------------------------------------------------
% remove the mean of each window
p1_window_nomean = p1_window - ( ones(w,1) * mean(p1_window)  ) ; 
p2_window_nomean = p2_window - ( ones(w,1) * mean(p2_window)  ) ; 
% loop to remove quadratic trend (attempt to reduce tidal leakage)
secs = [sample_rate^-1:sample_rate^-1:(w/sample_rate)]';
for q=1:windows
quadfitp1 = polyfit(secs,p1_window_nomean(:,q),2);
p1_window_detrend(:,q) = p1_window_nomean(:,q) - (quadfitp1(3) + secs.*quadfitp1(2) + (secs.^2).* quadfitp1(1) );
quadfitp2 = polyfit(secs,p2_window_nomean(:,q),2);
p2_window_detrend(:,q) = p2_window_nomean(:,q) - (quadfitp2(3) + secs.*quadfitp2(2) + (secs.^2).* quadfitp2(1) );
end
%------------------------------------------------------------------------------



% TAPER THE DATA (use a Hanning type window)-----------------------------------
% form taper matrix (columns of taper coef)
taper = sin ([1:w] * pi/w )' * ones(1,windows); 
% taper each window
p1_window_taper = p1_window_detrend .* taper;
p2_window_taper = p2_window_detrend .* taper;
% now find the correction factor (comparing old/new variance)
factp1 = sqrt( var(p1_window_detrend) ./ var(p1_window_taper) );
factp2 = sqrt( var(p2_window_detrend) ./ var(p2_window_taper) );
% and correct for the change in variance
% (mult each window by it's variance ratio factor)
p1_window_ready = (ones(w,1)*factp1).* p1_window_taper;
p2_window_ready = (ones(w,1)*factp2).* p2_window_taper;
% check & report
if abs(  var(p1_window_ready) - var(p1_window_detrend)  ) > 0.1,
  disp('******************************')
  disp('Problem preserving variance variance');
  disp('******************************')
  else end
%note that this non-uniform window has changed the mean to be non-zero. 
  %The variance was preserved, but what about this introduced mean???!?!?
%------------------------------------------------------------------------------



% SPECTRA (FFT)-----------------------------------------------------------------
% calculate Fourier coefs
P1_window = fft(p1_window_ready);
P2_window = fft(p2_window_ready);
% second half of fft is redundant, so throw it out
P1_window( (w/2+1):w, : ) = [];
P2_window( (w/2+1):w, : ) = [];
% throw out the mean (first coef) and add a zero (to make it the right length)  
P1_window(1,:)=[];  P2_window(1,:)=[];  
P1_window(w/2,:)=0; P2_window(w/2,:)=0;  
% POWER SPECTRA (auto-spectra)
PP1_window = ( P1_window .* conj(P1_window) ) ;
PP2_window = ( P2_window .* conj(P2_window) ) ;
% CROSS-SPECTRA  
P1P2_window = ( P1_window .* conj(P2_window)) ;
% definitions for later:
% Cospectra = real(cross-spectra), Quadrature = image(cross-spectra)
% -----------------------------------------------------------------------------



% MERGE FREQUENCY BANDS -------------------------------------------------------
% raw fft has (w/2)^-1  resolution of frequency before merging
% number of bands to merge is an input to function
for i = merge:merge:(w/2) 
	PP1_window_merged(i/merge,:) = mean( PP1_window((i-merge+1):i , : ) );
	PP2_window_merged(i/merge,:) = mean( PP2_window((i-merge+1):i , : ) );
	P1P2_window_merged(i/merge,:) = mean( P1P2_window((i-merge+1):i , : ) );
end
% freq range and bandwidth
n = (w/2) / merge;                         % number of f bands
Nyquist = .5 * sample_rate;                % sampling was 2 hz
bandwidth = Nyquist/n ;                    % freq (Hz) bandwitdh
% find middle of each freq band
freq= 1/(window_length) + bandwidth/2 + bandwidth.*[0:(n-1)] ; 
% -----------------------------------------------------------------------------



% ENSEMBLE AVERAGE THE WINDOWS -------------------------------------------------
% take the average of all windows at each freq-band
% create power density, divide by N*sample rate
PP1 = mean( PP1_window_merged' ) / (w/2 * sample_rate ) ;    
PP2 = mean( PP2_window_merged' ) / (w/2 * sample_rate );
P1P2 = mean( P1P2_window_merged' ) / (w/2 * sample_rate ) ;
%--------------------------------------------------------------------------



% DEPTH CORRECTION -------------------------------------------------------------
% find correction coefs at each f-band 
% to calc sea surface elevation 
for index = 1:n
   f = freq(index);
   fcutoff = 0.25 ;  % cutoff frequency 
		    % (beyond which correction will be too large to trust)
   if f < fcutoff,
   % find k for each f with function wavenumber.m (a Newton-Ralphson iteration)
   % use avg depth of all windows (inexact but not a big deal)
   k1 = wavenumber(f,mean(depth1));
   correction1(index) =  cosh( k1*mean(depth1) ) / cosh(k1*doffp);
   k2 = wavenumber(f,mean(depth2));
   correction2(index) =  cosh( k2*mean(depth2)) / cosh(k2*doffp);
   else correction1(index) = 0; correction2(index) = 0;   end  % beyond fcutoff
end   
% correct pressure for attentuation to get sea-surface elevation
SSE1 = PP1 .* (correction1.^2) ; 
SSE2 = PP2 .* (correction2.^2) ; 
%-------------------------------------------------------------------------------



% COHERENCE & PHASE, etc -------------------------------------------------------
% Definitions from Elgar, 1985, JFM ... gives same result as Wunsch's def.
% Cospectrum & Quadrature:
coP1P2 = real(P1P2);   quP1P2 = imag(P1P2);
% Coherence & Phase at each freq-band
% *** note that it's important to calc this AFTER all merging and ensemble avg.
cohP1P2 = sqrt( (coP1P2.^2 + quP1P2.^2) ./ (PP1.* PP2) );
phP1P2  = 180/pi .* atan2( quP1P2 , coP1P2 );  
% -----------------------------------------------------------------------------



% SPECTRAL WEIGHTED AVERAGES & STATS -----------------------------------------
% find indices of freq bands
ig = find(freq>0.005 & freq<0.03);
swell = find(freq>0.05 & freq<0.25);
% significant wave height (cm):
Hsig_swell_1 = 4 * sqrt( sum( SSE1(swell) * bandwidth ) ); 
Hsig_swell_2 = 4 * sqrt( sum( SSE2(swell) * bandwidth ) );       
Hsig_ig_1 = 4 * sqrt( sum( SSE1(ig) * bandwidth ) );  
Hsig_ig_2 = 4 * sqrt( sum( SSE2(ig) * bandwidth ) );  
% centriod frequency
centroid_ig_1 = sum ( freq(ig).* PP1(ig) ) / sum ( PP1(ig) ) ;
centroid_swell_1 = sum ( freq(swell).* PP1(swell) ) / sum ( PP1(swell) ) ;
centroid_ig_2 = sum ( freq(ig).* PP2(ig) ) / sum ( PP2(ig) ) ;
centroid_swell_2 = sum ( freq(swell).* PP2(swell) ) / sum ( PP2(swell) ) ;
% ----------------------------------------------------------------------------



% DEGREES OF FREEDOM and level of no significant coherence --------------------
% DOF = 2 * (# independent windows) * (# bands merged)
DOF = 2 * pts/w * merge;  
% 95% significance level for zero coherence
SIG = sqrt(6/DOF);
SIGindices = find(cohP1P2 > SIG);
% ------------------------------------------------------------------------------



% PLOTTING ---------------------------------------------------------------------
if plottoggle==1,
figure(1), % sea surface spectra
K(1) = subplot(2,1,1);
loglog ( freq, SSE1, 'b',  freq, SSE2, 'r' )
ylabel('SSE energy density (cm^2/Hz)') 
axis( [ 0.002 0.3 0 inf] ), hold on
legend(id1, id2, 2),
title([ mmddhhhh ', bandwidth = ' num2str(bandwidth,2) ',  DOF = ' num2str(DOF)])
K(2) = subplot(4,1,3);
semilogx( freq,cohP1P2,'g', [0.002 0.3], [SIG SIG],'k:')
ylabel('Coherence') 
axis( [ 0.002 0.3 0 1] ),
K(3) = subplot(4,1,4);
semilogx ( freq(SIGindices),phP1P2(SIGindices),'kx', [0.002 0.3], [0 0],'k--')
ylabel('Phase (deg)'), xlabel('frequency (Hz)') 
axis( [ 0.002 0.3 -180 180] )

else end

