% function T0 = sonic_integral_timescale(y,fs)
% function T0 = sonic_integral_timescale(y,fs,'method',methodstring)
% function T0 = sonic_integral_timescale(y,fs,'method',methodstring,DO_DEBUG,true/false)
%
% calculates the integral timescale for a turbulent signal y sampled at a
% frequency fs. The value returned depends on the method, which can be
% specified as methodstring
% 
% 'zc' = zero-crossing; The time is the integral of the autocorrelation
% function of the turbulent component from t=0 to time at which the
% correlation drops to zero.
%
% 'te' = time to reach 1/e; the time is the time at which the correlation
% reaches 1/e. This may be more robust.
%
% written by Andy Clifton, September 2010

function T = sonic_integral_timescale(y,fs,varargin)

% default values
DO_DEBUG = 0;
method = 'te';

% check the inputs
optargin = size(varargin,2);
stdargin = nargin - optargin;
if optargin > 0
        for i = 1:2:optargin
            eval(['optargin{i} = ' optargin{i+1} ';']);
        end
end

% get the correlation function Rxx;
[c,lags]= xcorr(y-nanmean(y),'coeff');

% find the max (or the middle of the spectra
[maxval,mi]=max(c);
% cut it down
c = c(mi:end);
lags = lags(mi:end);
switch lower(method)
    case 'zc'
        % find the first zero crossing
        deltac=diff(c./abs(c)); % look for a sign change
        nzero = min(find(deltac<0));
        if isempty(nzero)
            T = NaN;
        else
            T = trapz([1:nzero]./fs,c(1:nzero));
        end
    case 'te'
        % find the time at which the correlation drops to 1/e
        nte = min(find(c<(1/exp(1))));
        if isempty(nte)
            T = NaN;
        else
            T = nte / fs;
        end
end


%% debug
if DO_DEBUG
    figure
    plot(lags,c,'k')
    grid on
    hold on
    plot(xlim,[1/exp(1) 1/exp(1)],'k--')
    switch method
        case 'zc'
            plot([nzero nzero],ylim,'r--')
        case 'te'

            plot([nte nte],ylim,'r--')
    end
end