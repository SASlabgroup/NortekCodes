%% nortek codes
clear all; close all; clc

% Eq from Nortek:  TS = 10log10(10^(Pr/10)) + 40log10(r) + a*alpha*r + PL +
% Gcal  (Eq 1)

% For now I'm going to set alpha = 0 and PL = 0; Just correct these in your
% final code.

% Final equation for Gcal.  In this equation TS is known from the analtical
% value for the sphere. I will update you with a value but for now will use
% -41 dB (TS_sphere).  To give you a better estimate I need to know the pulse
% Gcal = TS - 10log10(10^(Pr/10)) + 40log10(r) + a*alpha*r + PL 

% load file
load('C:\Users\cbassett\Downloads\SWIFT23_SIG_19Mar2021_16_02.mat')
% check the following to make sure it fits with what Nortek is doing
r = echo.Blanking:echo.CellSize:echo.CellSize*(length(echo.EchoSounder(1,:))-1)+echo.Blanking;


%% see that, assuming r is calculated right, the sphere is at about 2.2 m
figure(1)
pcolor([1:344],-r,echo.EchoSounder'), shading flat
axis([0 350 -6 0])
set(gca,'clim',[4500 8300])

%% search data between 2 to 2.4 m to find maximum received power value
% indices associated with those ranges
inds = find(and(r<2.4, r>2));

% now loop over a pings and log max value
for j = 1:length(echo.EchoSounder(:,1))
[maxPr, maxind] = max(echo.EchoSounder(j,inds))
tarrange(j) = r(inds(1) + maxind -1); % record target max range
Pr(j) = maxPr./100; % note I have divided Pr by 100 because I believe the Nortek units
% are dB*100
end

%%
pings = 1:length(echo.EchoSounder(:,1));
figure(2)
subplot(211)
plot(pings, tarrange,'k')
xlabel('Ping'), ylabel('r [m]')
subplot(212)
plot(pings, Pr,'k')
xlabel('Ping'), ylabel('Pr')
axis([1 max(pings) 50 90])

%% set percentile limit
uplim = prctile(Pr,98); % get the 98th percentile target limit - this is somewhat arbtitrary
% in this case it actually doesn't matter at all but generally having the
% average of a few targets is better so we choose ones on the upper end
% where the taget is on-axis. Sometimes it is also a good idea to remove a
% couple at the higher end where perhaps an animal could have made it an
% outlier. 

useinds = find(Pr >= uplim);

goodranges = tarrange(useinds);
Prgood = Pr(useinds)
goodpings = pings(useinds);

figure(3)
subplot(211)
plot(goodpings, goodranges,'k')
xlabel('Ping'), ylabel('r [m]')
subplot(212)
plot(goodpings, Prgood,'k')
xlabel('Ping'), ylabel('Pr')
axis([1 max(pings) 50 90])

%% Now I'm going to lump all but Gcal in on the right hand side of Eq 1 to get a TS measured (well, it would be if Gcal = 0)

alpha = 0.31; % put in actual PL and alpha values here. Note we were in fresh water
PL = 0

for j = 1:length(Prgood)
TSmeas(j) = Prgood(j) + 40*log10(goodranges(j)) + 2*alpha*goodranges(j) + PL;    
end

% Calculate mean TS measure
meanTSmeas = 10.*log10(mean(10.^(TSmeas./10)));

%% Now calculate Gcal
% Gcal = TS_sphere - meanTSmeas
TS_sphere = -41;

% this is a Gcal
Gcal =  TS_sphere - meanTSmeas;