% process Vector data from TTT in daily blocks (continuous sampling)
% 32 Hz sampling
%

filebase = 'TTM_PNNLvector_Jun2012_sampleburst';
fpath = './';

data = load([filebase '.dat']);

u = data(:,3);
v = data(:,4);
w = data(:,5);

%%

mag = sqrt(u.^2 + v.^2);


%%
slowt = [0:60*32:length(mag)];
slowmag = interp1( [1:length(mag)], mag, slowt );


%%

figure(1), clf

subplot(2,1,1),
plot(slowt./(32*60*60),slowmag,'r-'), hold on,
for i = 1:10:(60*24), 
    plot(slowt(i)./(32*60*60),slowmag(i),'k.'),
end
set(gca,'FontSize',16,'FontWeight','demi')
axis([0 24 0 2]),
ylabel('V [m/s]'), xlabel('Time [hrs]'),
title('17 Feb 2011','FontSize',20,'FontWeight','bold')


subplot(2,1,2),
plot([0:(32*60*10)]./(32*60),mag([2.2e6:(2.2e6+32*60*10)]),'r-'), hold on
plot([0:(32*60*1)]./(32*60),mag([2.2e6:(2.2e6+32*60*1)]),'k-'), hold on
set(gca,'FontSize',16,'FontWeight','demi')
axis([0 (32*60*10)./(32*60) 0 2]),
ylabel('V [m/s]'), xlabel('Time [min]'),



