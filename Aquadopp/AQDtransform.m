% script to rotate Aquadopp data 
%
%% XYZ to ENU (or XYU)

% orientation data
H = 0;   % set to zero for XYZ to XYU (apply pitch and roll, but not heading, information), as in Duck FRF fixed Aquadopp, Sep 2010
P = pitch;
R = roll;

u = v1 .* [ (cos(H).*cos(P)) * ones(1,cells) ] - v2 .* [(cos(H).*sin(P).*sin(R)-sin(H).*cos(R)) * ones(1,cells)] - v3 .* [(cos(H).*cos(R).*sin(P)+sin(H).*sin(R)) * ones(1,cells)];

v = - v1 .* [ (sin(H).*cos(P)) * ones(1,cells) ] + v2 .* [(sin(H).*sin(P).*sin(R)+cos(H).*cos(R))* ones(1,cells) ] + v3 .* [(sin(H).*sin(P).*cos(R)-cos(H).*sin(R)) * ones(1,cells) ];

w =  v1 .* [ sin(P) * ones(1,cells) ]    + v2 .* [ (sin(R).*cos(P)) * ones(1,cells) ]     + v3 .* [ (cos(P).*cos(R)) * ones(1,cells) ];


%% XYZ to beam

T = [  6461 -3232 -3232;    0 -5596 5596 ;  1506 1506 1506; ] ./4096;

for i=1:cells,
    beambin = (inv(T) * [v1(:,i) v2(:,i) v3(:,i)]'); % vel x time at single bin
    b1(:,i) = beambin(1,:);
    b2(:,i) = beambin(2,:);
    b3(:,i) = beambin(3,:);
end

% project cell range to alongbeam axis
zbeam = z / cos(deg2rad(25));
