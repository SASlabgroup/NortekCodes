function M = CalcHPRMatrix(hh, pitch, roll)

hh = hh - 90;
H = [cosd(hh) sind(hh) 0; -sind(hh) cosd(hh) 0; 0 0 1];
T = CalcTiltMatrix(pitch, roll);
M=H*T;

end


% The reason for the -90 above (line 3) is because we use ENU coordinates.
% Since the X-axis points to North we need to rotate in order to have East
% as the first component, i.e. classical ADCP stuff since ENU seems to be
% the defacto standard around. NED would have been more logical.