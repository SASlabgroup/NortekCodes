function heading = CalcAD2CPHeading(HxHyHz_sensor, pitch, roll, Status)


% Assumes data exported from MIDAS

   % *** Start of magnetometer recalibration ***


   % *** End of magnetometer recalibration ***


%    orientation = bitand(bitshift(uint32(Status), -25),7);
   %HxHyHz_new = reorientAD2CPaxis(HxHyHz_sensor, orientation);
orientation = 5;    % Hard code orientation to ZDOWN to for down-looking Sig1000 from CB2018 deployment with AHRS
   
   if orientation == 4,
       % Upwards looking instrument
       HxHyHz_inst = HxHyHz_sensor;
   else
       % Downwards looking instrument (rotation around x-axis)
       HxHyHz_inst(:,1)   =  HxHyHz_sensor(:,1);
       HxHyHz_inst(:,2:3) = -HxHyHz_sensor(:,2:3);
   end

   % Calculate the magnetic vector in the Earth aligned coordinate system

   %{
   %Matrix versions - identical maths to below
   M_tilt = CalcTiltMatrix(pitch, roll);
   HxHyHz_Earth = M_tilt*HxHyHz_inst';
   %}


      sinPitch = sind(pitch);
      cosPitch = cosd(pitch);
      sinRoll = sind(roll);
      cosRoll = cosd(roll);
      sinPsinR = sinPitch.*sinRoll;
      sinPcosR = sinPitch.*cosRoll;

      earthHx = HxHyHz_inst(:,1).*cosPitch - HxHyHz_inst(:,2).*sinPsinR - HxHyHz_inst(:,3) .* sinPcosR;
      earthHy = HxHyHz_inst(:,2).*cosRoll - HxHyHz_inst(:,3).*sinRoll;

      heading = atan2(earthHy, earthHx)*180/pi;

   % Calculate the heading
   %heading = atan2(HxHyHz_Earth(2),HxHyHz_Earth(1))*180/pi;


%    if heading < 0,
%        heading = heading + 360;
%    end

for ii = 1:length(heading)
    if heading(ii) < 0
        heading(ii) = heading(ii) + 360;
    else
        heading(ii) = heading(ii);
    end
end

end