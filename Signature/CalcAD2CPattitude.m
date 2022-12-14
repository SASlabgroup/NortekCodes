function [pitch, roll] = CalcAD2CPattitude(Acceleration_sensor, Status)

% Assumes data exported from MIDAS!

%    orientation = bitand(bitshift(uint32(Status), -25),7);
   %HxHyHz_new = reorientAD2CPaxis(HxHyHz_sensor, orientation);
   orientation = 5; % Hard code orientation to ZDOWN to for down-looking Sig1000 from CB2018 deployment with AHRS

   if orientation == 4,
       % Upwards looking instrument
       Acceleration = Acceleration_sensor;
   else
       % Downwards looking instrument (rotation around x-axis)
       Acceleration(:,1)   =  Acceleration_sensor(:,1);
       Acceleration(:,2:3) = -Acceleration_sensor(:,2:3);
   end

   Acceleration = -Acceleration;
   
   roll = atan2(Acceleration(:,2), Acceleration(:,3));    

   cosr = cos(roll);
   sinr = sin(roll);

   d1 = Acceleration(:,2).*sinr;
   d2 = Acceleration(:,3).*cosr;


   pitch = atan(-Acceleration(:,1)./(d1+d2));

   roll = roll/pi*180;
   pitch = pitch/pi*180;
   
end