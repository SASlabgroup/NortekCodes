function [StndHeading, StndPitch, StndRoll] = CalcStndHPR(Data,DoPlot)

%   Compute HPR values from AHRS source data but disregarding the Gyro
% component.  This will make the new HPR values be the same as if they were
% recorded with a "standard" compass & 2-axis tilt sensor.  Assumes data
% was Burst and exported from MIDAS so structure is called "Data".
%   DoPlot is logical 1 or 0

HxHyHz_sensor = [Data.Burst_MagnetometerX ...
                 Data.Burst_MagnetometerY ...
                 Data.Burst_MagnetometerZ];
             
Acceleration_sensor = [Data.Burst_AccelerometerX ...
                       Data.Burst_AccelerometerY ...
                       Data.Burst_AccelerometerZ];

[StndPitch, StndRoll] = CalcAD2CPattitude(Acceleration_sensor, Data.Burst_Status);

StndHeading = CalcAD2CPHeading(HxHyHz_sensor, StndPitch, StndRoll, Data.Burst_Status);

if DoPlot
    figure
    plot(StndHeading)
    hold on
    plot(Data.Burst_Heading,'g')
    grid on
    legend('Standard Heading', 'AHRS Heading')
    ylabel('Degrees')
    xlabel('Sample Number')
  
%     Re-scaling Roll first, just for plotting
        pos = find(Data.Burst_Roll>100);
        neg = find(Data.Burst_Roll<100);
        OrigRoll = Data.Burst_Roll;
        OrigRoll(pos) = OrigRoll(pos) - 180;
        OrigRoll(neg) = OrigRoll(neg) + 180;
        
        pos = find(StndRoll>100);
        neg = find(StndRoll<100);
        OrigStndRoll = StndRoll;
        OrigStndRoll(pos) = OrigStndRoll(pos) - 180;
        OrigStndRoll(neg) = OrigStndRoll(neg) + 180;
    figure
    plot(OrigStndRoll)
    hold on
    plot(OrigRoll,'g')
    grid on
    legend('Standard Roll', 'AHRS Roll')
    ylabel('Degrees')
    xlabel('Sample Number')

end

end