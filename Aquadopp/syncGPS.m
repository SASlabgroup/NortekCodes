function GPS = syncGPS(time, GPS)

% Inputs:
% time: 1Hz timestamp of velocity measurements (Aquadopp, Vector, or 
% Signature), as read from the .sen file.
% GPS: GPS structure created by readGPS function
% 
% Output:
% GPS structure with additional "mean" field - the 1-second time synced and
% averaged output of all specified GPS units

% Loop through all fields of the GPS data structure
for GPSi = 1 : size(fields(GPS),1)
    
    fieldname   = ['GPS_' num2str(GPSi)];
    
    for ti = 1:(length(time)-1),
        
        % For every 1 second of the velocity instrument timestamp, identify
        % the corresponding GPS measurements and calculate the mean values.
        % The values for each GPS are stored in columns which correspond to 
        % the name of the structure field. A NaN is returned if no GPS
        % reading is available for that second.
        
        thissec     = find( GPS.(fieldname).gpstime >= time(ti) ...
            &  GPS.(fieldname).gpstime < time(ti+1) );
        
        % GPS
        lat_matrix(ti,GPSi)         = mean( GPS.(fieldname).lat(thissec) );
        lon_matrix(ti,GPSi)         = mean( GPS.(fieldname).lon(thissec) );
        speed_matrix(ti,GPSi)       = mean( GPS.(fieldname).speed(thissec) );
        GPSheading_matrix(ti,GPSi)  = mean( GPS.(fieldname).GPSheading(thissec) );
        stdspeed_matrix(ti,GPSi)    = std( GPS.(fieldname).speed(thissec) );
        
    end
       
    % Calculate the mean of each minute's measurements by all available GPS 
    % units. 
    GPS.mean.gpstime    = time;
    GPS.mean.lat        = mean(lat_matrix,2);
    GPS.mean.lon        = mean(lon_matrix,2);
    GPS.mean.GPSheading = mean(GPSheading_matrix,2);
    GPS.mean.speed      = mean(speed_matrix,2);
    GPS.mean.stdspeed   = mean(stdspeed_matrix,2);
    
end

