function [GPS] = readGPS(fpath, GPSfilelist, date, UTCtolocal, readraw)

% Input variables:
% fpath         = filepath to the GPS file (.csv, or .mat if exists)
% GPSfilelist   = list of filenames of the GPS files, excluding their
%                   extensions (e.g. DC01_GPSa_09Jul2015)
% date          = date of measurement (e.g. 09Jul2015)
% UTCtolocal    = correction from UTC to local time, in hours (for Igiugig,
%                   UTCtolocal = -8)
% readraw       = logic input which forces the re-parsing of the .csv data,
%                   regardless of the existance of the .mat version.
%
% Output structure fields for x GPS files: 
% GPS.GPS_x.gpstime
%          .lat 
%          .lon 
%          .speed 
%          .GPSheading 
%          .date
%

for GPSi = 1:length(GPSfilelist)
    
    GPSfname = GPSfilelist{GPSi};
    
    if isempty(dir([ fpath GPSfname '.mat'])) | readraw==1,
        
        [gi track utcdate utctime localdate localtime miliseconds valid lat northsouth lon eastwest altitude speed GPSheading distance] = textread([fpath GPSfname '.csv'],'%n%s%s%s%s%s%n%s%n%s%n%s%n%n%n%n','delimiter',',','headerlines',1,'bufsize',1600000);
        % Note that this will fail if the Qstarz restarted and has a line
        % of header info somewhere in the middle of the file
        
        % Generate the timestamp of the GPS measurement by concatenating the
        % date and time vectors
        for gi = 1:length(lat),
            gpstime(gi) = datenum([char(utcdate(gi)) ' ' char(utctime(gi)) ] );
        end
        
        gpstime = gpstime' + ( UTCtolocal /24);  % adjust to local time stamp
        speed   = speed * .2777;                % convert from Km/h to m/s
        
        save([fpath GPSfname '.mat'],'gpstime','lat','lon','speed','GPSheading','date')
        
    else
        save temp3
        load([ fpath GPSfname '.mat'])
        load temp3
    end
    
    % Save GPS measurements in a structure
    
    GPS.(['GPS_' num2str(GPSi)]).gpstime    = gpstime;
    GPS.(['GPS_' num2str(GPSi)]).lat        = lat;
    GPS.(['GPS_' num2str(GPSi)]).lon        = lon;
    GPS.(['GPS_' num2str(GPSi)]).speed      = speed;
    GPS.(['GPS_' num2str(GPSi)]).GPSheading = GPSheading;
    GPS.(['GPS_' num2str(GPSi)]).date       = date;
    
end
