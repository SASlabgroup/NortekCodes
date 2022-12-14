% Sample code for calling:
% 1) 'readGPS' to read GPS data from either .csv or .mat (if exists)
% 2) 'syncGPS' to synchronize the GPS data from a number of GPS units (of
% any sampling frequency, to the 1Hz timestamp of the .sen file output by
% Nortek acoustic instruments.

ID          = '01';
date        = '09Jul2015';  % date (and directory) to process
UTCtolocal  = -8;           % hours to adjust from UTC to local time stamps (signed)
readraw     = 1;            % binary flag (0 or 1) to force read in of raw data

fpath       = ['./' date '/'];          % path to directory with that days data
AQDfname    = ['DC' ID '_AQD_' date];   % filenaming convention (must be exact)
GPSfname_1  = ['DC' ID '_GPSa_' date];  % filename of GPS file #1
GPSfname_2  = ['DC' ID '_GPSb_' date];  % filename of GPS file #2


% Read in velocity instrument to get time stamp for sycronisation
sen         = load( [ fpath AQDfname '.sen' ]);
month       = sen(:, 1);
day         = sen(:, 2);
year        = sen(:, 3);
hour        = sen(:, 4);
minute      = sen(:, 5);
second      = sen(:, 6);
time        = datenum( year, month, day, hour, minute, second);

% List file names of GPS to be read in a cell array
GPSfilelist = {GPSfname_1, GPSfname_2}; % create cell array

% Call function to read GPS data from all GPS files specified in
% GPSfilelist
[GPS]       = readGPS(fpath, GPSfilelist, date, UTCtolocal, readraw);

% Call function to synchronize GPS readings with the 1Hz time stamp of
% velocity instrument
[GPS]       = syncGPS(time, GPS);

gpstime     = GPS.mean.gpstime;
lat         = GPS.mean.lat;
lon         = GPS.mean.lon;
speed       = GPS.mean.speed;
GPSheading  = GPS.mean.GPSheading;
