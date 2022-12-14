% status=bin2mat(PathAndFile[,ShowOutput])
%
% Convert binary files from Nortek instruments into Matlab *.mat files and return status.
% If status equals 1 the conversion is OK. Supported binary files:
% -Aquadopp files (*.aqd)
% -Aquadopp Profiler files (*.prf)
% -Vector (*.vec) [not if the Vector is set up to "sample on sync"]
%
% Example:
% bin2mat('C:\Nortek\Data\test.aqd')
% will read the Aquadopp file "test.aqd" and save it as "C:\Nortek\Data\test.mat".
%
% If bin2mat is called with 1 argument, there will be no screen output.
% If it is called with 2 arguments, a summary of the conversion will be displayed.
% The variable "ShowOutput" can have any value. Both
% -bin2mat('C:\Nortek\Data\test.aqd',134623) and
% -bin2mat('C:\Nortek\Data\test.aqd','abcd') will produce screen output.
%
% Running the program will produce the variable "README", which explains the other variables.
