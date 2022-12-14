bin2mat:
Converts binary files from Nortek instruments into MATLAB mat-files (*.mat).
These instruments (files) are supported:
-Aquadopp (*.aqd)
-Aquadopp Profiler (*.prf)
-Vector (*.vec) [not if the Vector is set up to "sample on sync"]

------------------------------------------------------------------------------------------------------------
To "install" bin2mat:
------------------------------------------------------------------------------------------------------------
MATLAB:
Copy the files bin2mat.dll and bin2mat.m to a directory on MATLAB's search path (e.g. {MATLAB}\bin), 
or copy them to any directory and add this directory to MATLAB's search path. 
(copy to {ANY_DIR} and then write "addpath {ANY_DIR}" in the MATLAB command prompt )

EXE-FILE:
Copy bin2mat.exe to any directory.

------------------------------------------------------------------------------------------------------------
To get help:
------------------------------------------------------------------------------------------------------------
MATLAB:
Write "help bin2mat" in MATLAB, or read the file "bin2mat.m".
EXE-FILE:
Type "bin2mat /?"

------------------------------------------------------------------------------------------------------------
BE AWARE:
------------------------------------------------------------------------------------------------------------
Processing a large file might take some time.
On a 1.5 GHz Pentium 4 with 512 MB ram these are examples of the time spent:
----------------------
file
size	time	speed
(MB)	(s)	(MB/s)
----------------------
84.9	76.4	1.11
23.0	18.4	1.26
19.8	16.8	1.18
11.6	 9.6	1.21
 0.3	 0.3	1.18

------------------------------------------------------------------------------------------------------------
Help texts:
------------------------------------------------------------------------------------------------------------
1)    MATLAB            » help bin2mat.m
2)    EXE               C:\Kristoffer\Cpp\bin2mat>bin2mat /?
3)    README
3a)   Aquadopp          » README
3b)   Aquadopp Profiler » README
3c)   Vector            » README

----------
1) MATLAB:
----------

» help bin2mat.m

  bin2mat(PathAndFile[,ShowOutput])
 
  Convert binary files from Nortek instruments into Matlab *.mat files.
  Supported binary files:
  -Aquadopp files (*.aqd)
  -Aquadopp Profiler files (*.prf)
  -Vector files (*.vec)
 
  Example:
  bin2mat('C:\Nortek\Data\test.aqd')
  will read the Aquadopp file "test.aqd" and save it as "C:\Nortek\Data\test.mat".
 
  If bin2mat is called with 1 argument, there will be no screen output.
  If it is called with 2 arguments, a summary of the conversion will be displayed.
  The variable "ShowOutput" can have any value. Both
  -bin2mat('C:\Nortek\Data\test.aqd',134623) and
  -bin2mat('C:\Nortek\Data\test.aqd','abcd') will produce screen output.
 
  Running the program will produce the variable "README", which explains the other variables.

------------
2) EXE-FILE:
------------

C:\Kristoffer\Cpp\bin2mat>bin2mat /?
-------------------------------------------------
bin2mat

(1)bin2mat
(2)bin2mat Input
(3)bin2mat Input Output

Read binary file 'Input' and save it as the mat-file 'Output'
'Input' and 'Output' are full or relative Path+Filename for the input and output files

(1):    The program will prompt for input and output files.
        Leaving output empty [just press enter] will make output=input, but with another postfix
        for the output file (*.mat)
(2)&(3) Start program from command line.
(2)     Read 'input' and save in input directory with same file name but with postfix '*.mat'
(3)     Read 'input' and save as 'output'

Example:
>bin2mat C:\Nortek\Data\afile.aqd
will read C:\Nortek\Data\afile.aqd and save it as
C:\Nortek\Data\afile.mat

Running the program will produce the variable 'ReadMe', which explains the other variables

C:\Kristoffer\Cpp\bin2mat>

-------------------------
3a) README (if Aquadopp):
-------------------------

» README

README =

 Aquadopp   
 --------
 Diagnostics: (D)

 Damp         Amplitude data (counts)     NumBeams x NumDiagnosticsSamples
 Dbattery     Battery voltage (Volt)      1 x NumDiagnosticsSamples
 Ddnum        Datenum (ref Matlab datenum)1 x NumDiagnosticsSamples
 Derror       Error code                  1 x NumDiagnosticsSamples
 Dheading     Compass heading (deg)       1 x NumDiagnosticsSamples
 Dpitch       Pitch (deg)                 1 x NumDiagnosticsSamples
 Dpressure    Pressure  (m)               1 x NumDiagnosticsSamples
 Drecords     Number of samples to follow 1 x NumDiagnosticsSeries
 Droll        Roll (deg)                  1 x NumDiagnosticsSamples
 Dsoundspeed  Speed of sound (m/s)        1 x NumDiagnosticsSamples
 Dstatus      Status code                 1 x NumDiagnosticsSamples
 Dtemperature Temperature (deg Celcius)   1 x NumDiagnosticsSamples
 Dtime        Time (YYYY,MM,DD,hh,mm,ss)  6 x NumDiagnosticsSamples
 Dvel         Velocity (m/s)              NumBeams x NumDiagnosticsSamples

 Noise (N):

 Namp         Amplitude data (counts)     NumBeams x NumDiagnosticsSeries
 Nbattery     Battery voltage (Volt)      1 x NumDiagnosticsSeries
 Ndnum        Datenum (ref Matlab datenum)1 x NumDiagnosticsSeries
 Nerror       Error code                  1 x NumDiagnosticsSeries
 Nheading     Compass heading (deg)       1 x NumDiagnosticsSeries
 Npitch       Pitch (deg)                 1 x NumDiagnosticsSeries
 Npressure    Pressure  (m)               1 x NumDiagnosticsSeries
 Nroll        Roll (deg)                  1 x NumDiagnosticsSeries
 Nsoundspeed  Speed of sound (m/s)        1 x NumDiagnosticsSeries
 Nstatus      Status code                 1 x NumDiagnosticsSeries
 Ntemperature Temperature (deg Celcius)   1 x NumDiagnosticsSeries
 Ntime        Time (YYYY,MM,DD,hh,mm,ss)  6 x NumDiagnosticsSeries
 Nvel         Velocity (m/s)              NumBeams x NumDiagnosticsSeries

 Velocity (V):

 Vamp         Amplitude data (counts)     NumBeams x NumVelocitySamples
 Vbattery     Battery voltage (Volt)      1 x NumVelocitySamples
 Vdnum        Datenum (ref Matlab datenum)1 x NumVelocitySamples
 Verror       Error code                  1 x NumVelocitySamples
 Vheading     Compass heading (deg)       1 x NumVelocitySamples
 Vpitch       Pitch (deg)                 1 x NumVelocitySamples
 Vpressure    Pressure  (m)               1 x NumVelocitySamples
 Vroll        Roll (deg)                  1 x NumVelocitySamples
 Vsoundspeed  Speed of sound (m/s)        1 x NumVelocitySamples
 Vstatus      Status code                 1 x NumVelocitySamples
 Vtemperature Temperature (deg Celcius)   1 x NumVelocitySamples
 Vtime        Time (YYYY,MM,DD,hh,mm,ss)  6 x NumVelocitySamples
 Vvel         Velocity (m/s)              NumBeams x NumVelocitySamples

 Other:

 avginterval  Average interval (s)   
 clockdeploy  Deployment time, DepDnum=datenum(clockdeploy)
 coordsystem  Coordinate system (ENU / XYZ / beam)
 diaginterval Diagnostics interval (s)
 frequency    Instrument frequency (kHz)
 serialno     Serial number
 transmatrix  Transformation matrix


----------------------------------
3b) README (if Aquadopp Profiler):
----------------------------------

» README

README =

 Aquadopp Profiler
 -----------------
 Velocity Data (V):

 Vamp1        Amplitude data (counts)     NumCells x NumSamples
 Vamp2        Amplitude data (counts)     NumCells x NumSamples
 Vamp3        Amplitude data (counts)     NumCells x NumSamples
 Vbattery     Battery voltage (Volt)      1 x NumSamples
 Vdnum        Datenum (ref Matlab datenum)1 x NumSamples
 Verror       Error code                  1 x NumSamples
 Vheading     Compass heading (deg)       1 x NumSamples
 Vpitch       Pitch (deg)                 1 x NumSamples
 Vpressure    Pressure  (m)               1 x NumSamples
 Vroll        Roll (deg)                  1 x NumSamples
 Vsoundspeed  Speed of sound (m/s)        1 x NumSamples
 Vstatus      Status code                 1 x NumSamples
 Vtemperature Temperature (deg Celcius)   1 x NumSamples
 Vtime        Time (YYYY,MM,DD,hh,mm,ss)  6 x NumSamples
 Vvel1        Velocity (m/s)              NumCells x NumSamples
 Vvel2        Velocity (m/s)              NumCells x NumSamples
 Vvel3        Velocity (m/s)              NumCells x NumSamples

 Other:

 avginterval  Average interval (s)
 clockdeploy  Deployment time, DepDnum=datenum(clockdeploy)
 coordsystem  Coordinate system (ENU / XYZ / beam)
 diaginterval Diagnostics interval (s)
 frequency    Instrument frequency (kHz)
 serialno     Serial number
 transmatrix  Transformation matrix

___1 / ___2 / ___3 means X/Y/Z, E/N/U or beam1/beam2/beam3, depending on 'coordsystem'.

-----------------------
3c) README (if Vector):
-----------------------

» README

README =

Vector
------
Velocity head (H):

Hcorr         Noise correlation              4 beams x NumHeaders
Hdnum         Datenum (ref Matlab datenum)   1 x NumHeaders
Hnoise        Noise amplitude (counts)       4 beams x NumHeaders
Hnrecords     Number of samples to follow    1 x NumHeaders
Htime         Time (YYYY,MM,DD,hh,mm,ss)     6 x NumHeaders

System, time series (ST):

STanain       Analog input                   1 or 2 x NumSamples
STbattery     Battery voltage (Volt)         1 x NumSamples
STdnum        Datenum (ref Matlab datenum)   1 x NumSamples
STheading     Compass heading (deg)          1 x NumSamples
STpitch       Pitch (deg)                    1 x NumSamples
STroll        Roll (deg)                     1 x NumSamples
STsoundspeed  Speed of sound (m/s)           1 x NumSamples
STtemperature Temperature (deg Celcius)      1 x NumSamples
STtime        Time (YYYY,MM,DD,hh,mm,ss)     6 x NumSamples

Velocity, burst-wise (VB):

VBamp1        Amplitude data (counts)        NumBurst x NumSamplesPerBurst
VBamp2        Amplitude data (counts)        NumBurst x NumSamplesPerBurst
VBamp3        Amplitude data (counts)        NumBurst x NumSamplesPerBurst
VBcorr1       Noise correlation              NumBurst x NumSamplesPerBurst
VBcorr2       Noise correlation              NumBurst x NumSamplesPerBurst
VBcorr3       Noise correlation              NumBurst x NumSamplesPerBurst
VBdnum        Datenum (ref Matlab datenum)   NumBurst x NumSamplesPerBurst
VBvel1        Velocity (m/s)                 NumBurst x NumSamplesPerBurst
VBvel2        Velocity (m/s)                 NumBurst x NumSamplesPerBurst
VBvel3        Velocity (m/s)                 NumBurst x NumSamplesPerBurst

Velocity, time series (VT):

VTamp         Amplitude data (counts)        NumBeams  x NumSamples
VTanain       Analog input                   1 or 2  x NumSamples
VTburst       Burst counter                  1 x NumSamples
VTcorr        Noise correlation              NumBeams x NumSamples
VTcount       Sample counter                 1 x NumSamples
VTdnum        Datenum (ref Matlab datenum)   1 x NumSamples
VTpressure    Pressure  (m)                  1 x NumSamples
VTtime        Time (YYYY,MM,DD,hh,mm,ss)     6 x NumSamples
VTvel         Velocity (m/s)                 NumBeams x NumSamples

Other:

avginterval   Average interval (s)
clockdeploy   Deployment time, DepDnum=datenum(clockdeploy)
coordsystem   Coordinate system (ENU / XYZ / beam)
diaginterval  Diagnostics interval (s)
frequency     Instrument frequency (kHz)
sampfreq      Sample frequency (Hz)
serialno      Serial number
transmatrix   Transformation matrix

___1 / ___2 / ___3 means X/Y/Z, E/N/U or beam1/beam2/beam3, depending on 'coordsystem'.
