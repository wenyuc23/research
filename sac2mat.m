function [data,header,kstnm] = sac2mat(sacfile);
%sac2mat  loads a SAC2000 file into Matlab, along with header information.
%  [data,header,kstnm] = sac2mat(sacfile);
%
%  sacfile must be a STRING, eg.  [fff,hhh,name]=sac2mat('dir1/file1.sac')
%
%Output:
%data=	the data from the file
%header=the first 105 header variables, as described on the SAC Home Page
%	http://www.llnl.gov/sac/
%	http://www.llnl.gov/sac/SAC_Manuals/FileFormatPt1.html
%kstnm=	KSTNM, the station name (string)
%	(header and kstnm are optional)
%
%Some common/useful headers:
%header(1:3) DELTA,DEPMIN,DEPMAX
%header(57) DEPMEN
%header(6:7) B,E
%header(8) Event origin time (seconds relative to reference time)
%header(32:39) STLA,STLO,STEL,STDP,EVLA,EVLO,EVEL,EVDP
%header(51:54) DIST,AZ,BAZ,GCARC        
%header(71:76)  NZYEAR,NZDAY,NZHOUR,NZMIN,NZSEC,NZMSEC
%header(80) NPTS
%						Case Bradford, 8 June 2004


sacfid = fopen(sacfile,'r','l');
if sacfid==-1;error='Not a valid path.'
return;end

%The first 70 header variables are floating point
header1 = fread(sacfid,70,'float32');

%The next 35 are integers (after that they're a mix of strings)
header2 = fread(sacfid,35,'int32');

%Outputs the number of points being read to the screen
header=[header1;header2];
npts = header(80);

%kstnm, the station name, is a string stored in the 110th header variable
%kstnm = KSTNM (First 3 letters of the station name)
fseek(sacfid,4*110,-1);
kstnm = char(fread(sacfid,4)');

%Now read the data...
fseek(sacfid,158*4,-1);
data(1:npts) = fread(sacfid,npts,'float32');
status = fclose(sacfid);
return