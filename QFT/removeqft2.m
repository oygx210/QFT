function removeqft2

CurrentDir = pwd;

% determine current directory for QFT 2.0
QFT2File = which('qelmtlistbox.m');
QFT2FilesDir = strrep(QFT2File, 'qelmtlistbox.m');
chdir(QFT2FilesDir);
chdir ..
chdir qftdemos
QFT2DemosDir = pwd;

% return to users previous directory
chdir(CurrentDir);

% create addqft2 file in /local directory
FileID = fopen([matlabroot,'/toolbox/local/addqft2.m'],'w');
fprintf(FileID, 'addpath %s -end\n', QFT2FilesDir);
fprintf(FileID, 'addpath %s -end\n', QFT2DemosDir);
fclose(FileID);
