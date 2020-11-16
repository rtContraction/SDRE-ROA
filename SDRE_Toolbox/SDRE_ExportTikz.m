function u=SDRE_ExportTikz(name,OS)
% filename         = 'C:/Users/ga48sam/Documents/flat-metrics-paper/IFAC2020/';
% dataPath         = 'C:/Users/ga48sam/flat-metrics-paper/data';
 fileWindows  = '.\..\..\SDRE_Paper\'; %..\..\ICMPaperGit\
 fileLinux    = './../../SDRE_Paper/';
 dataPath         = '.\..\..\SDRE_Paper\data';
relativeDataPath = 'data';
filename=fileLinux;
if strcmp('Windows',OS)
    filename=fileWindows;
end
filename=strcat(filename,name);
matlab2tikz(filename, 'relativeDataPath', ...    
    relativeDataPath, 'dataPath', dataPath,'height', '\fheight', 'width', '\fwidth' );
u=1;
end
%'.\..\..\Writting\ICM2021\'
%%%%%  Post processing %%%%%%
% Add the following commands %
% ylabel near ticks,
% xlabel near ticks
% instead of at={(0\fwidth,0\fheight)},

%%%%  Eliminate time to the figures they dont need it

%%%%%%%%%  For figure 1 %%%%%%%%%%%%%
%legend style={legend cell align=left, align=left, draw=white!15!black,at={(1,1)},anchor=north east}

%%%%%%%%%  For figure 2 %%%%%%%%%%%%%
%legend style={legend cell align=left, align=left, draw=white!15!black,at={(1,0)},anchor=south east}

%%%%%%%%%  For figure 3 %%%%%%%%%%%%%
