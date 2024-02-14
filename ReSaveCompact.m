%% setup
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*_sigObj.mat');
cd(codedir);

%% overwrite all files with compact version 
for f = 1:length(files)
    curfilename = [folder,'\',files(f).name]
    sig = loadSignalObj(curfilename);
    saveSignalCompact(curfilename);
    clear sig
end