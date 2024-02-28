%% setup
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*_sigObj.mat');
cd(codedir);

%% overwrite all files with compact version 
for f = 1:length(files)
    curfilename = [folder,filesep,files(f).name]
    sig = loadSignalObj(curfilename);
    curfilename = curfilename(1:(end-11));
    saveSignalCompact(curfilename, sig);
    clear sig
end