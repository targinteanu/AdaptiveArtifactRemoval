% TDTPATH = 'TDTMatlabSDK';
% addpath(genpath(TDTPATH));
% data = TDTbin2mat('Rodent SSEP Data/AC5-230830-130841');

folder = uigetdir;
files = dir(folder);
%check if mkdir already exists 
if isfolder([folder,'/','saved_signal_obj']) | isfolder([folder,'\','saved_signal_obj'])
    error('already saved')
else
    mkdir(folder,'saved_signal_obj');
end

codedir = cd;
savedir = [folder,filesep,'saved_signal_obj',filesep];

if contains(folder,'day') % for Kiara files
    for idx = 3:size(files,1)
        %        cd(folder); %brings loop back to the outer folder
        fname = [folder,filesep,files(idx,1).name];
        if isfolder(fname)
            try
                [sigN, sigV, trlName] = RodentNerveToSigWrapper(fname);
                for trl = 1:length(trlName)
                    sigNt = sigN(trl); sigVt = sigV(trl);
                    if ~isempty(sigNt.Times)
                        saveSignalObj([savedir,files(idx,1).name,'_',trlName{trl},'_Naive','_saved'],sigNt);
                    end
                    if ~isempty(sigVt.Times)
                        saveSignalObj([savedir,files(idx,1).name,'_',trlName{trl},'_VDMT','_saved'],sigVt);
                    end
                end
            catch ME
                warning([fname,' not saved: ',ME.message]);
            end
        end
    end
cd(codedir);

else %for Mingfeng files 
TDTPATH = 'TDTMatlabSDK';
addpath(genpath(TDTPATH));

for index = 3:size(files,1)
%    cd(folder); %brings loop back to the outer folder
    fname = files(index,1).name;
    if isfolder([folder,filesep,fname])
%        cd(codedir);

        sig = RodentSSEPtoSig([folder,filesep,fname]);

        saveSignalObj([savedir,files(index,1).name,'_saved'],sig);
    end
end

end