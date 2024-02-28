function [sigNaive, sigVDMT, trialNames] = RodentNerveToSigWrapper(foldername, truncateSig)
% Return signal object from Rodent nerve stim experiments (Kiara). 
% Wrapper for various subfolder cases 
% 
% Inputs
%   foldername: name of folder in which either trial folders or amplifier_data and trigger_data are located 
%               if empty or blank, folder will be selected with ui 
%   truncateSig: if true, will shorten the signal length and select one
%                channel. Use if previewing or debugging to save time by 
%                avoiding filtering the entire signal. Default = false
% 
% Outputs 
%   sig: signal object 

if nargin < 2
    truncateSig = false;
    if nargin < 1
        foldername = [];
    end
end
if isempty(foldername)
    foldername = uigetdir;
end

% removes '.' and '..'
removedots = @(mydir) mydir((~strcmp({mydir.name},'.')) & ...
                            (~strcmp({mydir.name},'..')));

%% organize data from file 
subfolders = dir(foldername);
subfolders = removedots(subfolders);
% subfolders may contain "trial1", "trial2", etc
if contains([subfolders.name], 'trial', 'IgnoreCase',true)
    % treat each trial separately 
    sigNaive = []; sigVDMT = []; trialNames = {};
    for trl = 1:length(subfolders)
        subfoldername = subfolders(trl).name;
        if contains(subfoldername, 'trial', 'IgnoreCase',true)
            [N,V] = RodentNerveToSig([foldername,filesep,subfoldername],truncateSig);
            sigNaive = [sigNaive; N]; sigVDMT = [sigVDMT; V];
            trialNames = [trialNames; subfoldername];
        end
    end
else
    % subfolders are not trials; assume they are amp and trig data
    [sigNaive, sigVDMT] = RodentNerveToSig(foldername, truncateSig);
    trialNames = {''};
end

end