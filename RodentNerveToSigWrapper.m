function [sigNaive, sigVDMT] = RodentNerveToSigWrapper(foldername, truncateSig)
% Return signal object from Rodent nerve stim experiments (Kiara). 
% Wrapper for various subfolder cases 
% 
% Inputs
%   foldername: name of folder in which amplifier_data and trigger_data are located 
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

end