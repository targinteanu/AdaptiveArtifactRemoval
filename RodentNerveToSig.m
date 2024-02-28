function [sigNaive, sigVDMT] = RodentNerveToSig(foldername, truncateSig)
% Return signal object from Rodent nerve stim experiments (Kiara). 
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
folder_amp  = [foldername,filesep,'amplifier_data'];
folder_trig = [foldername,filesep,'trigger_data'];
[folder_amp_N,  folder_amp_V]  = NVfolders(folder_amp);
[folder_trig_N, folder_trig_V] = NVfolders(folder_trig);

%% Naive 
if (~isempty(folder_amp_N)) & (~isempty(folder_trig_N))
    folderAmp = [folder_amp,filesep,folder_amp_N];
    folderTrig = [folder_trig,filesep,folder_trig_N];
    sigNaive = sigFromAmpTrigFolders(folderAmp, folderTrig, truncateSig);
else
    sigNaive = [];
end

%% VDMT 
if (~isempty(folder_amp_V)) & (~isempty(folder_trig_V))
    folderAmp = [folder_amp,filesep,folder_amp_V];
    folderTrig = [folder_trig,filesep,folder_trig_V];
    sigVDMT = sigFromAmpTrigFolders(folderAmp, folderTrig, truncateSig);
else
    sigVDMT = [];
end

%% main body helper: use on N and V

function sig = sigFromAmpTrigFolders(ampFolder, trigFolder, truncatesig)

%% load     
d0 = []; g0 = []; % init 

% amplifier/filt -> d
files = dir(ampFolder);
files = files(~[files.isdir]);
for f = 1:length(files)
    files(f).name
    load([ampFolder,'/',files(f).name]);
    d0 = [d0, filt];
    clear filt
end

% trigger/adc -> g
files = dir(trigFolder);
files = files(~[files.isdir]);
for f = 1:length(files)
    files(f).name
    load([trigFolder,'/',files(f).name]);
    g0 = [g0, board_adc_data];
    clear board_adc_data
end

%% truncate 
if truncatesig
    d0 = d0(:,1:3450000);
    g0 = g0(:,1:3450000);
end

%% define timing
Fs = 20000; % samples per second 
dt = 1/Fs; % time step (seconds) 
samples=1:length(d0);
samples = samples-1;
t = samples/Fs;
t = t';

%% fix noise reference 
% each amplitude: 10 pulses 
ampvals = {20:2:120, ...   % 2 uA => 20uA – 120uA
           125:5:300, ...  % 5 uA => 125uA – 300uA
           310:10:500, ... % 10 uA => 310uA – 500uA
           550:50:750};    % 50 uA => 550uA - 750uA
ampvalSz = arrayfun(@(i) length(ampvals{i}), 1:length(ampvals))';
ampvals = cell2mat(ampvals);

g = g0; % reset 
g = g./max(g);
g = g.*(g > .1);
% find actual recorded number of pulses 
[trigval, trigloc] = findpeaks(g, 'MinPeakProminence', .1*max(g));

% split pulses into groups of 10 
nInGrp = 10;
ngrp = length(trigloc)/nInGrp;
pulseFirstLast = zeros(nInGrp,ngrp);
pulseFirstLast(:) = trigloc(:);
pulseFirstLast = pulseFirstLast';
pulseFirstLast = pulseFirstLast(:,[1,end]);
pulsesGap0 = [pulseFirstLast(2:end,1), pulseFirstLast(1:(end-1),2)];
pulsesBnd = mean(pulsesGap0, 2)'; 

if sum(ampvalSz) == ngrp
    % There are the same number of expected and recorded pulses. Proceed
    % assigning amplitude from beginning. 
    startFromEnd = false; 
else
    % There is a mismatch between expected vs recorded number of pulses.
    % Decide how to assign amplitudes. 

pulsesGap = diff(pulsesGap0, [], 2);
[~,ord] = sort(pulsesGap);
ord = ord(1:3); ord = sort(ord);
groupFirstLast = [[0;ord]+1, [ord;ngrp]];
%groupSz = diff(groupFirstLast, [], 2);
groupSz = zeros(size(groupFirstLast,1),1);

for idx = 1:size(groupFirstLast,1)
    groupIdx = groupFirstLast(idx,:);
    pulsesIdx = pulseFirstLast(groupIdx(1):groupIdx(2), :);
    trigidx = (trigloc >= pulsesIdx(1,1)) & (trigloc <= pulsesIdx(end,end));
    groupSz(idx) = sum(trigidx)/10;
end

mismatchSz = ampvalSz - groupSz; 
if (mismatchSz(1) <= 0) & (mismatchSz(end) <= 0)
    % error is in the middle? 
    error('Unable to resolve discrepancy between expected and recorded stimuli.')
elseif mismatchSz(1) <= 0
    % end is cut off but beginning is not. 
    startFromEnd = false;
elseif mismatchSz(end) <= 0
    % beginning is cut off but end is not. 
    startFromEnd = true;
else
    % mismatch at both beginning and end
    if abs(mismatchSz(1)) <= abs(mismatchSz(end))
        startFromEnd = false;
        ampvals = ampvals(abs(mismatchSz(1)):end);
    else
        startFromEnd = true;
        ampvals = ampvals(1:(end-abs(mismatchSz(end))));
    end
end

end

if ~startFromEnd
    bnd1 = 1; idx = 1;
    pulsesBnd = [pulsesBnd, length(g)];
    for bnd2 = pulsesBnd
        g(bnd1:bnd2) = g(bnd1:bnd2) * ampvals(idx);
        idx = idx + 1; bnd1 = bnd2+1;
    end
else
    bnd2 = length(g); idx = 1; 
    ampvals = fliplr(ampvals);
    pulsesBnd = [1, pulsesBnd];
    for bnd1 = fliplr(pulsesBnd)
        g(bnd1:bnd2) = g(bnd1:bnd2) * ampvals(idx);
        idx = idx + 1; bnd2 = bnd1-1;
    end
end

%% processing loaded data into signal object
% "filt" is the neural recording data, which we want to change into "d_unfilt"
% "board_adc_data" is the noise reference data, which we want to change into "g"
d_unfilt = d0;

% introduce delay in signal compared to noise reference 
padT = .005; % delay dur (s); should be at least .0005
padL = floor(padT*Fs); 
dPad = d_unfilt(:,1).*ones(1,padL);
d_unfilt = [dPad, d_unfilt(:,1:(end-padL))];

d_unfilt = (d_unfilt-32768)*0.0003125; %filt in multiple of volt units
d_unfilt = d_unfilt';

%g = [board_adc_data;board_adc_data];
%g = board_adc_data;
g = g';

chA = buildChannelObj('Chan1', 0,0,0, 'Cartesian');
chB = buildChannelObj('Chan2', 0,0,0, 'Cartesian');

% (t and g should have the same number of rows as d_unfilt see repmat)
t = repmat(t, 1, size(d_unfilt,2));
g = repmat(g, 1, size(d_unfilt,2));

%% define parameters for filter and training 
%trainfrac = .02;
%splIdx = floor(trainfrac*size(t,1));
splIdx = 14e4;
tTrainBnd = [t(1), t(splIdx)];

sig = buildSignalObj([], d_unfilt, t, g, Fs, [chA; chB], ...
                     tTrainBnd, tTrainBnd, 2);

end

%% other helpers 
    function [nameN, nameV] = NVfolders(outerfolder)
        % returns the names of subfolders in outerDir corresponding to
        % Naive or VDMT. Will return empty cell if none are found. 

        outerDir = dir(outerfolder);
        outerDir = outerDir([outerDir.isdir]);
        names = {outerDir.name};
        nIdx = false(size(names)); vIdx = false(size(names));

        % Naive / N
        for i = 1:length(nIdx)
            nIdx(i) = contains(names{i}, 'naive', 'IgnoreCase', 'true');
        end
        if ~sum(nIdx)
            for i = 1:length(nIdx)
                nIdx(i) = (names{i}(1) == 'N') | (names{i}(1) == 'n');
            end
        end
        if ~sum(nIdx)
            for i = 1:length(nIdx)
                nIdx(i) = contains(names{i}, 'n', 'IgnoreCase', 'true');
            end
        end

        % VDMT / V
        for i = 1:length(vIdx)
            vIdx(i) = contains(names{i}, 'VDMT', 'IgnoreCase', 'true');
        end
        if ~sum(vIdx)
            for i = 1:length(vIdx)
                vIdx(i) = (names{i}(1) == 'V') | (names{i}(1) == 'v');
            end
        end
        if ~sum(vIdx)
            for i = 1:length(vIdx)
                vIdx(i) = contains(names{i}, 'v', 'IgnoreCase', 'true');
            end
        end

        nameN = names(nIdx); nameV = names(vIdx);
        if (length(nameN) > 1) | (length(nameV) > 1)
            warning('More than one naive or VDMT folder found.')
        end
        if ~isempty(nameN)
            nameN = nameN{1}; 
        end
        if ~isempty(nameV)
            nameV = nameV{1};
        end
    end

end