% compute the phase space of a   
% _filtered_sigObj.mat file.
% Save the statistics of each as a spreadsheet. 
% Meant to work with rodent SSEP data (Mingfeng). 

%% setup before loop 
clear
[fn,fp] = uigetfile; 
[~,fn,fe] = fileparts(fn);

%{
ERPtables = cell(2, 4); % rows = unfilt vs filt; cols = channel
for r = 1:size(ERPtables,1)
    for c = 1:size(ERPtables,2)
        ERPtables{r,c} = table;
    end
end
clear r c 
%}

subjname = shortenFileName(fn)
sig = loadSignalObj(fullfile(fp,[fn,fe]));

%{
Fs = sig.SampleRate;

% highpass filtering (baseline removal) 
hpFilt = designfilt('highpassiir', ...
                    'StopbandFrequency', .5, ...
                    'PassbandFrequency', 1.5, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'butter');

% lowpass filtering (noise removal)
lpFilt = designfilt('lowpassiir', ...
                    'StopbandFrequency', 3000, ...
                    'PassbandFrequency', 2500, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');

%% stitching (questionable) 
sigHPF = filter(hpFilt, sig.Data_Unfiltered); % forward HPF 
sigBPF = filter(lpFilt, sigHPF); % forward BPF 
%}

%% get all trials 

N = size(sig.Data_Unfiltered,1) - size(sig.Data_LMS_LPF,1) + 1;

tBeforeTrig = .05;

[t_PrePost, d_PrePost, e_PrePost, tStim] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_HPF, sig.Data_LMS_LPF, sig.SampleRate, sig.Channels, N);

%{
[forward_t, forward_HPF, forward_BPF] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sigHPF, sigBPF, sig.SampleRate, sig.Channels, N);
%}

%% measurements: phase space 
tPost = t_PrePost(2,:);
PSvals = cell(2, length(e_PrePost)); % rows = unfilt vs filt; cols = channel

% calculate values for all channels 
for ch = 1:size(PSvals,2)

    % get filtered and unfiltered waveforms at this channel for each trial
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};
    tStim_ch = tStim{ch};

    % get number of trials 
    nTrl = size(d_PrePost_ch,1);

    % setup tables for calculated values 
    tblFilt = table('Size',[nTrl,5],...
                    'VariableTypes',["double","double","double","double","double"], ...
                    'VariableNames',{'PSA','PSAsub','PSA09','PSA30','PSA50'});
    tblFilt.StimTime = tStim_ch; % s
    tblUnfilt = tblFilt;

    % calculate values for each trial 
    for trl = 1:nTrl
        % populate row trl of table 

        % unfiltered 
        tblUnfilt.PSA(trl)    = getPhaseSpace(tPost, d_PrePost_ch(trl,:,2), [.005, .03], 10, false);
        tblUnfilt.PSAsub(trl) = getPhaseSpace(tPost, d_PrePost_ch(trl,:,2), [.005,.009], 10, false);
        tblUnfilt.PSA09(trl)  = getPhaseSpace(tPost, d_PrePost_ch(trl,:,2), [   0,.009], 10, false);
        tblUnfilt.PSA30(trl)  = getPhaseSpace(tPost, d_PrePost_ch(trl,:,2), [   0, .03], 10, false);
        tblUnfilt.PSA50(trl)  = getPhaseSpace(tPost, d_PrePost_ch(trl,:,2), [   0, .05], 10, false);

        % filtered  
        tblFilt.PSA(trl)    = getPhaseSpace(tPost, e_PrePost_ch(trl,:,2), [.005, .03], 10, false);
        tblFilt.PSAsub(trl) = getPhaseSpace(tPost, e_PrePost_ch(trl,:,2), [.005,.009], 10, false);
        tblFilt.PSA09(trl)  = getPhaseSpace(tPost, e_PrePost_ch(trl,:,2), [   0,.009], 10, false);
        tblFilt.PSA30(trl)  = getPhaseSpace(tPost, e_PrePost_ch(trl,:,2), [   0, .03], 10, false);
        tblFilt.PSA50(trl)  = getPhaseSpace(tPost, e_PrePost_ch(trl,:,2), [   0, .05], 10, false);
    end

    % store tables in cell array 
    PSvals{1, ch} = tblUnfilt; PSvals{2, ch} = tblFilt;
    clear tblFilt tblUnfilt e_PrePost_ch d_PrePost_ch;
end
clear tPost t_PrePost d_PrePost e_PrePost

%% save after loop 
saveFileName = [fp,filesep,fn,'_PSA_Spreadsheet.xlsx'];

rnames = {'Unfilt', 'Filt'};
for r = 1:size(PSvals,1)
    for c = 1:size(PSvals,2)
        sheetname = ['Channel ',sig.Channels(c).labels,' ',rnames{r}]
        psTbl = PSvals{r,c};
        %erpTbl = [erpTbl.Properties.RowNames, erpTbl];
        %erpTbl.Properties.VariableNames{1} = 'SubjName';
        writetable(psTbl,saveFileName,'Sheet',sheetname)
        clear psTbl sheetname
    end
end