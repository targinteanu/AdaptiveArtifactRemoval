% Extract the ERP peak amplitudes/latencies and SNR from a  
% _filtered_sigObj.mat file.
% Save the statistics of each peak as a spreadsheet. 
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

%% get all trials 

N = size(sig.Data_Unfiltered,1) - size(sig.Data_LMS_LPF,1) + 1;

tBeforeTrig = .06;

[t_PrePost, d_PrePost, e_PrePost, tStim] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_HPF, sig.Data_LMS_LPF, sig.SampleRate, sig.Channels, N);

[forward_t, forward_HPF, forward_BPF] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sigHPF, sigBPF, sig.SampleRate, sig.Channels, N);

%% SNR (signal : noise ratio) - measure ERP peaks 
tPost = t_PrePost(2,:);
ERPvals = cell(2, length(e_PrePost)); % rows = unfilt vs filt; cols = channel

% calculate values for all channels 
for ch = 1:size(ERPvals,2)

    % get filtered and unfiltered waveforms at this channel for each trial
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};
    forward_HPF_ch = forward_HPF{ch};
    forward_BPF_ch = forward_BPF{ch};
    tStim_ch = tStim{ch};

    % get number of trials 
    nTrl = size(d_PrePost_ch,1);

    % setup tables for calculated values 
    tblFilt = table('Size',[nTrl,8],...
                    'VariableTypes',["double","double","double","double","double","double","double","double"], ...
                    'VariableNames',{'n07amp','n07lat','n15amp','n15lat','n40amp','n40lat', 'mean','SD'});
    tblFilt.StimTime = tStim_ch; % s
    tblUnfilt = tblFilt;

    % calculate values for each trial 
    for trl = 1:nTrl

        %figure; subplot(211); 

        % populate row trl of table with n20 & n40 amps and latencies and stats (mean & SD)
        [n20n40,~,stat] = measureERP(tPost, d_PrePost_ch(trl,:,2), [.007 .015 .04], [], [.004,.1], 100, false);
        tblUnfilt.n07amp(trl) = n20n40(1,1); tblUnfilt.n07lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblUnfilt.n15amp(trl) = n20n40(1,2); tblUnfilt.n15lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblUnfilt.n40amp(trl) = n20n40(1,3); tblUnfilt.n40lat(trl) = n20n40(2,3); % [amplitude; latency]
        tblUnfilt.mean(trl) = mean(forward_HPF_ch(trl,:,1)); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = prior to stim
            % noise = std thru all times
        clear n20n40 stat

        %title('unfilt');
        %subplot(212); 

        % do the same for filtered data 
        [n20n40,~,stat] = measureERP(tPost, e_PrePost_ch(trl,:,2), [.007 .015 .04], [], [.004,.1], 100, false);
        tblFilt.n07amp(trl) = n20n40(1,1); tblFilt.n07lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblFilt.n15amp(trl) = n20n40(1,2); tblFilt.n15lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblFilt.n40amp(trl) = n20n40(1,3); tblFilt.n40lat(trl) = n20n40(2,3); % [amplitude; latency]
        tblFilt.mean(trl) = mean(forward_BPF_ch(trl,:,1)); tblFilt.SD(trl) = std(e_PrePost_ch(trl,:,2)); 
        clear n20n40 stat

        %title('filt');
    end

    % subtract baseline from amp val 
    for V = {'n07amp', 'n15amp', 'n40amp'}
        v = V{:};
        eval(['tblUnfilt.',v,' = tblUnfilt.',v,' - tblUnfilt.mean;']);
        eval(['tblFilt.',v,' = tblFilt.',v,' - tblFilt.mean;']);
        clear v
    end
    clear V

    % store tables in cell array 
    ERPvals{1, ch} = tblUnfilt; ERPvals{2, ch} = tblFilt;
    clear tblFilt tblUnfilt e_PrePost_ch d_PrePost_ch;
end
clear tPost t_PrePost d_PrePost e_PrePost

%% save after loop 
saveFileName = [fp,filesep,fn,'_ERP_SNR_Spreadsheet.xlsx'];

rnames = {'Unfilt', 'Filt'};
for r = 1:size(ERPvals,1)
    for c = 1:size(ERPvals,2)
        sheetname = ['Channel ',sig.Channels(c).labels,' ',rnames{r}]
        erpTbl = ERPvals{r,c};
        %erpTbl = [erpTbl.Properties.RowNames, erpTbl];
        %erpTbl.Properties.VariableNames{1} = 'SubjName';
        writetable(erpTbl,saveFileName,'Sheet',sheetname)
        clear erpTbl sheetname
    end
end