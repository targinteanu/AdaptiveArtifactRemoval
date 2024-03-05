% Extract the SNR from a batch of _filtered_sigObj.mat files in a folder.
% Save the summary statistics as a spreadsheet. 
% Meant to work with rodent motor nerve data (Kiara). 
% Based on batchExtractEERP.m 

%% setup before loop 
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*_filtered_sigObj.mat');
cd(codedir);

SNRtables = cell(2, 4); % rows = unfilt vs filt; cols = channel
for r = 1:size(SNRtables,1)
    for c = 1:size(SNRtables,2)
        SNRtables{r,c} = table;
    end
end
clear r c 

%%
for f = 1:length(files)
%% get all trials 
subjname = files(f).name; subjname = subjname(1:(end-20))
sig = loadSignalObj([folder,filesep,files(f).name]);

N = size(sig.Data_Unfiltered,1) - size(sig.Data_LMS_LPF,1) + 1;

tBeforeTrig = .02;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_HPF, sig.Data_LMS_LPF, sig.SampleRate, sig.Channels, N);

%% SNR (signal : noise ratio) - measure ERP peaks 
tPost = t_PrePost(2,:);
SNRvals = cell(2, length(e_PrePost)); % rows = unfilt vs filt; cols = channel

% calculate values for all channels 
for ch = 1:size(SNRvals,2)

    % get filtered and unfiltered waveforms at this channel for each trial
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};

    % get number of trials 
    nTrl = size(d_PrePost_ch,1);

    % setup tables for calculated values 
    tblFilt = table('Size',[nTrl,6],...
                    'VariableTypes',["double","double","double","double","double","double"], ...
                    'VariableNames',{'n14amp','n14lat','p10amp','p10lat', 'mean','SD'});
    tblUnfilt = tblFilt;

    % calculate values for each trial 
    for trl = 1:nTrl

        %figure; subplot(211); 

        % populate row trl of table with n20 & n40 amps and latencies and stats (mean & SD)
        [n14,p10,stat] = measureERP(tPost, d_PrePost_ch(trl,:,2), .014, .01, [.008,.25], 10, false);
        tblUnfilt.n14amp(trl) = n14(1); tblUnfilt.n14lat(trl) = n14(2); % [amplitude; latency]
        tblUnfilt.p10amp(trl) = p10(1); tblUnfilt.p10lat(trl) = p10(2); % [amplitude; latency]
        tblUnfilt.mean(trl) = stat(1); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = only within selected time range
            % noise = std thru all times
        clear n14 p10 stat

        %title('unfilt');
        %subplot(212); 

        % do the same for filtered data 
        [n14,p10,stat] = measureERP(tPost, e_PrePost_ch(trl,:,2), .014, .01, [.008,.25], 10, false);
        tblFilt.n14amp(trl) = n14(1); tblFilt.n14lat(trl) = n14(2); % [amplitude; latency]
        tblFilt.p10amp(trl) = p10(1); tblFilt.p10lat(trl) = p10(2); % [amplitude; latency]
        tblFilt.mean(trl) = stat(1); tblFilt.SD(trl) = std(e_PrePost_ch(trl,:,2)); 
            % mean = only within selected time range
            % noise = std thru all times
        clear n14 p10 stat

        %title('filt');
    end

    % store tables in cell array 
    SNRvals{1, ch} = tblUnfilt; SNRvals{2, ch} = tblFilt;
    clear tblFilt tblUnfilt e_PrePost_ch d_PrePost_ch;
end
clear tPost t_PrePost d_PrePost e_PrePost

%% store SNR and ERP info in tables
for ch = 1:size(SNRvals,2)
    % get tables for this channel 
    tblUnfilt = SNRvals{1,ch}; tblFilt = SNRvals{2,ch}; 
    bigTblUnfilt = SNRtables{1,ch}; bigTblFilt = SNRtables{2,ch};


    bigTblUnfilt.n14ampMean(subjname) = mean(tblUnfilt.n14amp, 'omitnan');
    bigTblUnfilt.p10ampMean(subjname) = mean(tblUnfilt.p10amp, 'omitnan');
    bigTblUnfilt.n14latMean(subjname) = mean(tblUnfilt.n14lat, 'omitnan');
    bigTblUnfilt.p10latMean(subjname) = mean(tblUnfilt.p10lat, 'omitnan');
    bigTblUnfilt.snrMean(subjname) = mean( (tblUnfilt.n14amp-tblUnfilt.p10amp)./tblUnfilt.SD, 'omitnan' );
    bigTblUnfilt.noiseMean(subjname) = mean(tblUnfilt.SD); 
    
    bigTblUnfilt.n14ampSD(subjname) = std(tblUnfilt.n14amp, 'omitnan');
    bigTblUnfilt.p10ampSD(subjname) = std(tblUnfilt.p10amp, 'omitnan');
    bigTblUnfilt.n14latSD(subjname) = std(tblUnfilt.n14lat, 'omitnan');
    bigTblUnfilt.p10latSD(subjname) = std(tblUnfilt.p10lat, 'omitnan');
    bigTblUnfilt.snrSD(subjname) = std( (tblUnfilt.n14amp-tblUnfilt.p10amp)./tblUnfilt.SD, 'omitnan' );
    bigTblUnfilt.noiseSD(subjname) = std(tblUnfilt.SD); 

    bigTblUnfilt.n14num(subjname) = sum(~isnan(tblUnfilt.n14lat));
    bigTblUnfilt.p10num(subjname) = sum(~isnan(tblUnfilt.p10lat));


    bigTblFilt.n14ampMean(subjname) = mean(tblFilt.n14amp, 'omitnan');
    bigTblFilt.p10ampMean(subjname) = mean(tblFilt.p10amp, 'omitnan');
    bigTblFilt.n14latMean(subjname) = mean(tblFilt.n14lat, 'omitnan');
    bigTblFilt.p10latMean(subjname) = mean(tblFilt.p10lat, 'omitnan');
    bigTblFilt.snrMean(subjname) = mean( (tblFilt.n14amp-tblFilt.p10amp)./tblFilt.SD, 'omitnan' );
    bigTblFilt.noiseMean(subjname) = mean(tblFilt.SD); 
    
    bigTblFilt.n14ampSD(subjname) = std(tblFilt.n14amp, 'omitnan');
    bigTblFilt.p10ampSD(subjname) = std(tblFilt.p10amp, 'omitnan');
    bigTblFilt.n14latSD(subjname) = std(tblFilt.n14lat, 'omitnan');
    bigTblFilt.p10latSD(subjname) = std(tblFilt.p10lat, 'omitnan');
    bigTblFilt.snrSD(subjname) = std( (tblFilt.n14amp-tblFilt.p10amp)./tblFilt.SD, 'omitnan' );
    bigTblFilt.noiseSD(subjname) = std(tblFilt.SD); 

    bigTblFilt.n14num(subjname) = sum(~isnan(tblFilt.n14lat));
    bigTblFilt.p10num(subjname) = sum(~isnan(tblFilt.p10lat));


    SNRtables{1,ch} = bigTblUnfilt; SNRtables{2,ch} = bigTblFilt;
    clear bigTblFilt bigTblUnfilt tblFilt tblUnfilt
end

clear sig subjname
end

%% save after loop 
saveFileName = 'SNR_Spreadsheet.xlsx';
sig = loadSignalObj([folder,filesep,files(1).name]);

rnames = {'Unfilt', 'Filt'};
for r = 1:size(SNRtables,1)
    for c = 1:size(SNRtables,2)
        sheetname = ['Channel ',sig.Channels(c).labels,' ',rnames{r}]
        snrTbl = SNRtables{r,c};
        snrTbl = [snrTbl.Properties.RowNames, snrTbl];
        snrTbl.Properties.VariableNames{1} = 'SubjName';
        writetable(snrTbl,saveFileName,'Sheet',sheetname)
        clear snrTbl sheetname
    end
end