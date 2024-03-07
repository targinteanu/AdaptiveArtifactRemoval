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

tBeforeTrig = .025;
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
    tblFilt = table('Size',[nTrl,4],...
                    'VariableTypes',["double","double","double","double"], ...
                    'VariableNames',{'n13amp','n13lat', 'mean','SD'});
    tblUnfilt = tblFilt;

    % calculate values for each trial 
    for trl = 1:nTrl

        %figure; subplot(211); 

        % populate row trl of table with n20 & n40 amps and latencies and stats (mean & SD)
        [n13,~,~] = measureERP(tPost, d_PrePost_ch(trl,:,2), .013, [], [.008,.25], 10, false);
        tblUnfilt.n13amp(trl) = n13(1); tblUnfilt.n13lat(trl) = n13(2); % [amplitude; latency]
        tblUnfilt.mean(trl) = mean(d_PrePost_ch(trl,:,1)); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = prior to stim
            % noise = std thru all times
        clear n13 p10 stat

        %title('unfilt');
        %subplot(212); 

        % do the same for filtered data 
        [n13,~,~] = measureERP(tPost, e_PrePost_ch(trl,:,2), .013, [], [.008,.25], 10, false);
        tblFilt.n13amp(trl) = n13(1); tblFilt.n13lat(trl) = n13(2); % [amplitude; latency]
        tblFilt.mean(trl) = mean(e_PrePost_ch(trl,:,1)); tblFilt.SD(trl) = std(e_PrePost_ch(trl,:,2)); 
            % mean = prior to stim
            % noise = std thru all times
        clear n13 p10 stat

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


    bigTblUnfilt.n13ampMean(subjname) = mean(tblUnfilt.n13amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n13latMean(subjname) = mean(tblUnfilt.n13lat, 'omitnan');
    bigTblUnfilt.snrMean(subjname) = mean( (tblUnfilt.n13amp - tblUnfilt.mean)./tblUnfilt.SD, 'omitnan' );
    bigTblUnfilt.noiseMean(subjname) = mean(tblUnfilt.SD); 
    
    bigTblUnfilt.n13ampSD(subjname) = std(tblUnfilt.n13amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n13latSD(subjname) = std(tblUnfilt.n13lat, 'omitnan');
    bigTblUnfilt.snrSD(subjname) = std( (tblUnfilt.n13amp - tblUnfilt.mean)./tblUnfilt.SD, 'omitnan' );
    bigTblUnfilt.noiseSD(subjname) = std(tblUnfilt.SD); 

    bigTblUnfilt.n13num(subjname) = sum(~isnan(tblUnfilt.n13lat));


    bigTblFilt.n13ampMean(subjname) = mean(tblFilt.n13amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n13latMean(subjname) = mean(tblFilt.n13lat, 'omitnan');
    bigTblFilt.snrMean(subjname) = mean( (tblFilt.n13amp - tblFilt.mean)./tblFilt.SD, 'omitnan' );
    bigTblFilt.noiseMean(subjname) = mean(tblFilt.SD); 
    
    bigTblFilt.n13ampSD(subjname) = std(tblFilt.n13amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n13latSD(subjname) = std(tblFilt.n13lat, 'omitnan');
    bigTblFilt.snrSD(subjname) = std( (tblFilt.n13amp - tblFilt.mean)./tblFilt.SD, 'omitnan' );
    bigTblFilt.noiseSD(subjname) = std(tblFilt.SD); 

    bigTblFilt.n13num(subjname) = sum(~isnan(tblFilt.n13lat));


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