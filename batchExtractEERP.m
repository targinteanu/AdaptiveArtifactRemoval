%% setup before loop 
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*_filtered_sigObj.mat');
cd(codedir);

tblNcol = 19;
%{
tblBlank = table('Size',[length(files),tblNcol], ...
                 'VariableTypes',repmat("double",[1,tblNcol]), ...
                 'RowNames',{files.name});
%}
ERPtables = cell(2, 4); % rows = unfilt vs filt; cols = channel
for r = 1:size(ERPtables,1)
    for c = 1:size(ERPtables,2)
        %ERPtables{r,c} = tblBlank;
        ERPtables{r,c} = table;
    end
end
clear r c 

%%
for f = 1:length(files)
%% get all trials 
subjname = files(f).name; subjname = subjname(1:(end-37))
sig = loadSignalObj([folder,'\',files(f).name]);

N = size(sig.Data_Unfiltered,1) - size(sig.Data_LMS_LPF,1) + 1;

tBeforeTrig = .06;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_Unfiltered, sig.Data_LMS_LPF, sig.SampleRate, sig.Channels, N);
%PrePostAvgAll_v2(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%% SNR (signal : noise ratio) - measure ERP peaks 
tPost = t_PrePost(2,:);
ERPvals = cell(2, length(e_PrePost)); % rows = unfilt vs filt; cols = channel

% calculate values for all channels 
for ch = 1:size(ERPvals,2)

    % get filtered and unfiltered waveforms at this channel for each trial
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};

    % get number of trials 
    nTrl = size(d_PrePost_ch,1);

    % setup tables for calculated values 
    tblFilt = table('Size',[nTrl,6],...
                    'VariableTypes',["double","double","double","double","double","double"], ...
                    'VariableNames',{'n20amp','n20lat','n40amp','n40lat', 'mean','SD'});
    tblUnfilt = tblFilt;

    % calculate values for each trial 
    for trl = 1:nTrl

        %figure; subplot(211); 

        % populate row trl of table with n20 & n40 amps and latencies and stats (mean & SD)
        [n20n40,~,stat] = measureERP(tPost, d_PrePost_ch(trl,:,2), [.01 .02 .04], [], [.008,.1], 200, false);
        tblUnfilt.n10amp(trl) = n20n40(1,1); tblUnfilt.n10lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblUnfilt.n20amp(trl) = n20n40(1,2); tblUnfilt.n20lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblUnfilt.n40amp(trl) = n20n40(1,3); tblUnfilt.n40lat(trl) = n20n40(2,3); % [amplitude; latency]
        tblUnfilt.mean(trl) = stat(1); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = only within selected time range
            % noise = std thru all times
        clear n20n40 stat

        %title('unfilt');
        %subplot(212); 

        % do the same for filtered data 
        [n20n40,~,stat] = measureERP(tPost, e_PrePost_ch(trl,:,2), [.01 .02 .04], [], [.008,.1], 200, false);
        tblFilt.n10amp(trl) = n20n40(1,1); tblFilt.n10lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblFilt.n20amp(trl) = n20n40(1,2); tblFilt.n20lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblFilt.n40amp(trl) = n20n40(1,3); tblFilt.n40lat(trl) = n20n40(2,3); % [amplitude; latency]
        tblFilt.mean(trl) = stat(1); tblFilt.SD(trl) = std(e_PrePost_ch(trl,:,2)); 
        clear n20n40 stat

        %title('filt');
    end

    % store tables in cell array 
    ERPvals{1, ch} = tblUnfilt; ERPvals{2, ch} = tblFilt;
    clear tblFilt tblUnfilt e_PrePost_ch d_PrePost_ch;
end
clear tPost t_PrePost d_PrePost e_PrePost

%% store SNR and ERP info in tables
for ch = 1:size(ERPvals,2)
    % get tables for this channel 
    tblUnfilt = ERPvals{1,ch}; tblFilt = ERPvals{2,ch}; 
    bigTblUnfilt = ERPtables{1,ch}; bigTblFilt = ERPtables{2,ch};


    bigTblUnfilt.n10ampMean(subjname) = mean(tblUnfilt.n10amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n20ampMean(subjname) = mean(tblUnfilt.n20amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n40ampMean(subjname) = mean(tblUnfilt.n40amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n10latMean(subjname) = mean(tblUnfilt.n10lat, 'omitnan');
    bigTblUnfilt.n20latMean(subjname) = mean(tblUnfilt.n20lat, 'omitnan');
    bigTblUnfilt.n40latMean(subjname) = mean(tblUnfilt.n40lat, 'omitnan');
    bigTblUnfilt.snrMean(subjname) = mean( (tblUnfilt.n20amp-tblUnfilt.mean)./tblUnfilt.SD, 'omitnan' );
    bigTblUnfilt.noiseMean(subjname) = mean(tblUnfilt.SD); 
    
    bigTblUnfilt.n10ampSD(subjname) = std(tblUnfilt.n10amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n20ampSD(subjname) = std(tblUnfilt.n20amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n40ampSD(subjname) = std(tblUnfilt.n40amp - tblUnfilt.mean, 'omitnan');
    bigTblUnfilt.n10latSD(subjname) = std(tblUnfilt.n10lat, 'omitnan');
    bigTblUnfilt.n20latSD(subjname) = std(tblUnfilt.n20lat, 'omitnan');
    bigTblUnfilt.n40latSD(subjname) = std(tblUnfilt.n40lat, 'omitnan');
    bigTblUnfilt.snrSD(subjname) = std( (tblUnfilt.n20amp-tblUnfilt.mean)./tblUnfilt.SD, 'omitnan' );
    bigTblUnfilt.noiseSD(subjname) = std(tblUnfilt.SD); 

    bigTblUnfilt.n10num(subjname) = sum(~isnan(tblUnfilt.n10lat));
    bigTblUnfilt.n20num(subjname) = sum(~isnan(tblUnfilt.n20lat));
    bigTblUnfilt.n40num(subjname) = sum(~isnan(tblUnfilt.n40lat));


    bigTblFilt.n10ampMean(subjname) = mean(tblFilt.n10amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n20ampMean(subjname) = mean(tblFilt.n20amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n40ampMean(subjname) = mean(tblFilt.n40amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n10latMean(subjname) = mean(tblFilt.n10lat, 'omitnan');
    bigTblFilt.n20latMean(subjname) = mean(tblFilt.n20lat, 'omitnan');
    bigTblFilt.n40latMean(subjname) = mean(tblFilt.n40lat, 'omitnan');
    bigTblFilt.snrMean(subjname) = mean( (tblFilt.n20amp-tblFilt.mean)./tblFilt.SD, 'omitnan' );
    bigTblFilt.noiseMean(subjname) = mean(tblFilt.SD); 
    
    bigTblFilt.n10ampSD(subjname) = std(tblFilt.n10amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n20ampSD(subjname) = std(tblFilt.n20amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n40ampSD(subjname) = std(tblFilt.n40amp - tblFilt.mean, 'omitnan');
    bigTblFilt.n10latSD(subjname) = std(tblFilt.n10lat, 'omitnan');
    bigTblFilt.n20latSD(subjname) = std(tblFilt.n20lat, 'omitnan');
    bigTblFilt.n40latSD(subjname) = std(tblFilt.n40lat, 'omitnan');
    bigTblFilt.snrSD(subjname) = std( (tblFilt.n20amp-tblFilt.mean)./tblFilt.SD, 'omitnan' );
    bigTblFilt.noiseSD(subjname) = std(tblFilt.SD); 

    bigTblFilt.n10num(subjname) = sum(~isnan(tblFilt.n10lat));
    bigTblFilt.n20num(subjname) = sum(~isnan(tblFilt.n20lat));
    bigTblFilt.n40num(subjname) = sum(~isnan(tblFilt.n40lat));

    ERPtables{1,ch} = bigTblUnfilt; ERPtables{2,ch} = bigTblFilt;
    clear bigTblFilt bigTblUnfilt tblFilt tblUnfilt
end

clear sig subjname
end

%% save after loop 
saveFileName = 'ERP_SNR_Spreadsheet.xlsx';
sig = loadSignalObj([folder,'\',files(1).name]);

rnames = {'Unfilt', 'Filt'};
for r = 1:size(ERPtables,1)
    for c = 1:size(ERPtables,2)
        sheetname = ['Channel ',sig.Channels(c).labels,' ',rnames{r}]
        erpTbl = ERPtables{r,c};
        erpTbl = [erpTbl.Properties.RowNames, erpTbl];
        erpTbl.Properties.VariableNames{1} = 'SubjName';
        writetable(erpTbl,saveFileName,'Sheet',sheetname)
        clear erpTbl sheetname
    end
end