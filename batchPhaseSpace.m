%% setup before loop 
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*_filtered_sigObj.mat');
cd(codedir);

PStables = cell(2, 4); % rows = unfilt vs filt; cols = channel
for r = 1:size(PStables,1)
    for c = 1:size(PStables,2)
        PStables{r,c} = table;
    end
end
clear r c 

%%
for f = 1:length(files)
%% get all trials 
subjname = files(f).name; subjname = shortenFileName(subjname)
sig = loadSignalObj([folder,filesep,files(f).name]);

N = size(sig.Data_Unfiltered,1) - size(sig.Data_LMS_LPF,1) + 1;

tBeforeTrig = .05;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_HPF, sig.Data_LMS_LPF, sig.SampleRate, sig.Channels, N);
%PrePostAvgAll_v2(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%% measurements: phase space
tPost = t_PrePost(2,:);
PSvals = cell(2, length(e_PrePost)); % rows = unfilt vs filt; cols = channel

% calculate values for all channels 
for ch = 1:size(PSvals,2)

    % get filtered and unfiltered waveforms at this channel for each trial
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};

    % get number of trials 
    nTrl = size(d_PrePost_ch,1);

    % setup tables for calculated values 
    tblFilt = table('Size',[nTrl,5],...
                    'VariableTypes',["double","double","double","double","double"], ...
                    'VariableNames',{'PSA','PSAsub','PSA09','PSA30','PSA50'});
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

%% store SNR and ERP info in tables
for ch = 1:size(PSvals,2)
    % get tables for this channel 
    tblUnfilt = PSvals{1,ch}; tblFilt = PSvals{2,ch}; 
    bigTblUnfilt = PStables{1,ch}; bigTblFilt = PStables{2,ch};

    bigTblUnfilt.PSAmean(subjname) = mean(tblUnfilt.PSA, 'omitnan');
    bigTblUnfilt.PSAsubmean(subjname) = mean(tblUnfilt.PSAsub, 'omitnan');
    bigTblUnfilt.PSA09mean(subjname) = mean(tblUnfilt.PSA09, 'omitnan');
    bigTblUnfilt.PSA30mean(subjname) = mean(tblUnfilt.PSA30, 'omitnan');
    bigTblUnfilt.PSA50mean(subjname) = mean(tblUnfilt.PSA50, 'omitnan');

    bigTblFilt.PSAmean(subjname) = mean(tblFilt.PSA, 'omitnan');
    bigTblFilt.PSAsubmean(subjname) = mean(tblFilt.PSAsub, 'omitnan');
    bigTblFilt.PSA09mean(subjname) = mean(tblFilt.PSA09, 'omitnan');
    bigTblFilt.PSA30mean(subjname) = mean(tblFilt.PSA30, 'omitnan');
    bigTblFilt.PSA50mean(subjname) = mean(tblFilt.PSA50, 'omitnan');

    bigTblUnfilt.PSAstd(subjname) = std(tblUnfilt.PSA, 'omitnan');
    bigTblUnfilt.PSAsubstd(subjname) = std(tblUnfilt.PSAsub, 'omitnan');
    bigTblUnfilt.PSA09std(subjname) = std(tblUnfilt.PSA09, 'omitnan');
    bigTblUnfilt.PSA30std(subjname) = std(tblUnfilt.PSA30, 'omitnan');
    bigTblUnfilt.PSA50std(subjname) = std(tblUnfilt.PSA50, 'omitnan');

    bigTblFilt.PSAstd(subjname) = std(tblFilt.PSA, 'omitnan');
    bigTblFilt.PSAsubstd(subjname) = std(tblFilt.PSAsub, 'omitnan');
    bigTblFilt.PSA09std(subjname) = std(tblFilt.PSA09, 'omitnan');
    bigTblFilt.PSA30std(subjname) = std(tblFilt.PSA30, 'omitnan');
    bigTblFilt.PSA50std(subjname) = std(tblFilt.PSA50, 'omitnan');

    bigTblUnfilt.PSAnum(subjname) = sum(~isnan(tblUnfilt.PSA));
    bigTblUnfilt.PSAsubnum(subjname) = sum(~isnan(tblUnfilt.PSAsub));
    bigTblUnfilt.PSA09num(subjname) = sum(~isnan(tblUnfilt.PSA09));
    bigTblUnfilt.PSA30num(subjname) = sum(~isnan(tblUnfilt.PSA30));
    bigTblUnfilt.PSA50num(subjname) = sum(~isnan(tblUnfilt.PSA50));

    bigTblFilt.PSAnum(subjname) = sum(~isnan(tblFilt.PSA));
    bigTblFilt.PSAsubnum(subjname) = sum(~isnan(tblFilt.PSAsub));
    bigTblFilt.PSA09num(subjname) = sum(~isnan(tblFilt.PSA09));
    bigTblFilt.PSA30num(subjname) = sum(~isnan(tblFilt.PSA30));
    bigTblFilt.PSA50num(subjname) = sum(~isnan(tblFilt.PSA50));

    PStables{1,ch} = bigTblUnfilt; PStables{2,ch} = bigTblFilt;
    clear bigTblFilt bigTblUnfilt tblFilt tblUnfilt
end

clear sig subjname
end

%% save after loop 
saveFileName = 'PSA_Spreadsheet.xlsx';
sig = loadSignalObj([folder,filesep,files(1).name]);

rnames = {'Unfilt', 'Filt'};
for r = 1:size(PStables,1)
    for c = 1:size(PStables,2)
        sheetname = ['Channel ',sig.Channels(c).labels,' ',rnames{r}]
        psTbl = PStables{r,c};
        psTbl = [psTbl.Properties.RowNames, psTbl];
        psTbl.Properties.VariableNames{1} = 'SubjName';
        writetable(psTbl,saveFileName,'Sheet',sheetname)
        clear psTbl sheetname
    end
end