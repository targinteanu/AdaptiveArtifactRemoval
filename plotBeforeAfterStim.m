function fig = plotBeforeAfterStim(tBeforeTrig, g, d, e_t, Fs, uchan, N, nTrlInBatch, nUpdates)
% Plot the filtered and unfiltered signals and their power spectra some
% time before and after the stimulus. 
% 
% Inputs: 
%   tBeforeTrig: time before/after stimulus to plot (s) 
%   g: stimulus value over time, as columns 
%   d: unfiltered value over time, as columns (V)
%   e_t: filtered value over time, as columns (V)
%   Fs: sampling rate of g, d, and e_t (Hz) 
%   uchan: array of unique channels, same length as width of g, d, and e_t 
%   N: number of filter taps         
% 
% Outputs: 
%   fig: matlab figure, 3x2, showing time and spectrum before/after stim
%        and with/without filtering 

%% getting signals before and after stim  
%tBeforeTrig = .29; % s
nBeforeTrig = floor(tBeforeTrig*Fs); % samples
tBeforeTrig = nBeforeTrig/Fs;
t_PrePost = [-tBeforeTrig:(1/Fs):0; 0:(1/Fs):tBeforeTrig]; % [before; after]

d_PrePost           = cell(1, length(uchan));
%d_lpf_PrePost       = cell(size(d_PrePost));
%e_train_PrePost     = cell(size(d_PrePost));
%e_train_lpf_PrePost = cell(size(d_PrePost));
%e_test_PrePost      = cell(size(d_PrePost));
%e_test_lpf_PrePost  = cell(size(d_PrePost));
e_t_PrePost         = cell(size(d_PrePost));
%e_t_lpf_PrePost     = cell(size(d_PrePost));

for chIdx = 1:length(uchan)
    gch = g(:,chIdx);
    trig = [0; abs(diff(gch))];
    trig = trig > .1*max(trig); trig = find(trig);

    d_PrePost_ch           = nan(length(trig), nBeforeTrig+1, 2);
%    d_lpf_PrePost_ch       = nan(size(d_PrePost_ch));
%    e_train_PrePost_ch     = nan(size(d_PrePost_ch));
%    e_train_lpf_PrePost_ch = nan(size(d_PrePost_ch));
%    e_test_PrePost_ch      = nan(size(d_PrePost_ch));
%    e_test_lpf_PrePost_ch  = nan(size(d_PrePost_ch));
    e_t_PrePost_ch         = nan(size(d_PrePost_ch));
%    e_t_lpf_PrePost_ch     = nan(size(d_PrePost_ch));

    for trIdx = 1:length(trig)
        tr = trig(trIdx); % timepoint

        % consider a more robust solution; padding e_t or making a separate
        % time vector
        % consider aligning so that rows of e correspond to d exactly 
        if (tr -nBeforeTrig > 0) & (tr +nBeforeTrig <= size(d,1))
            d_PrePost_ch(trIdx,:,1) = d(tr + ((-nBeforeTrig):0), chIdx);
            d_PrePost_ch(trIdx,:,2) = d(tr + (  0:nBeforeTrig ), chIdx);
%            d_lpf_PrePost_ch(trIdx,:,1) = d_lpf(tr + ((-nBeforeTrig):0), idx);
%            d_lpf_PrePost_ch(trIdx,:,2) = d_lpf(tr + (  0:nBeforeTrig ), idx);
        end
        if (tr -N+1 -nBeforeTrig > 0) & (tr -N+1 +nBeforeTrig <= size(e_t,1))
            e_t_PrePost_ch(trIdx,:,1) = e_t(tr -N+1 + ((-nBeforeTrig):0), chIdx);
            e_t_PrePost_ch(trIdx,:,2) = e_t(tr -N+1 + (  0:nBeforeTrig ), chIdx);
%            e_t_lpf_PrePost_ch(trIdx,:,1) = e_t_lpf(tr -N+1 + ((-nBeforeTrig):0), idx);
%            e_t_lpf_PrePost_ch(trIdx,:,2) = e_t_lpf(tr -N+1 + (  0:nBeforeTrig ), idx);
        end
    end

    toRemove = sum(e_t_PrePost_ch,2); 
    toRemove = sum(toRemove,3);
    toRemove = isnan(toRemove);
    e_t_PrePost_ch = e_t_PrePost_ch(~toRemove,:,:);
%    e_t_lpf_PrePost_ch = e_t_lpf_PrePost_ch(~toRemove,:,:);
    toRemove = sum(d_PrePost_ch,2); 
    toRemove = sum(toRemove,3);
    toRemove = isnan(toRemove);
    d_PrePost_ch = d_PrePost_ch(~toRemove,:,:);
%    d_lpf_PrePost_ch = d_lpf_PrePost_ch(~toRemove,:,:);

    d_PrePost{chIdx} = d_PrePost_ch;
%    d_lpf_PrePost{idx} = d_lpf_PrePost_ch;
%    e_train_PrePost{idx} = e_train_PrePost_ch;
%    e_train_lpf_PrePost{idx} = e_train_lpf_PrePost_ch;
%    e_test_PrePost{idx} = e_test_PrePost_ch;
%    e_test_lpf_PrePost{idx} = e_test_lpf_PrePost_ch;
    e_t_PrePost{chIdx} = e_t_PrePost_ch;
%    e_t_lpf_PrePost{idx} = e_t_lpf_PrePost_ch;
end

%% cleanup 
clear d_PrePost_ch d_lpf_PrePost_che_train_PrePost_ch e_train_lpf_PrePost_ch 
clear e_test_PrePost_ch e_test_lpf_PrePost_ch e_t_PrePost_ch e_t_lpf_PrePost_ch 
clear gch trig trIdx tr toRemove

%% colors: 
dkBlue  = [  1,  50, 130] /255;
ltBlue  = [145, 190, 255] /255;
dkRed   = [110,   0,   0] /255;
ltRed   = [250, 150, 150] /255;
dkBlack = [  0,   0,   0] /255;
ltBlack = [110, 110, 110] /255;
dkPink  = [130,   1, 125] /255;
ltPink  = [255, 120, 245] /255;
dkAqua  = [  0, 100, 120] /255;
dkGreen = [  0,  85,  15] /255;

%% plotting averaged signals before and after stim 

sgramWinLen = .005; % s
sgramWinLen = ceil(max(32, sgramWinLen/Fs)); % samples 
sgramFmax = Fs/2.25; 
sgramFWinLen = 4/(tBeforeTrig); % Hz
sgramFQ = 0:sgramFWinLen:sgramFmax; 

for chIdx = 1:length(uchan)
    sigFiltCh = e_t_PrePost{chIdx};
    sigUnfiltCh = d_PrePost{chIdx};

    meanFiltBefore = mean(sigFiltCh(:,:,1));
    meanFiltAfter  = mean(sigFiltCh(:,:,2));
    errbFiltBefore =  std(sigFiltCh(:,:,1));
    errbFiltAfter  =  std(sigFiltCh(:,:,2));
    meanUnfiltBefore = mean(sigUnfiltCh(:,:,1));
    meanUnfiltAfter  = mean(sigUnfiltCh(:,:,2));
    errbUnfiltBefore =  std(sigUnfiltCh(:,:,1));
    errbUnfiltAfter  =  std(sigUnfiltCh(:,:,2));

    [~, sFiltBeforeF, sFiltBeforeT, sFiltBefore] = ...
        spectrogram(meanFiltBefore, sgramWinLen, [], sgramFQ, Fs);
    [~, sFiltAfterF, sFiltAfterT, sFiltAfter] = ...
        spectrogram(meanFiltAfter, sgramWinLen, [], sgramFQ, Fs);
    [~, sUnfiltBeforeF, sUnfiltBeforeT, sUnfiltBefore] = ...
        spectrogram(meanUnfiltBefore, sgramWinLen, [], sgramFQ, Fs);
    [~, sUnfiltAfterF, sUnfiltAfterT, sUnfiltAfter] = ...
        spectrogram(meanUnfiltAfter, sgramWinLen, [], sgramFQ, Fs);
    sFiltBefore = zeros([size(sFiltBefore),size(sigFiltCh,1)]);
    sFiltAfter = zeros([size(sFiltAfter),size(sigFiltCh,1)]);
    sUnfiltBefore = zeros([size(sUnfiltBefore),size(sigUnfiltCh,1)]);
    sUnfiltAfter = zeros([size(sUnfiltAfter),size(sigUnfiltCh,1)]);
    for trl = 1:size(sigFiltCh,1)
        [~,~,~,sFiltBefore(:,:,trl)] = spectrogram(sigFiltCh(trl,:,1), sgramWinLen, [], sgramFQ, Fs);
        [~,~,~,sFiltAfter(:,:,trl)] = spectrogram(sigFiltCh(trl,:,2), sgramWinLen, [], sgramFQ, Fs);
        if nUpdates
            if ~mod(trl, floor(size(sigFiltCh,1)/(nUpdates)))
                disp(['Obtaining Channel ',num2str(uchan(chIdx)),' Spectrograms: ',...
                      num2str(100*trl/(size(sigFiltCh,1)+size(sigUnfiltCh,1))),'%']);
            end
        end
    end
    for trl = 1:size(sigUnfiltCh,1)
        [~,~,~,sUnfiltBefore(:,:,trl)] = spectrogram(sigUnfiltCh(trl,:,1), sgramWinLen, [], sgramFQ, Fs);
        [~,~,~,sUnfiltAfter(:,:,trl)] = spectrogram(sigUnfiltCh(trl,:,2), sgramWinLen, [], sgramFQ, Fs);
        if nUpdates
            if ~mod(trl, floor(size(sigUnfiltCh,1)/(nUpdates)))
                disp(['Obtaining Channel ',num2str(uchan(chIdx)),' Spectrograms: ',...
                      num2str(100*(size(sigFiltCh,1)+trl)/(size(sigFiltCh,1)+size(sigUnfiltCh,1))),'%']);
            end
        end
    end

    [wFiltBefore, spectFiltBeforeCh] = PowerSpectrum(sigFiltCh(:,:,1), Fs);
    [wFiltAfter, spectFiltAfterCh] = PowerSpectrum(sigFiltCh(:,:,2), Fs);
    [wUnfiltBefore, spectUnfiltBeforeCh] = PowerSpectrum(sigUnfiltCh(:,:,1), Fs);
    [wUnfiltAfter, spectUnfiltAfterCh] = PowerSpectrum(sigUnfiltCh(:,:,2), Fs);

    %{
    spectFiltBeforeCh = abs(spectFiltBeforeCh);
    spectFiltAfterCh = abs(spectFiltAfterCh);
    spectUnfiltBeforeCh = abs(spectUnfiltBeforeCh); 
    spectUnfiltAfterCh = abs(spectUnfiltAfterCh);
    %}

    meanSpectFiltBefore = mean(spectFiltBeforeCh); 
    errbSpectFiltBefore =  std(spectFiltBeforeCh);
    meanSpectFiltAfter = mean(spectFiltAfterCh); 
    errbSpectFiltAfter =  std(spectFiltAfterCh);
    meanSpectUnfiltBefore = mean(spectUnfiltBeforeCh); 
    errbSpectUnfiltBefore =  std(spectUnfiltBeforeCh);
    meanSpectUnfiltAfter = mean(spectUnfiltAfterCh); 
    errbSpectUnfiltAfter =  std(spectUnfiltAfterCh);

    fig(chIdx,1) = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig(chIdx,1)); sgtitle(['Channel ',num2str(uchan(chIdx)),' Avg. Response to Stim']);

    figure(fig(chIdx,1)); subplot(5,2,1); 
           plotWithDistrib(t_PrePost(1,:), meanUnfiltBefore, errbUnfiltBefore, ltRed);
    yrng = plotWithDistrib(t_PrePost(1,:), meanFiltBefore, errbFiltBefore, dkBlue);
    title('Filtered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 
    ylim(yrng(2,:)); 
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig(chIdx,1)); subplot(5,2,2); 
           plotWithDistrib(t_PrePost(2,:), meanUnfiltAfter, errbUnfiltAfter, ltRed);
    yrng = plotWithDistrib(t_PrePost(2,:), meanFiltAfter, errbFiltAfter, dkBlue);
    title('Filtered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)');
    ylim(yrng(2,:));
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig(chIdx,1)); subplot(5,2,3);  
           plotWithDistrib(t_PrePost(1,:), meanFiltBefore, errbFiltBefore, ltBlue);
    yrng = plotWithDistrib(t_PrePost(1,:), meanUnfiltBefore, errbUnfiltBefore, dkRed);
    title('Unfiltered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)');
    ylim(yrng(2,:));
    legend('Filtered', '-1SD', '+1SD', 'Unfiltered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig(chIdx,1)); subplot(5,2,4);  
           plotWithDistrib(t_PrePost(2,:), meanFiltAfter, errbFiltAfter, ltBlue);
    yrng = plotWithDistrib(t_PrePost(2,:), meanUnfiltAfter, errbUnfiltAfter, dkRed);
    title('Unfiltered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)');
    ylim(yrng(2,:));
    legend('Filtered', '-1SD', '+1SD', 'Unfiltered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig(chIdx,1)); subplot(5,2,5); 
    semilogy(wUnfiltAfter, meanSpectUnfiltAfter, 'Color', ltRed); hold on; 
    semilogy(wFiltAfter, meanSpectFiltAfter, 'Color', ltBlue);
    plotWithDistrib(wUnfiltBefore, meanSpectUnfiltBefore, errbSpectUnfiltBefore, dkRed);
    plotWithDistrib(wFiltBefore, meanSpectFiltBefore, errbSpectFiltBefore, dkBlue);
    title('Spectrum Before'); grid on; 
    %set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)');
    legend('Unfiltered After', 'Filtered After', ...
        'Unfiltered Before', '-1SD', '+1SD', 'Filtered Before', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig(chIdx,1)); subplot(5,2,6); 
    semilogy(wUnfiltBefore, meanSpectUnfiltBefore, 'Color', ltRed); hold on; 
    semilogy(wFiltBefore, meanSpectFiltBefore, 'Color', ltBlue);
    plotWithDistrib(wUnfiltAfter, meanSpectUnfiltAfter, errbSpectUnfiltAfter, dkRed);
    plotWithDistrib(wFiltAfter, meanSpectFiltAfter, errbSpectFiltAfter, dkBlue);
    title('Spectrum After'); grid on; 
    %set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)');
    legend('Unfiltered Before', 'Filtered Before', ...
        'Unfiltered After', '-1SD', '+1SD', 'Filtered After', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig(chIdx,1)); subplot(5,2,7);
    imagesc(sFiltBeforeT, sFiltBeforeF, 10*log10(mean(sFiltBefore,3))); colorbar;
    title('Spectrogram Filtered Before (dB/Hz)'); 
    xlabel('time (s)'); ylabel('Frequency (Hz)');

    figure(fig(chIdx,1)); subplot(5,2,8);
    imagesc(sFiltAfterT, sFiltAfterF, 10*log10(mean(sFiltAfter,3))); colorbar;
    title('Spectrogram Filtered After (dB/Hz)'); 
    xlabel('time (s)'); ylabel('Frequency (Hz)');

    figure(fig(chIdx,1)); subplot(5,2,9);
    imagesc(sUnfiltBeforeT, sUnfiltBeforeF, 10*log10(mean(sUnfiltBefore,3))); colorbar;
    title('Spectrogram Unfiltered Before (dB/Hz)'); 
    xlabel('time (s)'); ylabel('Frequency (Hz)');

    figure(fig(chIdx,1)); subplot(5,2,10);
    imagesc(sUnfiltAfterT, sUnfiltAfterF, 10*log10(mean(sUnfiltAfter,3))); colorbar;
    title('Spectrogram Unfiltered After (dB/Hz)'); 
    xlabel('time (s)'); ylabel('Frequency (Hz)');

    pause(.25);
end

%% cleanup 
clear sFiltBeforeF sFiltBeforeT sFiltBefore sUnfiltBeforeF sUnfiltBeforeT sUnfiltBefore 
clear yrng 
clear sgramFQ sgramWinLen sgramFmax sgramFWinLen

%% Plotting averaged - in batches 
for chIdx = 1:length(uchan)
    sigFiltCh = e_t_PrePost{chIdx};
    sigUnfiltCh = d_PrePost{chIdx};

    [wFiltBefore, spectFiltBeforeCh] = PowerSpectrum(sigFiltCh(:,:,1), Fs);
    [wFiltAfter, spectFiltAfterCh] = PowerSpectrum(sigFiltCh(:,:,2), Fs);
    [wUnfiltBefore, spectUnfiltBeforeCh] = PowerSpectrum(sigUnfiltCh(:,:,1), Fs);
    [wUnfiltAfter, spectUnfiltAfterCh] = PowerSpectrum(sigUnfiltCh(:,:,2), Fs);

    % Filtered: means for each batch -------------------------------------
    nTrl = size(sigFiltCh,1);
    nBatch = floor(nTrl/nTrlInBatch);

    meanFiltBefore = zeros(nBatch, size(sigFiltCh,2));
    errbFiltBefore = zeros(size(meanFiltBefore));
    meanFiltAfter = zeros(size(meanFiltBefore));
    errbFiltAfter = zeros(size(meanFiltAfter));
    %%{
    meanSpectFiltBefore = zeros(nBatch, size(spectFiltBeforeCh,2));
    errbSpectFiltBefore = zeros(size(meanSpectFiltBefore));
    meanSpectFiltAfter = zeros(size(meanSpectFiltBefore));
    errbSpectFiltAfter = zeros(size(meanSpectFiltAfter));
    %}

    for batchIdx = 1:nBatch
        trl1 = 1 + nTrlInBatch*(batchIdx-1);
        curBatchTrlIdx = (0:(nTrlInBatch-1)) + trl1;

        % before
        curBatch = sigFiltCh(curBatchTrlIdx, :, 1);  
        meanFiltBefore(batchIdx,:) = mean(curBatch); 
        errbFiltBefore(batchIdx,:) =  std(curBatch);

        % after
        curBatch = sigFiltCh(curBatchTrlIdx, :, 2); 
        meanFiltAfter(batchIdx,:) = mean(curBatch); 
        errbFiltAfter(batchIdx,:) =  std(curBatch);

        % Spect 
        %%{
        curBatch = spectFiltBeforeCh(curBatchTrlIdx,:);
        meanSpectFiltBefore(batchIdx,:) = mean(curBatch); 
        errbSpectFiltBefore(batchIdx,:) =  std(curBatch);
        curBatch = spectFiltAfterCh(curBatchTrlIdx,:);
        meanSpectFiltAfter(batchIdx,:) = mean(curBatch); 
        errbSpectFiltAfter(batchIdx,:) =  std(curBatch);
        %}
    end
    % ====================================================================

    % Unfiltered: means for each batch -----------------------------------
    nTrl = size(sigUnfiltCh,1);
    nBatch = floor(nTrl/nTrlInBatch);

    meanUnfiltBefore = zeros(nBatch, size(sigUnfiltCh,2));
    errbUnfiltBefore = zeros(size(meanUnfiltBefore));
    meanUnfiltAfter = zeros(size(meanUnfiltBefore));
    errbUnfiltAfter = zeros(size(meanUnfiltAfter));
    %%{
    meanSpectUnfiltBefore = zeros(nBatch, size(spectUnfiltBeforeCh,2));
    errbSpectUnfiltBefore = zeros(size(meanSpectUnfiltBefore));
    meanSpectUnfiltAfter = zeros(size(meanSpectUnfiltBefore));
    errbSpectUnfiltAfter = zeros(size(meanSpectUnfiltAfter));
    %}

    for batchIdx = 1:nBatch
        trl1 = 1 + nTrlInBatch*(batchIdx-1);
        curBatchTrlIdx = (0:(nTrlInBatch-1)) + trl1;

        % before
        curBatch = sigUnfiltCh(curBatchTrlIdx, :, 1);  
        meanUnfiltBefore(batchIdx,:) = mean(curBatch); 
        errbUnfiltBefore(batchIdx,:) =  std(curBatch);

        % after
        curBatch = sigUnfiltCh(curBatchTrlIdx, :, 2); 
        meanUnfiltAfter(batchIdx,:) = mean(curBatch); 
        errbUnfiltAfter(batchIdx,:) =  std(curBatch);

        % Spect 
        %%{
        curBatch = spectUnfiltBeforeCh(curBatchTrlIdx,:);
        meanSpectUnfiltBefore(batchIdx,:) = mean(curBatch); 
        errbSpectUnfiltBefore(batchIdx,:) =  std(curBatch);
        curBatch = spectUnfiltAfterCh(curBatchTrlIdx,:);
        meanSpectUnfiltAfter(batchIdx,:) = mean(curBatch); 
        errbSpectUnfiltAfter(batchIdx,:) =  std(curBatch);
        %}
    end
    % ====================================================================

    fig(chIdx,2) = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig(chIdx,2)); sgtitle(['Channel ',num2str(uchan(chIdx)),' Batch Responses to Stim']);

    % Filtered: Plots ----------------------------------------------------
    nTrl = size(sigFiltCh,1);
    nBatch = floor(nTrl/nTrlInBatch);

    figure(fig(chIdx,2)); subplot(4,2,1); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(1,:), meanFiltBefore(batchIdx,:), errbFiltBefore(batchIdx,:), batchColor); 
    end
    title('Filtered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,2)); subplot(4,2,2); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(2,:), meanFiltAfter(batchIdx,:), errbFiltAfter(batchIdx,:), batchColor); 
    end
    title('Filtered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,2)); subplot(4,2,5); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(wFiltBefore, meanSpectFiltBefore(batchIdx,:), errbSpectFiltBefore(batchIdx,:), batchColor); 
    end
    title('Filtered Before'); set(gca, 'YScale', 'log'); grid on; 
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)'); 

    figure(fig(chIdx,2)); subplot(4,2,6); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(wFiltAfter, meanSpectFiltAfter(batchIdx,:), errbSpectFiltAfter(batchIdx,:), batchColor); 
    end
    title('Filtered After'); set(gca, 'YScale', 'log'); grid on; 
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)'); 
    % ====================================================================

    % Unfiltered: Plots --------------------------------------------------
    nTrl = size(sigUnfiltCh,1);
    nBatch = floor(nTrl/nTrlInBatch);

    figure(fig(chIdx,2)); subplot(4,2,3); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(1,:), meanUnfiltBefore(batchIdx,:), errbUnfiltBefore(batchIdx,:), batchColor); 
    end
    title('Unfiltered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,2)); subplot(4,2,4); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(2,:), meanUnfiltAfter(batchIdx,:), errbUnfiltAfter(batchIdx,:), batchColor); 
    end
    title('Unfiltered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,2)); subplot(4,2,7); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(wUnfiltBefore, meanSpectUnfiltBefore(batchIdx,:), errbSpectUnfiltBefore(batchIdx,:), batchColor); 
    end
    title('Unfiltered Before'); set(gca, 'YScale', 'log'); grid on; 
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)'); 

    figure(fig(chIdx,2)); subplot(4,2,8); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(wUnfiltAfter, meanSpectUnfiltAfter(batchIdx,:), errbSpectUnfiltAfter(batchIdx,:), batchColor); 
    end
    title('Unfiltered After'); set(gca, 'YScale', 'log'); grid on; 
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)'); 
    % ====================================================================

    pause(.25);
end

%% cleanup 
clear curBatchTrlIdx batchIdx curBatch theta batchColor

%% Plotting Stats
for chIdx = 1:length(uchan)
    sigFiltCh = e_t_PrePost{chIdx};
    sigUnfiltCh = d_PrePost{chIdx};

    sigFiltBeforeCh = sigFiltCh(:,:,1);
    sigUnfiltBeforeCh = sigUnfiltCh(:,:,1);

    [wFiltBefore, spectFiltBeforeCh] = PowerSpectrum(sigFiltCh(:,:,1), Fs);
    [wFiltAfter, spectFiltAfterCh] = PowerSpectrum(sigFiltCh(:,:,2), Fs);
    [wUnfiltBefore, spectUnfiltBeforeCh] = PowerSpectrum(sigUnfiltCh(:,:,1), Fs);
    [wUnfiltAfter, spectUnfiltAfterCh] = PowerSpectrum(sigUnfiltCh(:,:,2), Fs);

    % normalize power spectra: 
    spectFiltBeforeCh = spectFiltBeforeCh./sum(spectFiltBeforeCh,2);
    spectFiltAfterCh = spectFiltAfterCh./sum(spectFiltAfterCh,2);
    spectUnfiltBeforeCh = spectUnfiltBeforeCh./sum(spectUnfiltBeforeCh,2);
    spectUnfiltAfterCh = spectUnfiltAfterCh./sum(spectUnfiltAfterCh,2);

    nTrl = size(sigFiltCh,1);
    statsTimeFilt = zeros(nTrl,2);
    statsFreqFilt = zeros(nTrl,4);
    for trl = 1:nTrl
        [~,statsTimeFilt(trl,1)] = kstest2(sigFiltCh(trl,:,1), sigFiltCh(trl,:,2));
        [~,statsTimeFilt(trl,2)] =  ttest2(sigFiltCh(trl,:,1), sigFiltCh(trl,:,2), 'Vartype', 'unequal');
        
        statsFreqFilt(trl,1)  = signrank(spectFiltBeforeCh(trl,:) , spectFiltAfterCh(trl,:));
        [~,statsFreqFilt(trl,2)] =  corr(spectFiltBeforeCh(trl,:)', spectFiltAfterCh(trl,:)', 'Type', 'Pearson');
        [~,statsFreqFilt(trl,3)] =  corr(spectFiltBeforeCh(trl,:)', spectFiltAfterCh(trl,:)', 'Type', 'Spearman');
        [~,statsFreqFilt(trl,4)] = ttest(spectFiltBeforeCh(trl,:) , spectFiltAfterCh(trl,:));
    end

    nTrl = size(sigUnfiltCh,1);
    statsTimeUnfilt = zeros(nTrl,2);
    statsFreqUnfilt = zeros(nTrl,4);
    for trl = 1:nTrl
        [~,statsTimeUnfilt(trl,1)] = kstest2(sigUnfiltCh(trl,:,1), sigUnfiltCh(trl,:,2));
        [~,statsTimeUnfilt(trl,2)] =  ttest2(sigUnfiltCh(trl,:,1), sigUnfiltCh(trl,:,2), 'Vartype', 'unequal');
        
        statsFreqUnfilt(trl,1)  = signrank(spectUnfiltBeforeCh(trl,:) , spectUnfiltAfterCh(trl,:));
        [~,statsFreqUnfilt(trl,2)] =  corr(spectUnfiltBeforeCh(trl,:)', spectUnfiltAfterCh(trl,:)', 'Type', 'Pearson');
        [~,statsFreqUnfilt(trl,3)] =  corr(spectUnfiltBeforeCh(trl,:)', spectUnfiltAfterCh(trl,:)', 'Type', 'Spearman');
        [~,statsFreqUnfilt(trl,4)] = ttest(spectUnfiltBeforeCh(trl,:) , spectUnfiltAfterCh(trl,:));
    end

    statsTimeOverTimeFilt = zeros(2, size(t_PrePost,2));
    for tIdx = 1:size(t_PrePost,2)
%        [~,statsTimeOverTimeFilt(1,tIdx)] = kstest2(sigFiltBeforeCh(:), sigFiltCh(:,tIdx,2));
        [~,statsTimeOverTimeFilt(2,tIdx)] =  ttest2(sigFiltBeforeCh(:), sigFiltCh(:,tIdx,2), 'Vartype', 'unequal');
        if nUpdates
            if ~mod(tIdx, floor(size(t_PrePost,2)/(nUpdates)))
                disp(['Hypothesis-Testing Channel ',num2str(uchan(chIdx)),': ',...
                      num2str(50*tIdx/(size(t_PrePost,2))),'%']);
            end
        end
    end

    statsTimeOverTimeUnfilt = zeros(2, size(t_PrePost,2));
    for tIdx = 1:size(t_PrePost,2)
%        [~,statsTimeOverTimeUnfilt(1,tIdx)] = kstest2(sigUnfiltBeforeCh(:), sigUnfiltCh(:,tIdx,2));
        [~,statsTimeOverTimeUnfilt(2,tIdx)] =  ttest2(sigUnfiltBeforeCh(:), sigUnfiltCh(:,tIdx,2), 'Vartype', 'unequal');
        if nUpdates
            if ~mod(tIdx, floor(size(t_PrePost,2)/(nUpdates)))
                disp(['Hypothesis-Testing Channel ',num2str(uchan(chIdx)),': ',...
                      num2str(50 + 50*tIdx/(size(t_PrePost,2))),'%']);
            end
        end
    end

    fig(chIdx,3) = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig(chIdx,3)); sgtitle(['Channel ',num2str(uchan(chIdx)),' Before-After Comparison Statistics']);

    figure(fig(chIdx,3)); subplot(3,2,1);
    plot(statsTimeUnfilt(:,1));
    title('Unfiltered Time Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after'); 
    legend('KS Test', 'T Test');

    figure(fig(chIdx,3)); subplot(3,2,2);
    plot(statsTimeFilt(:,1));
    title('Filtered Time Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after');
    legend('KS Test', 'T Test');

    figure(fig(chIdx,3)); subplot(3,2,3);
    plot(statsFreqUnfilt(:,1)); 
    title('Unfiltered Frequency Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after');
    legend('Wilcoxon Rank', 'Pearson Corr.', 'Spearman Corr.', 'Paired T');

    figure(fig(chIdx,3)); subplot(3,2,4);
    plot(statsFreqFilt(:,1)); 
    title('Filtered Frequency Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after');
    legend('Wilcoxon Rank', 'Pearson Corr.', 'Spearman Corr.', 'Paired T');

    figure(fig(chIdx,3)); subplot(3,2,5);
    plot(t_PrePost(2,:), statsTimeOverTimeUnfilt);
    title('Unfiltered Time Domain'); grid on; 
    xlabel('time (s)'); ylabel('p value: before vs after'); 
    legend('KS Test', 'T Test');

    figure(fig(chIdx,3)); subplot(3,2,6);
    plot(t_PrePost(2,:), statsTimeOverTimeFilt);
    title('Filtered Time Domain'); grid on; 
    xlabel('time (s)'); ylabel('p value: before vs after'); 
    legend('KS Test', 'T Test');

    pause(.25);
end

end