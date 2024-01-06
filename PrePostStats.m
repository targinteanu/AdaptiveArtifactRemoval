function [fig] = PrePostStats(t_PrePost, d_PrePost, e_PrePost, Fs, uchan, bndlim, nUpdates)
% Plot the statistical difference between the signals and their spectra
% before and after the stimulus. 
% 
% Inputs: 
%   t_PrePost: time with respect to stim as [before; after]
%   d_PrePost: unfiltered value as cell array for each channel; 
%              d_PrePost{i}(r,c,1) is before stim, trial r, timepoint c, channel i;
%              d_PrePost{i}(r,c,2) is after stim, trial r, timepoint c, channel i
%   e_PrePost: filtered value in the same format as d_PrePost 
%   Fs: sampling rate of g, d, and e_t (Hz) 
%   uchan: array of unique channels, same length as width of g, d, and e
%   bndlim: frequency range to test: [lower, upper] limit (Hz) 
%   nUpdates: how many times to display progress. 0 = no output 
% 
% Outputs: 
%   fig: matlab figure, 3x2, showing with/without filtering: statistical
%        difference of each trial pair in the time domain, statistical
%        difference of each trial pair in the frequency domain, and
%        statistical difference of each timepoint after stimulus with the
%        signal before stimulus 


% TO DO: more outputs including the stats! as cell arrays? 

%% Plotting Stats
for chIdx = 1:length(uchan)
    chan = uchan(chIdx);
    sigFiltCh = e_PrePost{chIdx};
    sigUnfiltCh = d_PrePost{chIdx};

    sigFiltBeforeCh = sigFiltCh(:,:,1);
    sigUnfiltBeforeCh = sigUnfiltCh(:,:,1);

    [wFiltBefore, spectFiltBeforeCh] = PowerSpectrum(sigFiltCh(:,:,1), Fs);
    [wFiltAfter, spectFiltAfterCh] = PowerSpectrum(sigFiltCh(:,:,2), Fs);
    [wUnfiltBefore, spectUnfiltBeforeCh] = PowerSpectrum(sigUnfiltCh(:,:,1), Fs);
    [wUnfiltAfter, spectUnfiltAfterCh] = PowerSpectrum(sigUnfiltCh(:,:,2), Fs);

    if ~isempty(bndlim)
        % band-restrict power spectra:
        wIdx = (wFiltBefore >= bndlim(1)) & (wFiltBefore <= bndlim(2));
        spectFiltBeforeCh = spectFiltBeforeCh(:,wIdx);
        wIdx = (wFiltAfter >= bndlim(1)) & (wFiltAfter <= bndlim(2));
        spectFiltAfterCh = spectFiltAfterCh(:,wIdx);
        wIdx = (wUnfiltBefore >= bndlim(1)) & (wUnfiltBefore <= bndlim(2));
        spectUnfiltBeforeCh = spectUnfiltBeforeCh(:,wIdx);
        wIdx = (wUnfiltAfter >= bndlim(1)) & (wUnfiltAfter <= bndlim(2));
        spectUnfiltAfterCh = spectUnfiltAfterCh(:,wIdx);
    end

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
        baselineToTest = sigFiltBeforeCh(:,tIdx); % faster 
%        baselineToTest = sigFiltBeforeCh(:);      % probably more accurate 
        [~,statsTimeOverTimeFilt(1,tIdx)] = kstest2(baselineToTest, sigFiltCh(:,tIdx,2));
        [~,statsTimeOverTimeFilt(2,tIdx)] =  ttest2(baselineToTest, sigFiltCh(:,tIdx,2), 'Vartype', 'unequal');
        if nUpdates
            if ~mod(tIdx, floor(size(t_PrePost,2)/(nUpdates)))
                disp(['Hypothesis-Testing Channel ',chan.labels,': ',...
                      num2str(50*tIdx/(size(t_PrePost,2))),'%']);
            end
        end
    end

    statsTimeOverTimeUnfilt = zeros(2, size(t_PrePost,2));
    for tIdx = 1:size(t_PrePost,2)
        baselineToTest = sigUnfiltBeforeCh(:,tIdx); % faster 
%        baselineToTest = sigUnfiltBeforeCh(:);      % probably more accurate
        [~,statsTimeOverTimeUnfilt(1,tIdx)] = kstest2(baselineToTest, sigUnfiltCh(:,tIdx,2));
        [~,statsTimeOverTimeUnfilt(2,tIdx)] =  ttest2(baselineToTest, sigUnfiltCh(:,tIdx,2), 'Vartype', 'unequal');
        if nUpdates
            if ~mod(tIdx, floor(size(t_PrePost,2)/(nUpdates)))
                disp(['Hypothesis-Testing Channel ',chan.labels,': ',...
                      num2str(50 + 50*tIdx/(size(t_PrePost,2))),'%']);
            end
        end
    end

    fig(chIdx,1) = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig(chIdx,1)); sgtitle(['Channel ',chan.labels,' Before-After Comparison Statistics']);

    figure(fig(chIdx,1)); subplot(3,2,1);
    plot(statsTimeUnfilt(:,1:2));
    title('Unfiltered Time Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after'); 
    legend('KS Test', 'T Test');

    figure(fig(chIdx,1)); subplot(3,2,2);
    plot(statsTimeFilt(:,1:2));
    title('Filtered Time Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after');
    legend('KS Test', 'T Test');

    figure(fig(chIdx,1)); subplot(3,2,3);
    plot(statsFreqUnfilt(:,1)); 
    title('Unfiltered Frequency Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after');
    legend('Wilcoxon Rank', 'Pearson Corr.', 'Spearman Corr.', 'Paired T');

    figure(fig(chIdx,1)); subplot(3,2,4);
    plot(statsFreqFilt(:,1)); 
    title('Filtered Frequency Domain'); grid on; 
    xlabel('Trial #'); ylabel('p value: before vs after');
    legend('Wilcoxon Rank', 'Pearson Corr.', 'Spearman Corr.', 'Paired T');

    figure(fig(chIdx,1)); subplot(3,2,5);
    plot(t_PrePost(2,:), statsTimeOverTimeUnfilt);
    title('Unfiltered Time Domain'); grid on; 
    xlabel('time (s)'); ylabel('p value: before vs after'); 
    legend('KS Test', 'T Test');

    figure(fig(chIdx,1)); subplot(3,2,6);
    plot(t_PrePost(2,:), statsTimeOverTimeFilt);
    title('Filtered Time Domain'); grid on; 
    xlabel('time (s)'); ylabel('p value: before vs after'); 
    legend('KS Test', 'T Test');

    pause(.25);
end

end