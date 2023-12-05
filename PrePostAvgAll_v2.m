function [sgramUnfiltBefore, sgramUnfiltAfter, sgramFiltBefore, sgramFiltAfter, ...
    fig] = PrePostAvgAll_v2(tBeforeTrig, t_PrePost, d_PrePost, e_PrePost, Fs, uchan, nUpdates)
% Plot the filtered and unfiltered signals and their power spectra/spectrograms 
% before and after the stimulus, averaged over all available trials.  
% 
% Inputs: 
%   tBeforeTrig: time before/after stimulus to plot (s) 
%   t_PrePost: time with respect to stim as [before; after]
%   d_PrePost: unfiltered value as cell array for each channel; 
%              d_PrePost{i}(r,c,1) is before stim, trial r, timepoint c, channel i;
%              d_PrePost{i}(r,c,2) is after stim, trial r, timepoint c, channel i
%   e_PrePost: filtered value in the same format as d_PrePost 
%   Fs: sampling rate of g, d, and e_t (Hz) 
%   uchan: array of unique channels, same length as width of g, d, and e
%   nUpdates: how many times to display progress. 0 = no output 
% 
% Outputs: 
%   sgramUnfiltBefore: unfiltered spectrogram info before stim as cell array for each channel: 
%                      sgramUnfiltBefore{i,1} = channel i
%                      sgramUnfiltBefore{i,2} = frequency vals
%                      sgramUnfiltBefore{i,3} = time vals
%                      sgramUnfiltBefore{i,4} = spectrogram vals
%   sgramUnfiltAfter: unfiltered spectrogram info after stim in same format as above 
%   sgramFiltBefore: filtered spectrogram info before stim in same format as above 
%   sgramFiltAfter: filtered spectrogram info after stim in same format as above
%   fig: matlab figure, 5x2, showing time, spectrum, and spectrogram
%        before/after stim and with/without filtering, averaged over all
%        trials 


% TO DO: consider breaking spectrogram stuff into separate function or
% having option to disable it 

%% display: 
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

tlSz = 16;
lbSz = 16;
tkSz = 12;
lgSz = 14;

%% plotting averaged signals before and after stim 

for chIdx = 1:length(uchan)
    chan = uchan(chIdx);
    sigFiltCh = e_PrePost{chIdx};
    sigUnfiltCh = d_PrePost{chIdx};

    meanFiltBefore = mean(sigFiltCh(:,:,1),1);
    meanFiltAfter  = mean(sigFiltCh(:,:,2),1);
    errbFiltBefore =  std(sigFiltCh(:,:,1),[],1);
    errbFiltAfter  =  std(sigFiltCh(:,:,2),[],1);
    meanUnfiltBefore = mean(sigUnfiltCh(:,:,1),1);
    meanUnfiltAfter  = mean(sigUnfiltCh(:,:,2),1);
    errbUnfiltBefore =  std(sigUnfiltCh(:,:,1),[],1);
    errbUnfiltAfter  =  std(sigUnfiltCh(:,:,2),[],1);

    % calculate power spectra 
    [wFiltBefore, spectFiltBeforeCh] = PowerSpectrum(sigFiltCh(:,:,1), Fs);
    [wFiltAfter, spectFiltAfterCh] = PowerSpectrum(sigFiltCh(:,:,2), Fs);
    [wUnfiltBefore, spectUnfiltBeforeCh] = PowerSpectrum(sigUnfiltCh(:,:,1), Fs);
    [wUnfiltAfter, spectUnfiltAfterCh] = PowerSpectrum(sigUnfiltCh(:,:,2), Fs);

    %{
    % take magnitude of spectra (not needed if power spectra) 
    spectFiltBeforeCh = abs(spectFiltBeforeCh);
    spectFiltAfterCh = abs(spectFiltAfterCh);
    spectUnfiltBeforeCh = abs(spectUnfiltBeforeCh); 
    spectUnfiltAfterCh = abs(spectUnfiltAfterCh);
    %}

    % take mean/SD in time and frequency domains 
    meanSpectFiltBefore = mean(spectFiltBeforeCh,1); 
    errbSpectFiltBefore =  std(spectFiltBeforeCh,[],1);
    meanSpectFiltAfter = mean(spectFiltAfterCh,1); 
    errbSpectFiltAfter =  std(spectFiltAfterCh,[],1);
    meanSpectUnfiltBefore = mean(spectUnfiltBeforeCh,1); 
    errbSpectUnfiltBefore =  std(spectUnfiltBeforeCh,[],1);
    meanSpectUnfiltAfter = mean(spectUnfiltAfterCh,1); 
    errbSpectUnfiltAfter =  std(spectUnfiltAfterCh,[],1);

    fig(chIdx,1) = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig(chIdx,1)); 
    sgtitle(['Channel ',num2str(uchan(chIdx)),' Avg. Response to Stim'], ...
        'FontSize', tlSz);

    figure(fig(chIdx,1)); ax(1) = subplot(2,2,1); 
           plotWithDistrib(t_PrePost(1,:), meanUnfiltBefore, errbUnfiltBefore, ltRed);
    yrng = plotWithDistrib(t_PrePost(1,:), meanFiltBefore, errbFiltBefore, dkBlue);
    title('i) Pre-Stimulus Waveform', 'FontSize',tlSz); grid on; 
    xlabel('time (s)', 'FontSize',lbSz); ylabel('Signal (V)', 'FontSize',lbSz); 
    ylim(yrng(2,:)); 
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', ...
        'Location','eastoutside', 'FontSize',lgSz);
    set(gca, 'FontSize', tkSz);

    figure(fig(chIdx,1)); ax(2) = subplot(2,2,2); 
           plotWithDistrib(t_PrePost(2,:), meanUnfiltAfter, errbUnfiltAfter, ltRed);
    yrng = plotWithDistrib(t_PrePost(2,:), meanFiltAfter, errbFiltAfter, dkBlue);
    title('ii) Post-Stimulus Waveform', 'FontSize',tlSz); grid on; 
    xlabel('time (s)', 'FontSize',lbSz); ylabel('Signal (V)', 'FontSize',lbSz);
    ylim(yrng(2,:));
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', ...
        'Location','eastoutside', 'FontSize',lgSz);
    set(gca, 'FontSize', tkSz);

    linkaxes(ax,'y');
    clear ax;

    figure(fig(chIdx,1)); ax(1) = subplot(2,2,3); 
    %semilogy(wUnfiltAfter, meanSpectUnfiltAfter, 'Color', ltRed); hold on; 
    %semilogy(wFiltAfter, meanSpectFiltAfter, 'Color', ltBlue);
    plotWithDistrib(wUnfiltBefore, meanSpectUnfiltBefore, errbSpectUnfiltBefore, ltRed); hold on;
    plotWithDistrib(wFiltBefore, meanSpectFiltBefore, errbSpectFiltBefore, dkBlue);
    title('iii) Pre-Stimulus Spectrum', 'FontSize',tlSz); grid on; 
    set(gca, 'YScale', 'log');
    xlim([0 5000]);
    xlabel('Frequency (Hz)', 'FontSize',lbSz); ylabel('Power (V^2*s^2)', 'FontSize',lbSz);
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', ...
        'Location','eastoutside', 'FontSize',lgSz);
    set(gca, 'FontSize', tkSz);

    figure(fig(chIdx,1)); ax(2) = subplot(2,2,4); 
    %semilogy(wUnfiltBefore, meanSpectUnfiltBefore, 'Color', ltRed); hold on; 
    %semilogy(wFiltBefore, meanSpectFiltBefore, 'Color', ltBlue);
    plotWithDistrib(wUnfiltAfter, meanSpectUnfiltAfter, errbSpectUnfiltAfter, ltRed); hold on;
    plotWithDistrib(wFiltAfter, meanSpectFiltAfter, errbSpectFiltAfter, dkBlue);
    title('iv) Post-Stimulus Spectrum', 'FontSize',tlSz); grid on; 
    set(gca, 'YScale', 'log');
    xlim([0 5000]);
    xlabel('Frequency (Hz)', 'FontSize',lbSz); ylabel('Power (V^2*s^2)', 'FontSize',lbSz);
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', ...
        'Location','eastoutside', 'FontSize',lbSz);
    set(gca, 'FontSize', tkSz);

    linkaxes(ax,'y');
    clear ax;

    pause(.25);
end

end