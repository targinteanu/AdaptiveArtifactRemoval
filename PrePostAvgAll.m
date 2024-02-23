function [sgramUnfiltBefore, sgramUnfiltAfter, sgramFiltBefore, sgramFiltAfter, ...
    fig] = PrePostAvgAll(tBeforeTrig, t_PrePost, d_PrePost, e_PrePost, Fs, uchan, nUpdates)
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

% spectrogram parameters ------------------------------ make input to function? 
sgramWinLen = .0001; % s
%sgramWinLen = ceil(max(32, sgramWinLen/Fs)); % samples 
sgramWinLen = 10; 
sgramFmax = Fs/2.25; % -------------------------------- also use for power spectrum plots?
sgramFWinLen = .1/(tBeforeTrig); % Hz
sgramFQ = 0:sgramFWinLen:sgramFmax; 

% spectrogram outputs --------------------------------- can be organized better with multi-d matrices?
sgramFiltBefore   = cell(length(uchan), 4);
sgramFiltAfter    = cell(size(sgramFiltBefore));
sgramUnfiltBefore = cell(size(sgramFiltBefore));
sgramUnfiltAfter  = cell(size(sgramUnfiltBefore));

for chIdx = 1:length(uchan)
    chan = uchan(chIdx);
    sigFiltCh = e_PrePost{chIdx};
    sigUnfiltCh = d_PrePost{chIdx};

    meanFiltBefore = mean(sigFiltCh(:,:,1));
    meanFiltAfter  = mean(sigFiltCh(:,:,2));
    errbFiltBefore =  std(sigFiltCh(:,:,1));
    errbFiltAfter  =  std(sigFiltCh(:,:,2));
    meanUnfiltBefore = mean(sigUnfiltCh(:,:,1));
    meanUnfiltAfter  = mean(sigUnfiltCh(:,:,2));
    errbUnfiltBefore =  std(sigUnfiltCh(:,:,1));
    errbUnfiltAfter  =  std(sigUnfiltCh(:,:,2));


    % Spectrogram Calculation --------------------------------------------
    [~, sFiltBeforeF,   sFiltBeforeT,   sFiltBefore] = ...
        spectrogram(meanFiltBefore,   sgramWinLen, [], sgramFQ, Fs);

    [~, sFiltAfterF,    sFiltAfterT,    sFiltAfter] = ...
        spectrogram(meanFiltAfter,    sgramWinLen, [], sgramFQ, Fs);

    [~, sUnfiltBeforeF, sUnfiltBeforeT, sUnfiltBefore] = ...
        spectrogram(meanUnfiltBefore, sgramWinLen, [], sgramFQ, Fs);

    [~, sUnfiltAfterF,  sUnfiltAfterT,  sUnfiltAfter] = ...
        spectrogram(meanUnfiltAfter,  sgramWinLen, [], sgramFQ, Fs);

    sFiltBefore   = zeros([size(sFiltBefore),   size(sigFiltCh,1)]);
    sFiltAfter    = zeros([size(sFiltAfter),    size(sigFiltCh,1)]);
    sUnfiltBefore = zeros([size(sUnfiltBefore), size(sigUnfiltCh,1)]);
    sUnfiltAfter  = zeros([size(sUnfiltAfter),  size(sigUnfiltCh,1)]);

    for trl = 1:size(sigFiltCh,1)
        [~,~,~,sFiltBefore(:,:,trl)] = spectrogram(sigFiltCh(trl,:,1), sgramWinLen, [], sgramFQ, Fs);
        [~,~,~,sFiltAfter(:,:,trl)] = spectrogram(sigFiltCh(trl,:,2), sgramWinLen, [], sgramFQ, Fs);
        if nUpdates
            if ~mod(trl, floor(size(sigFiltCh,1)/(nUpdates)))
                disp(['Obtaining Channel ',chan.labels,' Spectrograms: ',...
                      num2str(100*trl/(size(sigFiltCh,1)+size(sigUnfiltCh,1))),'%']);
            end
        end
    end

    for trl = 1:size(sigUnfiltCh,1)
        [~,~,~,sUnfiltBefore(:,:,trl)] = spectrogram(sigUnfiltCh(trl,:,1), sgramWinLen, [], sgramFQ, Fs);
        [~,~,~,sUnfiltAfter(:,:,trl)] = spectrogram(sigUnfiltCh(trl,:,2), sgramWinLen, [], sgramFQ, Fs);
        if nUpdates
            if ~mod(trl, floor(size(sigUnfiltCh,1)/(nUpdates)))
                disp(['Obtaining Channel ',chan.labels,' Spectrograms: ',...
                      num2str(100*(size(sigFiltCh,1)+trl)/(size(sigFiltCh,1)+size(sigUnfiltCh,1))),'%']);
            end
        end
    end

    sgramFiltBefore(chIdx,:)   = {chan, sFiltBeforeF,   sFiltBeforeT,   sFiltBefore};
    sgramFiltAfter(chIdx,:)    = {chan, sFiltAfterF,    sFiltAfterT,    sFiltAfter};
    sgramUnfiltBefore(chIdx,:) = {chan, sUnfiltBeforeF, sUnfiltBeforeT, sUnfiltBefore};
    sgramUnfiltAfter(chIdx,:)  = {chan, sUnfiltAfterF,  sUnfiltAfterT,  sUnfiltAfter};
    % ====================================================================


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
    meanSpectFiltBefore = mean(spectFiltBeforeCh); 
    errbSpectFiltBefore =  std(spectFiltBeforeCh);
    meanSpectFiltAfter = mean(spectFiltAfterCh); 
    errbSpectFiltAfter =  std(spectFiltAfterCh);
    meanSpectUnfiltBefore = mean(spectUnfiltBeforeCh); 
    errbSpectUnfiltBefore =  std(spectUnfiltBeforeCh);
    meanSpectUnfiltAfter = mean(spectUnfiltAfterCh); 
    errbSpectUnfiltAfter =  std(spectUnfiltAfterCh);

    fig(chIdx,1) = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig(chIdx,1)); sgtitle(['Channel ',chan.labels,' Avg. Response to Stim']);

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

end