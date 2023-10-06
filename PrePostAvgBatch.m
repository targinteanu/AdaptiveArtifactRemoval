function fig = PrePostAvgBatch(nTrlInBatch, t_PrePost, d_PrePost, e_PrePost, Fs, uchan)
% Plot the filtered and unfiltered signals and their power spectra/spectrograms 
% before and after the stimulus, averaged in batches of trials of specified size.   
% 
% Inputs: 
%   nTrlInBatch: number of trials in each batch 
%   t_PrePost: time with respect to stim as [before; after]
%   d_PrePost: unfiltered value as cell array for each channel; 
%              d_PrePost{i}(r,c,1) is before stim, trial r, timepoint c, channel i;
%              d_PrePost{i}(r,c,2) is after stim, trial r, timepoint c, channel i
%   e_PrePost: filtered value in the same format as d_PrePost 
%   Fs: sampling rate of g, d, and e_t (Hz) 
%   uchan: array of unique channels, same length as width of g, d, and e
% 
% Outputs: 
%   fig: matlab figure, 4x2, showing time and spectrum
%        before/after stim and with/without filtering, averaged over 
%        batches of trials 


% TO DO: more outputs?  

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

%% Plotting averaged - in batches 
for chIdx = 1:length(uchan)
    sigFiltCh = e_PrePost{chIdx};
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

    fig(chIdx,1) = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig(chIdx,1)); sgtitle(['Channel ',num2str(uchan(chIdx)),' Batch Responses to Stim']);

    % Filtered: Plots ----------------------------------------------------
    nTrl = size(sigFiltCh,1);
    nBatch = floor(nTrl/nTrlInBatch);

    figure(fig(chIdx,1)); subplot(4,2,1); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(1,:), meanFiltBefore(batchIdx,:), errbFiltBefore(batchIdx,:), batchColor); 
    end
    title('Filtered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,1)); subplot(4,2,2); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(2,:), meanFiltAfter(batchIdx,:), errbFiltAfter(batchIdx,:), batchColor); 
    end
    title('Filtered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,1)); subplot(4,2,5); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(wFiltBefore, meanSpectFiltBefore(batchIdx,:), errbSpectFiltBefore(batchIdx,:), batchColor); 
    end
    title('Filtered Before'); set(gca, 'YScale', 'log'); grid on; 
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)'); 

    figure(fig(chIdx,1)); subplot(4,2,6); hold on;
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

    figure(fig(chIdx,1)); subplot(4,2,3); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(1,:), meanUnfiltBefore(batchIdx,:), errbUnfiltBefore(batchIdx,:), batchColor); 
    end
    title('Unfiltered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,1)); subplot(4,2,4); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(t_PrePost(2,:), meanUnfiltAfter(batchIdx,:), errbUnfiltAfter(batchIdx,:), batchColor); 
    end
    title('Unfiltered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 

    figure(fig(chIdx,1)); subplot(4,2,7); hold on;
    for batchIdx = 1:nBatch
        theta = batchIdx/nBatch;
        batchColor = theta*ltPink + (1-theta)*dkGreen;
        plotWithDistrib(wUnfiltBefore, meanSpectUnfiltBefore(batchIdx,:), errbSpectUnfiltBefore(batchIdx,:), batchColor); 
    end
    title('Unfiltered Before'); set(gca, 'YScale', 'log'); grid on; 
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)'); 

    figure(fig(chIdx,1)); subplot(4,2,8); hold on;
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

end