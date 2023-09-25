function plotBeforeAfterStim(tBeforeTrig, g, d, e_t, Fs, uchan, N)
% Plot the filtered and unfiltered signals and their power spectra some
% time before and after the stimulus. 
% 
% tBeforeTrig: time before/after stimulus to plot (s) 
% g: stimulus value over time, as columns 
% d: unfiltered value over time, as columns (V)
% e_t: filtered value over time, as columns (V)
% Fs: sampling rate of g, d, and e_t (Hz) 
% uchan: array of unique channels, same length as width of g, d, and e_t 
% N: number of filter taps                                                 <--- remove by making more robust(?)

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

for idx = 1:length(uchan)
    gch = g(:,idx);
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
            d_PrePost_ch(trIdx,:,1) = d(tr + ((-nBeforeTrig):0), idx);
            d_PrePost_ch(trIdx,:,2) = d(tr + (  0:nBeforeTrig ), idx);
%            d_lpf_PrePost_ch(trIdx,:,1) = d_lpf(tr + ((-nBeforeTrig):0), idx);
%            d_lpf_PrePost_ch(trIdx,:,2) = d_lpf(tr + (  0:nBeforeTrig ), idx);
        end
        if (tr -N+1 -nBeforeTrig > 0) & (tr -N+1 +nBeforeTrig <= size(e_t,1))
            e_t_PrePost_ch(trIdx,:,1) = e_t(tr -N+1 + ((-nBeforeTrig):0), idx);
            e_t_PrePost_ch(trIdx,:,2) = e_t(tr -N+1 + (  0:nBeforeTrig ), idx);
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

    d_PrePost{idx} = d_PrePost_ch;
%    d_lpf_PrePost{idx} = d_lpf_PrePost_ch;
%    e_train_PrePost{idx} = e_train_PrePost_ch;
%    e_train_lpf_PrePost{idx} = e_train_lpf_PrePost_ch;
%    e_test_PrePost{idx} = e_test_PrePost_ch;
%    e_test_lpf_PrePost{idx} = e_test_lpf_PrePost_ch;
    e_t_PrePost{idx} = e_t_PrePost_ch;
%    e_t_lpf_PrePost{idx} = e_t_lpf_PrePost_ch;
end

%% cleanup 
clear d_PrePost_ch d_lpf_PrePost_che_train_PrePost_ch e_train_lpf_PrePost_ch 
clear e_test_PrePost_ch e_test_lpf_PrePost_ch e_t_PrePost_ch e_t_lpf_PrePost_ch 
clear gch trig trIdx tr toRemove

%% plotting averaged signals before and after stim 
% colors: 
dkBlue  = [  1,  50, 130] /255;
ltBlue  = [145, 190, 255] /255;
dkRed   = [110,   0,   0] /255;
ltRed   = [250, 150, 150] /255;
dkBlack = [  0,   0,   0] /255;
ltBlack = [110, 110, 110] /255;

for idx = 1:length(uchan)
    sigFiltCh = e_t_PrePost{idx};
    sigUnfiltCh = d_PrePost{idx};

    meanFiltBefore = mean(sigFiltCh(:,:,1));
    meanFiltAfter  = mean(sigFiltCh(:,:,2));
    errbFiltBefore =  std(sigFiltCh(:,:,1));
    errbFiltAfter  =  std(sigFiltCh(:,:,2));
    meanUnfiltBefore = mean(sigUnfiltCh(:,:,1));
    meanUnfiltAfter  = mean(sigUnfiltCh(:,:,2));
    errbUnfiltBefore =  std(sigUnfiltCh(:,:,1));
    errbUnfiltAfter  =  std(sigUnfiltCh(:,:,2));

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

    fig = figure('Units','normalized', 'Position',[.1 .1 .8 .8]); 
    figure(fig); sgtitle(['Chennel ',num2str(uchan(idx)),' Avg. Response to Stim']);

    figure(fig); subplot(321); 
           plotWithDistrib(t_PrePost(1,:), meanUnfiltBefore, errbUnfiltBefore, ltRed);
    yrng = plotWithDistrib(t_PrePost(1,:), meanFiltBefore, errbFiltBefore, dkBlue);
    title('Filtered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)'); 
    ylim(yrng(2,:)); 
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig); subplot(322); 
           plotWithDistrib(t_PrePost(2,:), meanUnfiltAfter, errbUnfiltAfter, ltRed);
    yrng = plotWithDistrib(t_PrePost(2,:), meanFiltAfter, errbFiltAfter, dkBlue);
    title('Filtered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)');
    ylim(yrng(2,:));
    legend('Unfiltered', '-1SD', '+1SD', 'Filtered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig); subplot(323);  
           plotWithDistrib(t_PrePost(1,:), meanFiltBefore, errbFiltBefore, ltBlue);
    yrng = plotWithDistrib(t_PrePost(1,:), meanUnfiltBefore, errbUnfiltBefore, dkRed);
    title('Unfiltered Before'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)');
    ylim(yrng(2,:));
    legend('Filtered', '-1SD', '+1SD', 'Unfiltered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig); subplot(324);  
           plotWithDistrib(t_PrePost(2,:), meanFiltAfter, errbFiltAfter, ltBlue);
    yrng = plotWithDistrib(t_PrePost(2,:), meanUnfiltAfter, errbUnfiltAfter, dkRed);
    title('Unfiltered After'); grid on; 
    xlabel('time (s)'); ylabel('Signal (V)');
    ylim(yrng(2,:));
    legend('Filtered', '-1SD', '+1SD', 'Unfiltered', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig); subplot(325); 
    semilogy(wUnfiltAfter, meanSpectUnfiltAfter, 'Color', ltRed); hold on; 
    semilogy(wFiltAfter, meanSpectFiltAfter, 'Color', ltBlue);
    plotWithDistrib(wUnfiltBefore, meanSpectUnfiltBefore, errbSpectUnfiltBefore, dkRed);
    plotWithDistrib(wFiltBefore, meanSpectFiltBefore, errbSpectFiltBefore, dkBlue);
    title('Spectrum Before'); grid on; 
    %set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)');
    legend('Unfiltered After', 'Filtered After', ...
        'Unfiltered Before', '-1SD', '+1SD', 'Filtered Before', '-1SD', '+1SD', 'Location','eastoutside');

    figure(fig); subplot(326); 
    semilogy(wUnfiltBefore, meanSpectUnfiltBefore, 'Color', ltRed); hold on; 
    semilogy(wFiltBefore, meanSpectFiltBefore, 'Color', ltBlue);
    plotWithDistrib(wUnfiltAfter, meanSpectUnfiltAfter, errbSpectUnfiltAfter, dkRed);
    plotWithDistrib(wFiltAfter, meanSpectFiltAfter, errbSpectFiltAfter, dkBlue);
    title('Spectrum After'); grid on; 
    %set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)'); ylabel('Power Spectrum (V^2*s^2)');
    legend('Unfiltered Before', 'Filtered Before', ...
        'Unfiltered After', '-1SD', '+1SD', 'Filtered After', '-1SD', '+1SD', 'Location','eastoutside');
end

%% helper functions 
function range = plotWithDistrib(x, y, dist, colr)
    % plot y with a dashed +- distribution surrounding y. 
    % y and dist must be row vectors
    plot(x, y, 'Color', colr); 
    hold on; 
    Y = y + [1;-1].*dist;
    plot(x, Y, ':', 'Color', colr);
    range = [min(Y(:)), max(Y(:))];
    range = [range; 1.25*[-1,1]*.5*diff(range) + mean(range)];
end

function [wP, P, w, Y] = PowerSpectrum(y, Fs)
    % y: a row vector/matrix in time domain 
    % wP: one-sided frequency 
    % P: power spectrum (1-sided) of Y
    % w: two-sided frequency 
    % Y: frequency spectrum (2-sided, complex)
    L = size(y, 2);
    y = y';
    Y = fft(y);  
    w = Fs/L*(-L/2:L/2-1);
    P = abs(Y/L)'; P = P(:, 1:L/2+1);
    P(:, 2:end-1) = 2*P(:, 2:end-1);
    P = P.^2;
    wP = Fs/L*(0:(L/2));
    Y = fftshift(Y)';
end

end