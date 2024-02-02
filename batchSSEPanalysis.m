%% SNR (signal : noise ratio)
tPost = t_PrePost(2,:);
ERPvals = cell(2, length(e_PrePost)); % rows = unfilt vs filt; cols = channel

% calculate values for all channels 
for ch = 1:length(size(ERPvals,2))

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

        figure; subplot(211); 

        % populate row trl of table with n20 & n40 amps and latencies and stats (mean & SD)
        [n20n40,~,stat] = measureERP(tPost, d_PrePost_ch(trl,:,2), [.02 .04], [], [.01,.1], 200, true);
        tblUnfilt.n20amp(trl) = n20n40(1,1); tblUnfilt.n20lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblUnfilt.n40amp(trl) = n20n40(1,2); tblUnfilt.n40lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblUnfilt.mean(trl) = stat(1); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = only within selected time range
            % noise = std thru all times
        clear n20n40 stat

        title('unfilt');
        subplot(212); 

        % do the same for filtered data 
        [n20n40,~,stat] = measureERP(tPost, e_PrePost_ch(trl,:,2), [.02 .04], [], [.01,.1], 200, true);
        tblFilt.n20amp(trl) = n20n40(1,1); tblFilt.n20lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblFilt.n40amp(trl) = n20n40(1,2); tblFilt.n40lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblFilt.mean(trl) = stat(1); tblFilt.SD(trl) = std(e_PrePost_ch(trl,:,2)); 
        clear n20n40 stat

        title('filt');
    end

    % store tables in cell array 
    ERPvals{1, ch} = tblUnfilt; ERPvals{2, ch} = tblFilt;
    clear tblFilt tblUnfilt e_PrePost_ch d_PrePost_ch;
end

%%
figure('Units','normalized', 'Position',[.1,.1,.8,.8])
for ch = 1:size(ERPvals,2)
    % get tables for this channel 
    tblUnfilt = ERPvals{1,ch}; tblFilt = ERPvals{2,ch}; 
    % "Signal" = n20 amplitude - mean for [unfilt, filt]
    n20amp = [tblFilt.n20amp] - [tblFilt.mean];
    n40amp = [tblFilt.n40amp] - [tblFilt.mean];
    amp = [n20amp, n40amp];
    % loc = n20, n40 latency for unfilt, filt 
    loc = [tblUnfilt.n20lat, tblFilt.n20lat, tblUnfilt.n40lat, tblFilt.n40lat];
    numfound = sum(~isnan(loc));
    % hypothesis test whether unfilt and filt have different locs 
    [~,ploc20] = ttest(loc(:,1), loc(:,2), 'Tail', 'both');
    [~,ploc40] = ttest(loc(:,3), loc(:,4), 'Tail', 'both');

    % box plot results with p values displayed  
    ax1(ch) = subplot(size(ERPvals,2),2,2*(ch-1)+1); boxplot(amp);
    grid on;
    title('amp n20, n40'); ylabel(['Amplitude ',sig.Channels(ch).labels]);
    xticklabels({'n20', 'n40'});
    ax2(ch) = subplot(size(ERPvals,2),2,2*ch); boxplot(loc);
    grid on;
    title(['Latency: p = ',num2str(ploc20),' (n20), ',num2str(ploc40),' (n40)']); 
    xtl = {'Unfilt n20', 'Filt n20', 'Unfilt n40', 'Filt n40'};
    xtl = arrayfun(@(i) [xtl{i},': N=',num2str(numfound(i))], 1:length(xtl), 'UniformOutput',false);
    xticklabels(xtl);
end
linkaxes(ax1,'y'); linkaxes(ax2,'y');
%}