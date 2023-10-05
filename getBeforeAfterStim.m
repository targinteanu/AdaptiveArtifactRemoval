function [t_PrePost, d_PrePost, e_PrePost] = getBeforeAfterStim(tBeforeTrig, g, d, e, Fs, uchan, N)
% Plot the filtered and unfiltered signals and their power spectra some
% time before and after the stimulus. 
% 
% Inputs: 
%   tBeforeTrig: time before/after stimulus to plot (s) 
%   g: stimulus value over time, as columns 
%   d: unfiltered value over time, as columns (V)
%   e: filtered value over time, as columns (V)
%   Fs: sampling rate of g, d, and e_t (Hz) 
%   uchan: array of unique channels, same length as width of g, d, and e 
%   N: number of filter taps         
% 
% Outputs: 
%   t_PrePost: time with respect to stim as [before; after]
%   d_PrePost: unfiltered value as cell array for each channel; 
%              d_PrePost{i}(r,c,1) is before stim, trial r, timepoint c, channel i;
%              d_PrePost{i}(r,c,2) is after stim, trial r, timepoint c, channel i
%   e_PrePost: filtered value in the same format as d_PrePost 

%% getting signals before and after stim  
%tBeforeTrig = .29; % s
nBeforeTrig = floor(tBeforeTrig*Fs); % samples
tBeforeTrig = nBeforeTrig/Fs;
t_PrePost = [-tBeforeTrig:(1/Fs):0; 0:(1/Fs):tBeforeTrig]; % [before; after]

d_PrePost           = cell(1, length(uchan));
e_PrePost           = cell(size(d_PrePost));

for chIdx = 1:length(uchan)
    gch = g(:,chIdx);
    trig = [0; abs(diff(gch))];
    trig = trig > .1*max(trig); trig = find(trig);

    d_PrePost_ch           = nan(length(trig), nBeforeTrig+1, 2);
    e_PrePost_ch           = nan(size(d_PrePost_ch));

    for trIdx = 1:length(trig)
        tr = trig(trIdx); % timepoint

        % consider a more robust solution; padding e_t or making a separate
        % time vector
        % consider aligning so that rows of e correspond to d exactly 
        if (tr -nBeforeTrig > 0) & (tr +nBeforeTrig <= size(d,1))
            d_PrePost_ch(trIdx,:,1) = d(tr + ((-nBeforeTrig):0), chIdx);
            d_PrePost_ch(trIdx,:,2) = d(tr + (  0:nBeforeTrig ), chIdx);
        end
        if (tr -N+1 -nBeforeTrig > 0) & (tr -N+1 +nBeforeTrig <= size(e,1))
            e_PrePost_ch(trIdx,:,1) = e(tr -N+1 + ((-nBeforeTrig):0), chIdx);
            e_PrePost_ch(trIdx,:,2) = e(tr -N+1 + (  0:nBeforeTrig ), chIdx);
        end
    end

    toRemove = sum(e_PrePost_ch,2); 
    toRemove = sum(toRemove,3);
    toRemove = isnan(toRemove);
    e_PrePost_ch = e_PrePost_ch(~toRemove,:,:);
    toRemove = sum(d_PrePost_ch,2); 
    toRemove = sum(toRemove,3);
    toRemove = isnan(toRemove);
    d_PrePost_ch = d_PrePost_ch(~toRemove,:,:);

    d_PrePost{chIdx} = d_PrePost_ch;
    e_PrePost{chIdx} = e_PrePost_ch;
end

end