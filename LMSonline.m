function [e_t, w, fig] = LMSonline(t, g, d, stepsize, N, uchan, w, nUpdates)
% Perform online LMS adaptive filtering. 
% 
% Inputs: 
%   t: time of g and d in columns 
%   g: noise reference signal in columns 
%   d: unfiltered signal in columns
%   stepsize: LMS step size / learning rate
%   N: number of filter taps 
%   uchan: array of unique channels, same length as width of g, d, and t
%   w: starting weights 
%      Default = zeros
%   nUpdates: how many times to display progress and whether or not to plot weights 
%             and smoothed error. 0 = no output and no plots 
%             Default = 100
% 
% Outputs: 
%   e_t: LMS error signal 
%   w: final weights 
%   fig: matlab figure showing the resulting weights and error signal 

if nargin < 7
    nUpdates = 100;
    if nargin < 6
        w = zeros(N,size(d,2));
        if nargin < 5
            uchan = 1:size(d,2);
        end
    end
end

if nUpdates
    fig = figure('Units','normalized', 'Position',[.1 .1 .8 .8]);
end

e_t = nan(size(t,1)-N+1, length(uchan));

for idx = 1:length(uchan)
    ch = uchan(idx);
    % train w: iterate grad descent
    if nUpdates
        figure(fig);
        subplot(length(uchan),2,2*idx-1); wplot = stem(w(:,idx));    grid on;
        title(['Channel ',ch.labels,' online']);
        xlabel('tap'); ylabel('weight');
        subplot(length(uchan),2,2*idx);   eplot = semilogy(e_t(:,idx)); grid on;
        title(['Channel ',ch.labels,' online']);
        xlabel('timepoint'); ylabel('e^2');
        pause(.5);
    end
    for ep = (N:size(t,1))-N+1
        Gidx = g((1:N)+ep-1, idx)';
        E = d(ep+N-1,idx) - Gidx*w(:,idx);
        e_t(ep, idx) = E;
        dw = E*Gidx';
        w(:,idx) = w(:,idx) + stepsize*dw;
        if nUpdates
            if ~mod(ep, floor(size(t,1)/nUpdates))
                wplot.YData = w(:,idx); eplot.YData = movmean(e_t(:,idx).^2, 5000);
                disp(['Online Channel ',ch.labels,': ',num2str(100*ep/size(t,1)),'%'])
                pause(eps);
            end
        end
    end
end

end