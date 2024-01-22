function [e_t, w, fig] = LMSonline(t, g, d, stepsize, N, uchan, w, nUpdates, dLMS)
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
%   dLMS: if false, does LMS; if true, does dLMS. Default false
% 
% Outputs: 
%   e_t: LMS error signal 
%   w: final weights 
%   fig: matlab figure showing the resulting weights and error signal 

if nargin < 9
    dLMS = false;
    if nargin < 8
        nUpdates = 100;
        if nargin < 7
            w = zeros(N,size(d,2));
            if nargin < 6
                uchan = 1:size(d,2); % need to fix
            end
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
        xlabel('timepoint'); 
        if ~dLMS
            ylabel('e^2');
        else
            ylabel('(de/dt)^2');
        end
        pause(.5);
    end
    Gprev = zeros(1,N); Eprev = zeros(1,N);
    for ep = (N:size(t,1))-N+1
        Gidx = g((1:N)+ep-1, idx)';
        E = d(ep+N-1,idx) - Gidx*w(:,idx);
        e_t(ep, idx) = E;
        if ~dLMS
            dw = E*Gidx'; % LMS
        else
            dw = (E-Eprev) * (Gidx-Gprev)'; % dLMS
        end
        w(:,idx) = w(:,idx) + stepsize*dw;
        Gprev = Gidx; Eprev = E;
        if nUpdates
            if ~mod(ep, floor(size(t,1)/nUpdates))
                wplot.YData = w(:,idx); 
                edata = e_t(:,idx);
                if dLMS
                    edata = [0,diff(edata)];
                end
                eplot.YData = movmean(edata.^2, 5000);
                disp(['Online Channel ',ch.labels,': ',num2str(100*ep/size(t,1)),'%'])
                pause(eps);
            end
        end
    end
end

end