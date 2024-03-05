function [e_t, w, fig] = LMSonline(t, g, d, stepsize, N, uchan, w, nUpdates, dLMS, nLMS, resetThresh)
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
%   nLMS: if true, does NLMS; if false, does not normalize. Default false
% 
% Outputs: 
%   e_t: LMS error signal 
%   w: final weights 
%   fig: matlab figure showing the resulting weights and error signal 

if nargin < 11
    resetThresh = 1e100;
    %resetThresh = [];
    if nargin < 10
        nLMS = false;
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
    end
end

if nUpdates
    fig = figure('Units','normalized', 'Position',[.1 .1 .8 .8]);
end

if isempty(resetThresh)
    allowReset = false;
else
    allowReset = true;
    if length(resetThresh) == 1
        resetThresh = resetThresh*ones(size(uchan));
    end
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
    lastStableWts = w(:,idx); stabilityCounter = 0;
    for ep = (N:size(t,1))-N+1
        Gidx = g((1:N)+ep-1, idx)';
        E = d(ep+N-1,idx) - Gidx*w(:,idx);
        e_t(ep, idx) = E;
        if allowReset
            if abs(E) - abs(d(ep+N-1,idx)) > resetThresh(idx)
                % blowup 
                w(:,idx) = 0*lastStableWts; % fix
                E = d(ep+N-1,idx) - Gidx*w(:,idx);
                Eprev = d(ep+N-2,idx) - Gprev*w(:,idx);
                e_t(ep, idx) = E;
                if nUpdates
                    disp(['Weights reset at time ',num2str(t(ep,1)),' when e(t)=',num2str(E)])
                end
            end
        end
        if ~dLMS
            dw = E*Gidx'; % LMS
            if nLMS
                dw = dw./(Gidx*Gidx' + eps);
            end
            if ep > 1
                if abs(E) <= abs(Eprev)
                    % error is going down
                    stabilityCounter = stabilityCounter + 1;
                    if stabilityCounter >= N
                        % wts are stable
                        lastStableWts = w(:,idx);
                    end
                else
                    stabilityCounter = 0;
                end
            end
        else
            dG = Gidx-Gprev;
            dw = (E-Eprev) * dG'; % dLMS
            if nLMS
                dw = dw./(dG*dG' + eps);
            end
            if ep > 2
                if abs(e_t(ep,idx) - e_t(ep-1,idx)) <= abs(e_t(ep-1,idx) - e_t(ep-2,idx))
                    % d(error) is going down
                    stabilityCounter = stabilityCounter + 1;
                    if stabilityCounter >= N
                        % wts are stable
                        lastStableWts = w(:,idx);
                    end
                else
                    stabilityCounter = 0;
                end
            end
        end
        w(:,idx) = w(:,idx) + stepsize*dw;
        Gprev = Gidx; Eprev = E;
        if nUpdates
            if ~mod(ep, floor(size(t,1)/nUpdates))
                wplot.YData = w(:,idx); 
                %wplot.YData = lastStableWts;
                edata = e_t(:,idx);
                if dLMS
                    edata = [0;diff(edata)];
                end
                %eplot.YData = movmean(edata.^2, 5000);
                eplot.YData = edata.^2;
                disp(['Online Channel ',ch.labels,': ',num2str(100*ep/size(t,1)),'%'])
                pause(eps);
            end
        end
    end
end

end