function [w, e_train, op_train, fig] = ...
    preTrainWts(t_train, g_train, d_train, N, uchan, nUpdates, dLMS)
% Use a subset of the signal to train optimal LMS weights. 
%
% Inputs: 
%   t_train: training time vectors as columns (s)
%   g_train: training noise reference signal as columns 
%   d_train: training unfiltered signal as columns 
%   N: number of filter taps 
%   uchan: array of unique channels, same length as width of g, d, and t
%   nUpdates: how many times to display progress and whether or not to plot weights. 
%             0 = no output and no plots 
%             Default = 10
%   dLMS: if false, does LMS; if true, does dLMS. Default false
% 
% Outputs: 
%   w: trained weights, as columns, from oldest --> current time 
%   e_train: LMS error signal for training period 
%   op_train: LMS output signal for training period 
%   fig: matlab figure with the trained weights of each channel

if nargin < 7
    dLMS = false;
end
%{
if nargin < 7
    nUpdates = 10;
    if nargin < 6
        uchan = 1:size(d,2);
    end
end
%}

%% organize training epochs 
G = zeros(size(t_train,1)-N+1, N, length(uchan)); 
T = zeros(size(G)); 
D = zeros(size(t_train,1)-N+1, length(uchan));
for idx = 1:length(uchan)
    ch = uchan(idx);
    D(:,idx) = d_train(N:size(t_train,1), idx);
    for nf = 1:(size(t_train,1)-N+1)
        if nUpdates
            if ~mod(nf, floor(size(t_train,1)/(nUpdates)))
                disp(['Building Channel ',ch.labels,' Training Matrix: ',num2str(100*nf/size(t_train,1)),'%']);
            end
        end
        G(nf,:,idx) = g_train(nf:(nf+N-1), idx);
        T(nf,:,idx) = t_train(nf:(nf+N-1), idx);
    end
end

%% training  
if nUpdates
    fig = figure('Units','normalized', 'Position',[.1 .1 .4 .8]);
end
w = zeros(N, length(uchan));
for idx = 1:length(uchan)
    ch = uchan(idx);
    Gidx = G(:,:,idx); Didx = D(:,idx);
    dG = diff(Gidx); dD = diff(Didx);
    if ~dLMS
        w(:,idx) = (((Gidx'*Gidx)^-1)*Gidx')*Didx;
    else
        w(:,idx) = (((dG'*dG)^-1)*dG')*dD;
    end
    if nUpdates
        figure(fig); subplot(length(uchan), 1, idx); stem(w(:,idx)); grid on;
        title(['Channel ',ch.labels,' training']);
        xlabel('tap'); ylabel('weight');
        pause(eps);
    end
end
if nUpdates
    pause(.5);
end

wnan = isnan(w);
if sum(wnan(:))
    w(wnan) = 0;
    wnan = sum(wnan(:))/numel(wnan);
    warning(['Weights are ',num2str(wnan*100),'% undefined; reverting to zeros.'])
end

%% post-processing  
op_train = zeros([size(t_train,1)-N+1,size(t_train,2)]); 
for idx = 1:length(uchan)
    op_train(:,idx) = G(:,:,idx)     *w(:,idx);
end

e_train = d_train; e_train(N:end,:) = e_train(N:end,:) - op_train;

end