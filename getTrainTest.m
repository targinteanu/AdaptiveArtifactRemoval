function [t_train, g_train, d_train, t_test, g_test, d_test, T, G, D] = ...
    getTrainTest(t, g, d, trainTimeBounds, N, uchan, nUpdates)
% organize into testing and training 
% 
% Inputs: 
%   t: time of g and d in columns (s) 
%   g: noise reference signal in columns 
%   d: unfiltered signal in columns 
%   trainTimeBounds: training window as [Start_Time, End_time] (s) 
%   N: number of filter taps 
%   uchan: array of unique channels, same length as width of g, d, and t
%   nUpdates: how many times to display progress and whether or not to plot weights. 
%             0 = no output and no plots 
%             Default = 10
% 
% Outputs: 
%   t_train: training time vectors as columns (s)
%   g_train: training noise reference signal as columns 
%   d_train: training unfiltered signal as columns 
%   t_test: testing time vectors as columns (s). 
%           testing defined as all time that is not training time. 
%   g_test: testing noise reference signal as columns 
%   d_test: testing unfiltered signal as columns 
%   T: matrix of training time epochs  
%   G: matrix of training noise reference epochs 
%   D: matrix of training unfiltered epochs 

% split testing and training 
trainIdx = (t >= trainTimeBounds(1)) & (t <= trainTimeBounds(2)) ;
if sum(diff(sum(trainIdx)))
    error('Training windows are different sizes for different channels.')
end
testIdx = ~trainIdx;

t_train = zeros(sum(trainIdx(:,1)),size(trainIdx,2));
t_test  = zeros(sum(testIdx(:,1)), size(testIdx,2));
g_train = zeros(size(t_train)); g_test = zeros(size(t_test));
d_train = zeros(size(t_train)); d_test = zeros(size(t_test));

for ch = 1:length(uchan)
    trainIdxCh = trainIdx(:,ch); testIdxCh = testIdx(:,ch);
    t_train(:,ch) = t(trainIdxCh, ch); t_test(:,ch) = t(testIdxCh, ch);
    g_train(:,ch) = g(trainIdxCh, ch); g_test(:,ch) = g(testIdxCh, ch);
    d_train(:,ch) = d(trainIdxCh, ch); d_test(:,ch) = d(testIdxCh, ch);
end

% organize training epochs 
G = zeros(size(t_train,1)-N+1, N, length(uchan)); 
T = zeros(size(G)); 
D = zeros(size(t_train,1)-N+1, length(uchan));
for idx = 1:length(uchan)
    D(:,idx) = d_train(N:size(t_train,1), idx);
    for nf = 1:(size(t_train,1)-N+1)
        if nUpdates
            if ~mod(nf, floor(size(t_train,1)/(nUpdates)))
                disp(['Building Channel ',num2str(uchan(idx)),' Training Matrix: ',num2str(100*nf/size(t_train,1)),'%']);
            end
        end
        G(nf,:,idx) = g_train(nf:(nf+N-1), idx);
        T(nf,:,idx) = t_train(nf:(nf+N-1), idx);
    end
end

end