function [w, e_train, op_train, t_train, g_train, d_train, t_test, g_test, d_test, T, G, D] = ...
    preTrainWts(t, g, d, trainTimeBounds, N, uchan, dispTrainProgress)

%% organize into testing and training 

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
        G(nf,:,idx) = g_train(nf:(nf+N-1), idx);
        T(nf,:,idx) = t_train(nf:(nf+N-1), idx);
    end
end

%% training  
if dispTrainProgress
    fig = figure('Units','normalized', 'Position',[.1 .1 .4 .8]);
end
w = zeros(N, length(uchan));
for idx = 1:length(uchan)
    Gidx = G(:,:,idx); Didx = D(:,idx);
    w(:,idx) = (((Gidx'*Gidx)^-1)*Gidx')*Didx;
    if dispTrainProgress
        figure(fig); subplot(length(uchan), 1, idx); stem(w(:,idx)); grid on;
        title(['Channel ',num2str(uchan(idx)),' training']);
        xlabel('tap'); ylabel('weight');
        pause(eps);
    end
end
if dispTrainProgress
    pause(.5);
end

%% post-processing  
op_train = zeros([size(t_train,1)-N+1,size(t_train,2)]); 
for idx = 1:length(uchan)
    op_train(:,idx) = G(:,:,idx)     *w(:,idx);
end

e_train = d_train; e_train(N:end,:) = e_train(N:end,:) - op_train;

end