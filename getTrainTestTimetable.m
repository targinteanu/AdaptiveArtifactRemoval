function [TTtrain,TTtest,TT] = getTrainTestTimetable(sig, filtObj)

N = filtObj.Num_Taps;

tTr = sig.Times_Train(:,1);
tTe = sig.Times_Test(:,1);
t = sig.Times(:,1);
u = sig.Noise_Reference(:,1); % change if there are unique inputs possible?

if ~numel(tTr)
    error('No training time specified.')
end
if ~numel(tTe)
    warning('No testing time specified.')
    TTtest = [];
end

if numel(sig.Data_LMS_LPF)
    y = sig.Data_LMS_LPF;
    yTr = sig.Data_LMS_LPF_Train;
    yTe = sig.Data_LMS_LPF_Test;
elseif numel(sig.Data_LMS)
    warning('Signal has not been post-filtered.')
    y = sig.Data_LMS;
    yTr = sig.Data_LMS_Train;
    yTe = sig.Data_LMS_Test;
elseif numel(sig.Data_BPF)
    warning('Signal has not been adaptive filtered.')
    N = 1; % do not remove first N datapoints if no LMS used 
    y = sig.Data_BPF;
    yTr = sig.Data_BPF_Train;
    yTe = sig.Data_BPF_Test;
elseif numel(sig.Data_HPF)
    warning('Signal has not been adaptive filtered or post-filtered.')
    N = 1;
    y = sig.Data_HPF;
    yTr = sig.Data_HPF_Train;
    yTe = sig.Data_HPF_Test;
elseif numel(sig.Data_Unfiltered)
    warning('Signal has not been filtered.')
    N = 1;
    y = sig.Data_Unfiltered;
    yTr = sig.Data_Unfiltered_Train;
    yTe = sig.Data_Unfiltered_Test;
else
    error('No data found.')
end

t = t(N:end,:);
if tTr(1) <= t(1)
    tTr = tTr(N:end,:);
end
u = u(N:end,:);
if numel(tTe)
    if tTe(1) <= t(1)
        tTe = tTe(N:end,:);
    end
    uTe = interp1(t, u, tTe);
end
uTr = interp1(t, u, tTr);

makeTT = @(t,u,y) [timetable(seconds(t),u), array2table(y)];

TT = makeTT(t, u, y);
VN = TT.Properties.VariableNames;
TTtrain = makeTT(tTr, uTr, yTr);
TTtrain.Properties.VariableNames = VN;
if numel(tTe)
    TTtest = makeTT(tTe, uTe, yTe);
    TTtest.Properties.VariableNames = VN;
end

end