function [sigObj, filtObj] = preTrainWtsWrapper(filtObj, sigObj, nUpdates)

if numel(sigObj.Data_LMS_Train)
    disp('Overwriting Previous Pretrained Data');
end

sigUnfilt = sigObj.Data_HPF; 
if ~numel(sigUnfilt)
    sigUnfilt = sigObj.Data_Unfiltered; 
    if ~numel(sigUnfilt)
        error('No data to filter.');
    else
        warning('Signal has not been highpass filtered. Using Unfiltered instead (not recommended).');
    end
end

t = sigObj.Times; 
g = sigObj.Noise_Reference; 
uchan = sigObj.Channels;
trainTimeBounds = sigObj.Train_Time_Bounds;

N = filtObj.Num_Taps; 

[w, e_train, op_train, t_train, g_train, d_train, t_test, g_test, d_test, T, G, D] = ...
    preTrainWts(t, g, sigUnfilt, trainTimeBounds, N, uchan, nUpdates);

filtObj.Current_Weights = w;
sigObj.Data_LMS_Train = e_train; 
sigObj.Times_Train = t_train; 
sigObj.Noise_Reference_Train = g_train; 
sigObj.Data_Unfiltered_Train = d_train; 
sigObj.Times_Test = t_test; 
sigObj.Noise_Reference_Test = g_test; 
sigObj.Data_Unfiltered_Test = d_test; 

end