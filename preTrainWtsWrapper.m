function [sigObj, filtObj] = preTrainWtsWrapper(filtObj, sigObj, nUpdates)

if numel(sigObj.Data_LMS_Train)
    disp('Overwriting Previous Pretrained Data');
end

sigUnfilt = sigObj.Data_HPF_Train; 
if ~numel(sigUnfilt)
    sigUnfilt = sigObj.Data_Unfiltered_Train; 
    if ~numel(sigUnfilt)
        error('No data to filter.');
    else
        warning('Signal has not been highpass filtered. Using Unfiltered instead (not recommended).');
    end
end

uchan = sigObj.Channels;
t_train = sigObj.Times_Train;
g_train = sigObj.Noise_Reference_Train;

N = filtObj.Num_Taps; 

[w, e_train, op_train] = ...
    preTrainWts(t_train, g_train, sigUnfilt, N, uchan, nUpdates);

filtObj.Current_Weights = w;
sigObj.Data_LMS_Train = e_train; 

end