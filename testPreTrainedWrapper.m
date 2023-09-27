function sigObj = testPreTrainedWrapper(filtObj, sigObj, nUpdates)

if numel(sigObj.Data_LMS_Test)
    disp('Re-Testing Data');
end

sigUnfilt = sigObj.Data_HPF_Test; 
if ~numel(sigUnfilt)
    sigUnfilt = sigObj.Data_Unfiltered_Test; 
    if ~numel(sigUnfilt)
        error('No data to test.');
    else
        warning('Signal has not been highpass filtered. Using Unfiltered instead (not recommended).');
    end
end

w = filtObj.Current_Weights;
if ~numel(w)
    w = filtObj.Starting_Weights;
    if ~numel(w)
        error('No weights provided.')
    end
end

t_test = sigObj.Times_Test; 
g_test = sigObj.Noise_Reference_Test; 
uchan = sigObj.Channels;

N = filtObj.Num_Taps; 

e_test = testPreTrained(w, t_test, g_test, sigUnfiltTest, N, uchan, nUpdates);

sigObj.Data_LMS_Test = e_test;

end