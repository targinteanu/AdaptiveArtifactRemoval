function sigObj = getTrainTestWrapper(sigObj)

t = sigObj.Times; 
g = sigObj.Noise_Reference; 
trainTimeBounds = sigObj.Train_Time_Bounds;
uchan = sigObj.Channels;

if numel(sigObj.Times_Train)
    disp('Overwriting previous training and testing data')
end

d = sigObj.Data_Unfiltered;
if ~numel(d)
    error('No data to filter.')
else

    [t_train, g_train, d_train, t_test, g_test, d_test] = ...
        getTrainTest(t, g, d, trainTimeBounds, uchan); 
    
    sigObj.Times_Train = t_train; 
    sigObj.Noise_Reference_Train = g_train; 
    sigObj.Data_Unfiltered_Train = d_train; 
    sigObj.Times_Test = t_test; 
    sigObj.Noise_Reference_Test = g_test; 
    sigObj.Data_Unfiltered_Test = d_test; 

end

d = sigObj.Data_HPF;
if ~numel(d)
    warning('Signal has not been highpass filtered.');
else

    [t_train, g_train, d_train, t_test, g_test, d_test] = ...
        getTrainTest(t, g, d, trainTimeBounds, uchan); 
    
    sigObj.Data_HPF_Train = d_train; 
    sigObj.Data_HPF_Test = d_test; 

end

d = sigObj.Data_BPF;
if ~numel(d)
    warning('Signal has not been bandpass filtered.');
else

    [t_train, g_train, d_train, t_test, g_test, d_test] = ...
        getTrainTest(t, g, d, trainTimeBounds, uchan); 
    
    sigObj.Data_BPF_Train = d_train; 
    sigObj.Data_BPF_Test = d_test; 

end

end