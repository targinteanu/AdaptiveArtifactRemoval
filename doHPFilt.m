function sigObj = doHPFilt(filtObj, sigObj)

useFF = filtObj.useFiltFilt;
if useFF
    filtFcn = @(f, x) filtfilt(f, x);
else
    filtFcn = @(f, x) filter(f, x);
end
hpFilt = filtObj.HPF;

d_unfilt = sigObj.Data_Unfiltered;
d_unfilt_train = sigObj.Data_Unfiltered_Train;
d_unfilt_test = sigObj.Data_Unfiltered_Test;

if numel(d_unfilt)
    d = filtFcn(hpFilt, d_unfilt);
    sigObj.Data_HPF = d;
else
    warning('No data to filter.')
end

if numel(d_unfilt_train)
    d_train = filtFcn(hpFilt, d_unfilt_train);
    sigObj.Data_HPF_Train = d_train;
else
    disp('No training data.')
end

if numel(d_unfilt_test)
    d_test = filtFcn(hpFilt, d_unfilt_test);
    sigObj.Data_HPF_Test = d_test;
else
    disp('No testing data.')
end

end