function sigObj = doLPFilt(filtObj, sigObj)

useFF = filtObj.useFiltFilt;
if useFF
    filtFcn = @(f, x) filtfilt(f, x);
else
    filtFcn = @(f, x) filter(f, x);
end
lpFilt = filtObj.LPF;

%%

d = sigObj.Data_HPF; 
d_train = sigObj.Data_HPF_Train; 
d_test = sigObj.Data_HPF_Test; 

e = sigObj.Data_LMS; 
e_train = sigObj.Data_LMS_Train; 
e_test = sigObj.Data_LMS_Test;

%%

if numel(d)
    d_BPF = filtFcn(lpFilt, d);
    sigObj.Data_BPF = d_BPF;
else
    warning('Signal has not been highpass filtered.')
end

if numel(d_train)
    d_BPF_train = filtFcn(lpFilt, d_train);
    sigObj.Data_BPF_Train = d_BPF_train;
else
    disp('No highpass-filtered training data.')
end

if numel(d_test)
    d_BPF_test = filtFcn(lpFilt, d_test);
    sigObj.Data_BPF_Test = d_BPF_test;
else
    disp('No highpass-filtered testing data.')
end


if numel(e)
    e_LPF = filtFcn(lpFilt, e); 
    sigObj.Data_LMS_LPF = e_LPF;
else
    warning('Signal has not been adaptive filtered.')
end

if numel(e_train)
    e_LPF_train = filtFcn(lpFilt, e_train);
    sigObj.Data_LMS_LPF_Train = e_LPF_train;
else
    disp('No adaptive-filtered training data.')
end

if numel(e_test)
    e_LPF_test = filtFcn(lpFilt, e_test);
    sigObj.Data_LMS_LPF_Test = e_LPF_test;
else
    disp('No adaptive-filtered testing data.')
end

end