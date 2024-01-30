function sigObj = doPreFilt(filtObj, sigObj)

filtDirs = filtObj.preFiltDir;
preFilts = filtObj.preFilt;

d_unfilt = sigObj.Data_Unfiltered;
d_unfilt_train = sigObj.Data_Unfiltered_Train;
d_unfilt_test = sigObj.Data_Unfiltered_Test;

if numel(d_unfilt)
    d = d_unfilt;
else
    warning('No data to filter.')
end

if numel(d_unfilt_train)
    d_train = d_unfilt_train; 
else
    disp('No training data.')
end

if numel(d_unfilt_test)
    d_test = d_unfilt_test; 
else
    disp('No testing data.')
end

for idx = 1:length(preFilts)
    preFilt = preFilts(idx); 
    filtDir = filtDirs(idx);

    if ~filtDir
        filtFcn = @(f, x) filtfilt(f, x);
    elseif filtDir < 0
        filtFcn = @(f, x) flipud(filter(f, flipud(x))); 
    else
        filtFcn = @(f, x) filter(f, x);
    end
    
    if numel(d_unfilt)
        d = filtFcn(preFilt, d);        
    end
    
    if numel(d_unfilt_train)
        d_train = filtFcn(preFilt, d_train);        
    end
    
    if numel(d_unfilt_test)
        d_test = filtFcn(preFilt, d_test);
    else
        
    end

end

if numel(d_unfilt)
    sigObj.Data_HPF = d;
end
if numel(d_unfilt_train)
    sigObj.Data_HPF_Train = d_train;
end
if numel(d_unfilt_test)
    sigObj.Data_HPF_Test = d_test;
end

end