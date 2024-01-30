function sigObj = doPostFilt(filtObj, sigObj)

filtDirs = filtObj.postFiltDir;
postFilts = filtObj.postFilt;

d = sigObj.Data_HPF; 
d_train = sigObj.Data_HPF_Train; 
d_test = sigObj.Data_HPF_Test; 

e = sigObj.Data_LMS; 
e_train = sigObj.Data_LMS_Train; 
e_test = sigObj.Data_LMS_Test;

%%

if numel(d)
    d_PrePost = d;
else
    warning('Signal has not been pre-filtered.')
end

if numel(d_train)
    d_PrePost_train = d_train;
else
    disp('No pre-filtered training data.')
end

if numel(d_test)
    d_PrePost_test = d_test;
else
    disp('No pre-filtered testing data.')
end


if numel(e)
    e_Post = e; 
else
    warning('Signal has not been adaptive filtered.')
end

if numel(e_train)
    e_Post_train = e_train;
else
    disp('No adaptive-filtered training data.')
end

if numel(e_test)
    e_Post_test = e_test;
else
    disp('No adaptive-filtered testing data.')
end

%%
for idx = 1:length(postFilts)
    lpFilt = postFilts(idx);
    filtDir = filtDirs(idx);

    if ~filtDir
        filtFcn = @(f, x) filtfilt(f, x);
    elseif filtDir < 0
        filtFcn = @(f, x) flipud(filter(f, flipud(x))); 
    else
        filtFcn = @(f, x) filter(f, x);
    end


    if numel(d)
        d_PrePost = filtFcn(lpFilt, d_PrePost);
    end

    if numel(d_train)
        d_PrePost_train = filtFcn(lpFilt, d_PrePost_train);
    end

    if numel(d_test)
        d_PrePost_test = filtFcn(lpFilt, d_PrePost_test);
    end


    if numel(e)
        e_Post = filtFcn(lpFilt, e_Post);
    end

    if numel(e_train)
        e_Post_train = filtFcn(lpFilt, e_Post_train);
    end

    if numel(e_test)
        e_Post_test = filtFcn(lpFilt, e_Post_test);
    end

end

%%
if numel(d)
    sigObj.Data_BPF = d_PrePost;
end
if numel(d_train)
    sigObj.Data_BPF_Train = d_PrePost_train;
end
if numel(d_test)
    sigObj.Data_BPF_Test = d_PrePost_test;
end

if numel(e)
    sigObj.Data_LMS_LPF = e_Post;
end
if numel(e_train)
    sigObj.Data_LMS_LPF_Train = e_Post_train;
end
if numel(e_test)
    sigObj.Data_LMS_LPF_Test = e_Post_test;
end

end