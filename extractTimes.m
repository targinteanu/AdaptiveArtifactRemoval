function sigOut = extractTimes(sigIn, timeBounds)
% gets signal object within specified time bounds. 
% 
% uses messy eval statements and functions. 
% currently does not work (debugging). 

nCh = length(sigIn.Channels);

% make timeBounds 3D
tb(1,1,2) = timeBounds(2); tb(1,1,1) = timeBounds(1);
timeBounds = tb;

%% full time window 

% get indexes 
timeSearch = sigIn.Times - timeBounds;
timeSearch(:,:,2) = -timeSearch(:,:,2);
idx = zeros(2,nCh);
for ch = 1:nCh
    for b = 1:2
        t = timeSearch(:,ch,b);
        t(t < 0) = inf;
        [~,idx(b,ch)] = min(t);
    end
end

% selection 
for ch = 1:nCh
    eval(asgnIdx('Times'));

    varnameslist = {'Data_Ground_Truth', 'Data_Unfiltered', 'Noise_Reference', ...
                    'Data_HPF', 'Data_LMS', 'Data_LMS_LPF', 'Data_BPF'};
    for v = 1:length(varnameslist)
        [cond, exif, exelse] = condAssign(varnameslist{v});
        if eval(cond)
            eval(exif);
        else
            eval(exelse);
        end
    end
end

%% training time window 

% get indexes 
if numel(sigIn.Times_Train)
timeSearch = sigIn.Times_Train - timeBounds;
timeSearch(:,:,2) = -timeSearch(:,:,2);
idx = zeros(2,nCh);
for ch = 1:nCh
    for b = 1:2
        t = timeSearch(:,ch,b);
        t(t < 0) = inf;
        [~,idx(b,ch)] = min(t);
    end
end
else
    idx = [];
end

% selection 
for ch = 1:nCh
    varnameslist = {'Times_Train', ...
                    'Data_Unfiltered_Train', 'Noise_Reference_Train', ...
                    'Data_HPF_Train', 'Data_LMS_Train', 'Data_LMS_LPF_Train', 'Data_BPF_Train'};
    for v = 1:length(varnameslist)
        [cond, exif, exelse] = condAssign(varnameslist{v});
        if eval(cond)
            eval(exif);
        else
            eval(exelse);
        end
    end
end

%% testing time window 

% get indexes 
if numel(sigIn.Times_Test)
timeSearch = sigIn.Times_Test - timeBounds;
timeSearch(:,:,2) = -timeSearch(:,:,2);
idx = zeros(2,nCh);
for ch = 1:nCh
    for b = 1:2
        t = timeSearch(:,ch,b);
        t(t < 0) = inf;
        [~,idx(b,ch)] = min(t);
    end
end
else
    idx = [];
end

% selection 
for ch = 1:nCh
    varnameslist = {'Times_Test', ...
                    'Data_Unfiltered_Test', 'Noise_Reference_Test', ...
                    'Data_HPF_Test', 'Data_LMS_Test', 'Data_LMS_LPF_Test', 'Data_BPF_Test'};
    for v = 1:length(varnameslist)
        [cond, exif, exelse] = condAssign(varnameslist{v});
        if eval(cond)
            eval(exif);
        else
            eval(exelse);
        end
    end
end

%% vars to copy
varnameslist = {'SampleRate', 'Channels', 'Preview_Channel_Index'};
for v = 1:length(varnameslist)
    eval(asgn(varnameslist{v}));
end

%% vars to adjust
varnameslist = {'Train_Time_Bounds', 'Preview_Time_Bounds'};
for v = 1:length(varnameslist)
    eval(handleTimeBound(varnameslist{v}));
end

%% helper functions - strings to eval

    function str = handleTimeBound(varname)
        str = [varname,' = sigIn.',varname,';', ...
               varname,'(1) = max(',varname,'(1), timeBounds(1);', ...
               varname,'(2) = min(',varname,'(2), timeBounds(2);', ...
               'sigOut.',varname,' = ',varname,';'];
    end

    function str = asgn(varname)
        str = ['sigOut.',varname,' = sigIn.',varname,';'];
    end

    function str = asgnIdx(varname)
        str = ['sigOut.',varname,'(:,ch) = sigIn.',varname,'(idx(1,ch):idx(2,ch),ch);'];
    end

    function [cond, exif, exelse] = condAssign(varname)
        cond = ['numel(sigIn.',varname,');'];
        exif = asgnIdx(varname);
        exelse = asgn(varname);
    end

end 