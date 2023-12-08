function sigOut = extractChannel(sigIn, chan)
% Get signal object of just one channel, chan.  
% chan can be an index or a channel object. 
% 
% uses eval statements, which is messy and not preferred 

chlist = sigIn.Channels;

if isa(chan, 'struct')
    sigOut.Channels = chan;
    chanIdx = find(strcmp(chan.labels, {chlist.labels}));
    chan = chanIdx;
else
    sigOut.Channels = chlist(chan);
end
% chan is now an index. 
sigOut.Preview_Channel_Index = chan;

strSelectChan = @(varname) ['sigOut.',varname,' = sigIn.',varname,'(:,chan);'];
eval(strSelectChan('Times'));

strCond = @(varname) ['numel(sigIn.',varname,');'];

strAsgn = @(varname) ['sigOut.',varname,' = sigIn.',varname,';'];
varnameslist = {'SampleRate', 'Train_Time_Bounds', 'Preview_Time_Bounds'};
for v = 1:length(varnameslist)
    eval(strAsgn(varnameslist{v}));
end

varnameslist = {'Data_Ground_Truth', ...
                'Data_Unfiltered', 'Noise_Reference', ...
                'Data_HPF', 'Data_LMS', 'Data_LMS_LPF', 'Data_BPF', ...
                'Times_Train', 'Times_Test', 'Noise_Reference_Train', 'Noise_Reference_Test', ...
                'Data_Unfiltered_Train', 'Data_Unfiltered_Test', ...
                'Data_HPF_Train', 'Data_HPF_Test', 'Data_LMS_Train', 'Data_LMS_Test', ...
                'Data_LMS_LPF_Train', 'Data_LMS_LPF_Test', 'Data_BPF_Train', 'Data_BPF_Test'};
for v = 1:length(varnameslist)
    varname = varnameslist{v};
    if eval(strCond(varname))
        eval(strSelectChan(varname));
    else
        eval(strAsgn(varname));
    end
end

end