function [sigObj, w_end] = LMSonlineWrapper(filtObj, sigObj, nUpdates)

if numel(sigObj.Data_LMS)
    disp('Overwriting Previous LMS-filtered Data');
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

N = filtObj.Num_Taps; 
stepsize = filtObj.Step_Size;

w_start = filtObj.Current_Weights;
if ~numel(w_start)
    w_start = filtObj.Starting_Weights;
    if ~numel(w_start)
        w_start = zeros(N, length(uchan));
        warning('No starting weights set. Using zeros.')
    end
end

[e_t, w_end] = LMSonline(t, g, sigUnfilt, stepsize, N, uchan, w_start, nUpdates, true);

sigObj.Data_LMS = e_t;

end