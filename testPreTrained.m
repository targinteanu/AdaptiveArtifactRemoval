function [e_test, op_test] = testPreTrained(w, t_test, g_test, d_test, N, uchan, nUpdates)
% Test a given set of weights on a specified part of the signal without
% adapting the weights.
% 
% Inputs:
%   w: weights, as columns, from oldest --> current time 
%   t_test: time vectors to use for testing as columns 
%   g_test: noise reference signal as columns 
%   d_test: unfiltered signal as columns 
%   N: number of taps 
%   uchan: array of unique channels, same length as width of g, d, and t
%   nUpdates: how many times to display progress. 0 = no output 
% 
% Outputs: 
%   e_test: LMS error signal
%   op_test: LMS output signal 

if nargin < 7
    nUpdates = 10;
    if nargin < 6
        uchan = 1:size(d,2);
    end
end

%% testing  
op_test = zeros(size(t_test,1)-N+1, length(uchan));
for idx = 1:length(uchan)
    ch = uchan(idx);
    for ep = (N:size(t_test,1))-N+1 
        if nUpdates
            if ~mod(ep, floor(size(t_test,1)/(nUpdates)))
                disp(['Testing Channel ',ch.labels,': ',num2str(100*ep/size(t_test,1)),'%']);
            end
        end
        Gidx = g_test((1:N)+ep-1, idx)';
        op_test(ep,idx) = Gidx*w(:,idx);
    end
end

e_test = d_test; e_test(N:end,:) = e_test(N:end,:) - op_test;

end