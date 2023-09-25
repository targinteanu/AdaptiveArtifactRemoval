function [e_test, op_test] = testPreTrained(w, t_test, g_test, d_test, N, uchan, nUpdates)

%% testing  
op_test = zeros(size(t_test,1)-N+1, length(uchan));
for idx = 1:length(uchan)
    for ep = (N:size(t_test,1))-N+1 
        if nUpdates
            if ~mod(ep, floor(size(t_test,1)/(.1*nUpdates)))
                disp(['Testing Channel ',num2str(uchan(idx)),': ',num2str(100*ep/size(t_test,1)),'%']);
            end
        end
        Gidx = g_test((1:N)+ep-1, idx)';
        op_test(ep,idx) = Gidx*w(:,idx);
    end
end

e_test = d_test; e_test(N:end,:) = e_test(N:end,:) - op_test;

end