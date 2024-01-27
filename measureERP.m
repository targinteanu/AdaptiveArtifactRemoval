function [nOut, pOut, xstats] = measureERP(t, x, nPk, pPk, tRng, avgwinlen, showPlot)
% 
% Find peak latencies and properties from an ERP evoked potential signal.
% Input desired/approximate latencies; returns found peaks prioritizing
% both larger amplitude and closer priximity to desired latency. 
% 
% Inputs: 
%   t: time series 
%   x: signal corresponding to t
%   nPk: approximate latencies of upward-convex peaks 
%   pPk: approximate latencies of downward-convex peaks
%   tRng: range of possible latencies to search, in same units as t. 
%         Leave blank to use full range of input t 
%   avgwinlen: Window length for moving-average. If 0 or empty, no
%              averaging. Default = 0
%   showPlot: if true, displays plot of identified peaks. (default = false)
% 
% Outputs: 
%   nOut: found nPk as [amplitudes; latencies] 
%   pOut: found pPk as [amplitudes; latencies]
%   xstats: [mean; SD] of x within tRng 

if nargin < 7
    % default: if showPlot not specified, do not show 
    showPlot = false;
    if nargin < 6
        avgwinlen = 0;
    end
end

if ~isempty(tRng) 
    % if specified (nonempty), limit x and t to the allowed points 
    tIdx = (t >= tRng(1)) & (t <= tRng(2));
    t = t(tIdx);
    x = x(tIdx);
end
tRng = max(t) - min(t); % domain of all times 
xRng = max(x) - min(x); % range of signal 

xstats = [mean(x); std(x)];

if ~( isempty(avgwinlen) | (avgwinlen==0) )
    x = movmean(x, avgwinlen);
end

pOut = zeros(size(pPk)); 
nOut = zeros(size(nPk));

% get n peak candidates that exceed promience threshold compared to signal
% range (excludes very small peaks/oscillations) 
[npk, nlc, nw, npr] = findpeaks( x, t, 'MinPeakProminence', .01*xRng, 'MinPeakWidth', .001);
% flip the signal to get p peak candidates 
[ppk, plc, pw, ppr] = findpeaks(-x, t, 'MinPeakProminence', .01*xRng, 'MinPeakWidth', .001);

if ~isempty(pPk)
pPksList = candidatePeaks(pPk, ppk, plc, pw, ppr);
pOut = selectPeaks(pPk, pPksList); 
pOut(1,:) = -pOut(1,:);
end

if ~isempty(nPk)
nPksList = candidatePeaks(nPk, npk, nlc, nw, npr);
nOut = selectPeaks(nPk, nPksList);
end

if showPlot
% plot (for debugging)  
plot(t, x); grid on; 
hold on; 
if numel(nOut)
    plot(nOut(2,:),nOut(1,:),'^');
end
if numel(pOut)
    plot(pOut(2,:),pOut(1,:),'v');
end
xlabel('time'); ylabel('amplitude');
legend('signal', 'n peak(s)', 'p peak(s)');
end

    function pksList = candidatePeaks(locDesired, pk, lc, w, pr)
        % match found peaks to the nearest desired latency 
        % return a cell array corresponding to each desired latency 
        %   each cell contains a list of found peak candidates as rows;
        %   columns are [amplitude, latency, width, prominence] 
        pksList = cell(size(locDesired));
        for i = 1:length(pk)
            tdiff = abs(lc(i) - locDesired);
            [~,sel] = min(tdiff);
            pksList{sel} = [pksList{sel}; pk(i), lc(i), w(i), pr(i)];
        end
    end

    function lcOut = selectPeaks(des, cand)
        % for each desired latency, pick the "best" candidate peak 
        % return list of peaks as columns; rows are [amplitude; latency]
        lcOut = zeros(2, length(des));
        for i = 1:length(des)
            desloc = des(i); 
            candi = cand{i};
            if isempty(candi)
                % return nan amplitide and latency if no suitable candidate
                % near this desired latency was found. 
                lcOut(1,i) = nan;
                lcOut(2,i) = nan;
            else
                % "best" candidate defined by score preferring larger value
                % normalized to amplitude range (pscore) and exponentially 
                % penalizing latency distance normalized to time range (tscore)
                tscore = abs(desloc - candi(:,2))./tRng;
                pscore = candi(:,4)./xRng; 
                score = pscore .* exp(-4*tscore);
                [~,sel] = max(score);
                lcOut(1,i) = candi(sel,1);
                lcOut(2,i) = candi(sel,2);
            end
        end
    end

end