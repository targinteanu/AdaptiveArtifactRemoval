function [nOut, pOut] = measureERP(t, x, nPk, pPk, tRng)

if ~isempty(tRng) 
    tIdx = (t >= tRng(1)) & (t <= tRng(2));
    t = t(tIdx);
    x = x(tIdx);
end
tRng = max(t) - min(t);
xRng = max(x) - min(x);

pOut = zeros(size(pPk)); 
nOut = zeros(size(nPk));

[npk, nlc, nw, npr] = findpeaks( x, t, 'MinPeakProminence', .1*xRng);
[ppk, plc, pw, ppr] = findpeaks(-x, t, 'MinPeakProminence', .1*xRng);

if ~isempty(pPk)
pPksList = candidatePeaks(pPk, ppk, plc, pw, ppr);
pOut = selectPeaks(pPk, pPksList); 
pOut(1,:) = -pOut(1,:);
end

if ~isempty(nPk)
nPksList = candidatePeaks(nPk, npk, nlc, nw, npr);
nOut = selectPeaks(nPk, nPksList);
end

% plotting: add input to enable 
%{
figure; plot(t, x); grid on; 
hold on; 
plot(nOut(2,:),nOut(1,:),'^');
plot(pOut(2,:),pOut(1,:),'v');
%}

    function pksList = candidatePeaks(locDesired, pk, lc, w, pr)
        pksList = cell(size(locDesired));
        for i = 1:length(pk)
            tdiff = abs(lc(i) - locDesired);
            [~,sel] = min(tdiff);
            pksList{sel} = [pksList{sel}; pk(i), lc(i), w(i), pr(i)];
        end
    end

    function lcOut = selectPeaks(des, cand)
        lcOut = zeros(2, length(des));
        for i = 1:length(des)
            desloc = des(i); 
            candi = cand{i};
            if isempty(candi)
                lcOut(1,i) = nan;
                lcOut(2,i) = nan;
            else
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