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
%{
figure; subplot(211); findpeaks(-x, t, 'MinPeakProminence', .1*xRng, 'Annotate', 'extents');
subplot(212); stem(ppr/xRng);
%}

pPksList = candidatePeaks(pPk, ppk, plc, pw, ppr);
nPksList = candidatePeaks(nPk, npk, nlc, nw, npr);

pOut = selectPeaks(pPk, pPksList); nOut = selectPeaks(nPk, nPksList);
pOut(1,:) = -pOut(1,:);

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
            tscore = abs(desloc - candi(:,2))./tRng;
            pscore = candi(:,4)./xRng; 
            score = pscore .* exp(-7*tscore);
            [~,sel] = max(score);
            lcOut(1,i) = candi(sel,1);
            lcOut(2,i) = candi(sel,2);
        end
    end

end