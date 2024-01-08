function [pkOut, lcOut] = measureERP(t, x, tPk, sPk, tRng)

if ~(numel(tPk)==numel(sPk))
    error('tPk and sPk must be same size.')
end

if ~isempty(tRng) 
    tIdx = (t >= tRng(1)) & (t <= tRng(2));
    t = t(tIdx);
    x = x(tIdx);
end
tRng = max(t) - min(t);
xRng = max(x) - min(x);

lcOut = zeros(size(tPk)); 
pkOut = zeros(size(sPk));

[npk, nlc, nw, npr] = findpeaks( x, t, 'MinPeakProminence', .1*xRng);
%{
figure; subplot(211); findpeaks( x, t, 'MinPeakProminence', .1*xRng, 'Annotate', 'extents');
subplot(212); stem(npr/xRng);
%}
[ppk, plc, pw, ppr] = findpeaks(-x, t, 'MinPeakProminence', .1*xRng);

pkslist = cell(size(tPk));
for idx = 1:length(tPk)
    if sPk(idx) > 0
        pk = npk;  lc = nlc; w = nw; pr = npr;
    elseif sPk(idx) < 0
        pk = -ppk; lc = plc; w = pw; pr = ppr;
    end
    tscore = abs(lc - tPk(idx))./tRng;
    pscore = abs(pk)./xRng; 
    [~,selIdx] = min(tscore);
    clear pk lc w pr
end
%}

end