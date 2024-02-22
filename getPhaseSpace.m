function [av, phaseSpace] = getPhaseSpace(t, x, tRng, window, showPlotBool, filelabel)

if nargin < 6
    filelabel = '';
    if nargin < 5
        showPlotBool = false;
        if nargin < 4
            window = 10;
            if nargin < 3
                tRng = [];
            end
        end
    end
end

if ~isempty(tRng)
    tIdx = (t >= tRng(1)) & (t <= tRng(2));
    t = t(tIdx);
    x = x(tIdx);
end

    x=movmean(x,window);
    deriv = diff(x)./diff(t);
    phaseSpace = [x(1:(end-1));deriv]';
    [k,av] = convhull(phaseSpace, 'Simplify', true);

    if showPlotBool == true
       plot(phaseSpace(:,1),phaseSpace(:,2),".")
       title(filelabel); xlabel('x'); ylabel('dx/dt');
       hold on
       plot(phaseSpace(k,1),phaseSpace(k,2))
    end

end