function [av, phaseSpace] = getPhaseSpace(t, x, window, showPlotBool, filelabel)

if nargin < 5
    filelabel = '';
    if nargin < 4
        showPlotBool = false;
        if nargin < 3
            window = 10;
        end
    end
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