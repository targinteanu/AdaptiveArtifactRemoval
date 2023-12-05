function range = plotWithDistrib(x, y, dist, colr)
    % plot y with a dashed +- distribution surrounding y. 
    % y and dist must be row vectors
    plot(x, y, 'Color', colr, 'LineWidth',1); 
    hold on; 
    Y = y + [1;-1].*dist;
    plot(x, Y, ':', 'Color', colr, 'LineWidth',1);
    range = [min(Y(:)), max(Y(:))];
    range = [range; 1.25*[-1,1]*.5*diff(range) + mean(range)];
end