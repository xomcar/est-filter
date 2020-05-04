function b = testSchur(poly1, poly2)
    %convolution
    prod = conv(poly1, poly2);
    
    %schur evalutation
    if abs(roots(prod)) < 1
        b = true;
    else
        b = false;
    end
    
    %plotting
    t = -10:0.01:10;
    hold on
    plot(t, polyval(prod, t), ...
        'LineStyle' ,'-', ...
        'LineWidth', 2,...
        'Color', [1 0 0] ...
    );    
    %axes management
    axes = gca;
    axes.FontUnits = 'points';
    axes.FontSize = 18;
    axes.Title.String = 'Product of p(x) and q(x)';
    axes.YLabel.Interpreter = 'latex';
    axes.YLabel.String = '$p(x) * q(x)$';
    axes.XLabel.Interpreter = 'latex';
    axes.XLabel.String = '$x$';
    axes.XLim = [-12 12];
    axes.XGrid = 'on';
    axes.YGrid = 'on';
    axes.GridLineStyle = '-.';
end