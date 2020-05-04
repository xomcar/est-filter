function plotSpectrum(tranFun)
    %plots spectrum of tranFun if the function D is Schur
    %check for stability of transfer function
    a = checkTFStability(tranFun);
    if (~a)
        %exit if function is not stable
        fprintf('Function is not stable! Exiting...\n')
        return
    end
    fprintf('Function is stable! Plotting spectrum...\n')
    %computing transfer function with z = 1/z
    [n, d] = tfdata(tranFun, 'v');
    %flipping num and denum coefficients
    n = flip(n);
    d = flip(d);
    %regenerating transfer function
    tranFunMinus = tf(n, d, -1);
    %generate a gaussian noise
    sigma = 1;
    %compute spectral density by transforming the variance
    Sy = tranFun * sigma^2 * tranFunMinus;
    %x axys for plotting
    t = [-pi:0.001:pi];
    %corresponding unit circle
    z = exp(1j*t);
    %evaluation of spectral density in unit circle
    ev = zeros(1, length(z));
    for i = 1:length(z)
        ev(i) = evalfr(Sy, z(i));
    end
    %estracting absolute value
    ev = abs(ev);
    %plotting
    figure
    plot(t,ev, ...
        'LineStyle' ,'-', ...
        'LineWidth', 2,...
        'Color', [1 0 0] ...
    );

    %axes management
    axes = gca;
    axes.FontUnits = 'points';
    axes.FontSize = 18;
    axes.Title.String = 'Spectral density'; 
    axes.YLabel.Interpreter = 'latex';
    axes.YLabel.String = '';
    axes.XLabel.Interpreter = 'latex';
    axes.XLabel.String = '$\theta \in [-\pi,+\pi]$';
    axes.XLim = [-3.14 3.14];
    axes.XGrid = 'on';
    axes.YGrid = 'on';
    axes.GridLineStyle = '-.';
end