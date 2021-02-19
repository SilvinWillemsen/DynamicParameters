close all;
alf = 0:0.001:1; 
maxEps = 30;
mat = zeros(length(alf), maxEps);
epsRange = 1:1:maxEps;
for eps = epsRange
    mat(:, eps) = (1-alf).^(eps);
end
for eps = epsRange
    plot(alf, mat(:, eps), 'b', 'Linewidth', 2);
    text(0.6, 0.6, "$\epsilon =\ $" + num2str(eps), ...
        "horizontalAlignment", "center", ...
        'interpreter', 'latex', 'Fontsize', 20)
    set(gca, 'Fontsize', 16, 'Linewidth', 2)
    grid on;
    xlabel("$\alpha$", 'interpreter', 'latex')
    ylabel("Displacement correction")
    pause(0.25);
    drawnow;
end