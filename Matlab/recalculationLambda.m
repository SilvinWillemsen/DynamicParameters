clear all;
close all;
clc
figure('Position', [100, 100, 800, 50])
k = 1/44100;
h = 1470 * k;
 N = floor(1/h);
for c = 1470:1:1500
    cla
    h = c * k;
%     N = floor(1/h);
    hRecalc = 1/N;
    lambda = c * k / hRecalc;
    text(0, 0, "$c = " + c + ",\quad h := " + num2str(h, 5) + ",$", 'verticalAlignment', 'middle', 'Fontsize', 18, 'interpreter', 'latex')
    text(0.36, 0, "$N := " + N + ",\quad h := "+ hRecalc + ",\quad \lambda := " + lambda + "$", 'verticalAlignment', 'middle', 'Fontsize', 18, 'interpreter', 'latex')
    axis off;
    
    ylim([-1, 1])
    set(gcf, 'color', 'w')
    pause(0.1)
end
