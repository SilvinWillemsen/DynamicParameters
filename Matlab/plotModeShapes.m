% to be used in the modalAnalysis.m file

[~, order] = sort(D, 'descend');
if floor(Ninit) ~= Ninit
    rangeToPlot = ([0:floor(Ninit), Ninit]) / Ninit;
    endLoop = Ninit;
    zeroPad = 0;
else
    rangeToPlot = (0:Ninit) / Ninit;
    endLoop = length(D~=0);
    zeroPad = [];
end
if modeToPlot == -1
    loopRange = 1:endLoop;
    loopRange = loopRange(D~=0);
else
    loopRange = modeToPlot;
end
subplotIdx = 1;
hold off;
for J = loopRange
%     if length(loopRange) <= limSubplots
%         subplot (length(loopRange), 1, subplotIdx);
%         subplotIdx = subplotIdx + 1;
%     end
    subplot(1, 2, 1)
    plot(rangeToPlot, [0; sign(W(1, order(J))) * W(:,order(J)); zeroPad], 'Linewidth', 2)
    grid on
%     if length(loopRange) > limSubplots
% %         hold on;
%     end
    ylim([-0.6, 0.6])
    set(gca, 'Linewidth', 2, 'Fontsize', 14)
    title("Modal Shapes: N = " + num2str(Ninit) + "   Eigen value " + num2str(J) + " = " + ...
        num2str(round(1/(2 * pi * k) * acos (1/2 * D(order(J)))) ...
        + " Hz"), 'Fontsize', 14);

    subplot(1, 2, 2)
    hold off;

    plot(([i, i]-1) * 10000 / loopAmount, [0, fs/2], '--', 'color', [0.5, 0.5, 0.5], 'Linewidth', 2)
    hold on;
    plotModesSave;
    hold on;
    scatter((i-1) * 10000 / loopAmount, 1/(2*pi*k) * acos(1/2 * D(order(J))), 40, 'k', 'Marker', 'o', 'Linewidth', 2)
    
    drawnow;
    
%     if length(loopRange) > limSubplots || modeToPlot == -1
%         pause(0.5);
%     end
end