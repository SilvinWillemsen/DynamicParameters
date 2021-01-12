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
    subplot(1, 2, 1)
    hold off;
    plot(rangeToPlot, [0; sign(W(1, order(J))) * W(:,order(J)); zeroPad], '-o', 'Linewidth', 2)
    hold on;
    scatter(rangeToPlot(end-1), sign(W(1, order(J))) * W(end - (1 - length(zeroPad)), order(J)), 200, 'r', 'Linewidth', 1)
    grid on

    ylim([-1, 1]);
    set(gca, 'Linewidth', 2, 'Fontsize', 14)
    title("Modal Shape of Eigenvalue " + num2str(J), 'Fontsize', 14);

    subplot(1, 2, 2)
    hold off;

    plot(([i, i]-1) * 10000 / loopAmount, [1400, fs/2], '--', 'color', [0.5, 0.5, 0.5], 'Linewidth', 2)
    hold on;
    plotModesSave;
    hold on;
    scatter((i-1) * 10000 / loopAmount, 1/(2*pi*k) * acos(1/2 * D(order(J))), 40, 'k', 'Marker', 'o', 'Linewidth', 2)
    text(i * 10000 / loopAmount, 600, "$N = " + num2str(Ninit, 4) + "$",  "horizontalAlignment", 'center', 'interpreter', 'latex', 'Fontsize', 16);
    drawnow;
    
%     if modeToPlot == -1
%         pause(0.5);
%     end
end