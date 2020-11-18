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
    if length(loopRange) <= limSubplots
        subplot (length(loopRange), 1, subplotIdx);
        subplotIdx = subplotIdx + 1;
    end
    plot(rangeToPlot, [0; W(:,order(J)); zeroPad])
    if length(loopRange) > limSubplots
        hold on;
    end
    title("N = " + num2str(Ninit) + " Eigen value " + num2str(J) + " = " + ...
        num2str(round(1/(2 * pi * k) * acos (1/2 * D(order(J)))) ...
        + " Hz"));
    ylim([-1, 1])
    drawnow;
    
    if length(loopRange) > limSubplots || modeToPlot == -1
        pause(0.5);
    end
end