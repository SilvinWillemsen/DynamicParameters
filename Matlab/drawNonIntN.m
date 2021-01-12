close all

figure1 = figure('Color',[1 1 1]);
% Create axes
axes1 = axes('Parent',figure1,'Position',[0 0 1 1]);

markerWidth = 0.05;
height = 0;
curColor = 'r';

jRange = 6:0.25:7;

plot([0, 0], [0.2, -length(jRange)+1.15], '--', 'color', [0.5, 0.5, 0.5])
hold on;
plot([1, 1], [0.2, -length(jRange)+1.15], '--', 'color', [0.5, 0.5, 0.5])

addingAt = "middle";
for j = jRange
    range = (0:j)/j;
    if range(end) ~= 1
        if  addingAt == "boundary"
            range = [range, 1];
        else
            range = [range(1:ceil(length(range)/2)), range(ceil(length(range)/2):end) + (1-range(end))];
        end
    end
    curPlot = plot(range, ones(ceil(j)+1, 1) * height, curColor, 'Linewidth', 2);

    for i = range
        plot([i,i], [height-markerWidth, height+markerWidth], 'Color', curColor, 'Linewidth', 2);
    end
    
    text(0.5, height + 0.2, "$N ="+ j +"$", 'interpreter', 'latex', 'Fontsize', 18 ,'horizontalAlignment', 'center')
    text(range(2) / 2, height - 0.2, "$h$", 'interpreter', 'latex', 'Fontsize', 18 ,'horizontalAlignment', 'center', 'color', curColor)
    
    if curColor == 'r'
        curColor = 'b'
    else
        curColor = 'r';
    end
    
    height = height - 0.8;
end
text(0.5, height+0.3, '$L$', 'interpreter', 'latex', 'Fontsize', 18 ,'horizontalAlignment', 'center', 'color', [0.5, 0.5, 0.5])
plot([0, 1], [height + 0.2, height + 0.2],'color', [0.5, 0.5, 0.5])

xlim([-0.1, 1.1])
ylim([-4, 0.5]);
axis off;
set(gcf, 'color', 'w')