clear all;
interpOrder = 12;
symsVec = {};
for i = 0:interpOrder
    string = "I" + num2str(i);
    
    syms(string);
%     eval("subs(I" + num2str(i) + ", 1)")
    eval("symsVec = [symsVec;" + string + "];");
end

syms alf xLoc xLocs;

% xLoc = (alf + 1);
% xLocs = (0:interpOrder) - floor((interpOrder - 1) / 2);
%% custom ip
% xLoc = 2 + 0*alf;
% xLocs = (0:interpOrder) + 0 * alf;
% xLocs((length(xLocs) / 2) + 1 : end) = xLocs((length(xLocs) / 2) + 1 : end) + alf;

%% other ip
% xLoc = ((interpOrder + 1) / 2) - 1 + alf;
% xLocs = (0:interpOrder) + 0*alf;
% 
% xLocs((length(xLocs) / 2) + 2 : end) = xLocs((length(xLocs) / 2) + 2 : end) - 1 + alf;

%% test 7 jan
% halfI = (interpOrder + 1) / 2;
% xLoc = halfI + 0*alf;
% xLocs = 0:interpOrder + 0*alf;
% xLocs = xLocs + [zeros(1, halfI) * alf, alf * ones(1, halfI)];

if mod(interpOrder,2) == 0
    %% even
    xLoc = interpOrder*0.5 + 0*alf;
    xLocs = (0:interpOrder) + 0*alf;
    xLocs = xLocs + [zeros(1, interpOrder * 0.5), ones(1, interpOrder * 0.5 + 1) * (alf - 1)];
else
    xLoc = 2 + 0*alf;
    xLocs = 0:interpOrder + 0*alf;    
    xLocs = xLocs + [zeros(1, ceil(interpOrder * 0.5)), ones(1, ceil(interpOrder * 0.5)) * (alf - 1)];
end
for i = 1:interpOrder+1
    for k = 1 : interpOrder+1
        if k ~= i
            symsVec(i) = symsVec(i) * (xLoc - xLocs(k)) / (xLocs(i) - xLocs(k));
        end
    end
    eval("symsVec(i) = symsVec(i) / I" + num2str(i-1) + ";")
end
for alf = 0:0.01:1
    plot(subs(symsVec,alf))
    drawnow
    pause(0.1);
end