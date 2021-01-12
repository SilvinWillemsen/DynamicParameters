clear all;
interpOrder = 3;
symsVec = {};
for i = 0:interpOrder
    string = "I" + num2str(i);
    
    syms(string);
%     eval("subs(I" + num2str(i) + ", 1)")
    eval("symsVec = [symsVec;" + string + "];");
end

syms alpha xLoc xLocs;

% xLoc = (alpha + 1);
% xLocs = (0:interpOrder) - floor((interpOrder - 1) / 2);
%% custom ip
% xLoc = 2 + 0*alpha;
% xLocs = (0:interpOrder) + 0 * alpha;
% xLocs((length(xLocs) / 2) + 1 : end) = xLocs((length(xLocs) / 2) + 1 : end) + alpha;

%% other ip
% xLoc = ((interpOrder + 1) / 2) - 1 + alpha;
% xLocs = (0:interpOrder) + 0*alpha;
% 
% xLocs((length(xLocs) / 2) + 2 : end) = xLocs((length(xLocs) / 2) + 2 : end) - 1 + alpha;

%% test 7 jan
halfI = (interpOrder + 1) / 2;
xLoc = halfI + 0*alpha;
xLocs = 0:interpOrder + 0*alpha;
xLocs = xLocs + [zeros(1, halfI) * alpha, alpha * ones(1, halfI)];

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