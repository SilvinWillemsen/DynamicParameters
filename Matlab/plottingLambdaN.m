clear all
idx = 1;
k = 1/44100;

cRange = 1/(k*10):-0.01:(1/(k*30));
N = zeros(length(cRange), 1);
lambdaSave = zeros(length(cRange), 1);

for c = cRange
    N(idx) = floor(1/(c*k));
    if N == 0
        break;
    end
    h = 1/N(idx);
    lambdaSave(idx) = c*k/h;
    idx = idx + 1;
end
plot(cRange, lambdaSave)
set(gca, 'xdir', 'reverse')
clear all;
k = 1/44100;


NRange = 1:0.01:1000;
idx = 1;
for N = NRange
    cSave(idx) = 1/(N * k);

    h = 1/floor(N);
    lambdaSave(idx) = cSave(idx)*k/h;
    idx = idx + 1;
end
plot(lambdaSave)
