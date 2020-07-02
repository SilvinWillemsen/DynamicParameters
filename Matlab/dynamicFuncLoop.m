range = 0:0.001:4;
results = zeros(length(range), 1);
for i = 1:length(range)
    results(i) = dynamicBoundFunc(30.75, range(i));
end`