clear all
close all

N = 50;
xx = 1:N;
% yy = rand(N, 1);
yy = sin(2 * pi * 10 * (1:N)' / N);

%% Direct-form
syms f(x);
mat = ones(1, N);
for n = 1:N-1
    mat = [mat; xx.^n];
end
mat = mat';
invA = inv(mat);
a = invA * yy;
f = a(1);
for n = 2:N
    f = f + a(n) * x.^(n-1);
end
fplot(f, [xx(1), xx(end)])
figure(1)
hold on;
scatter(xx, yy)

%% Newton
syms f(x);


%% Lagrange
syms f(x)
f = 0;
for i = 1:N
    Li = 1;
    for m = 1:N
        if m == i
            continue
        end
        Li = Li * (x - xx(m)) / (xx(i) - xx(m));
    end
    f = f + Li * yy(i);
end
fplot(f, [xx(1), xx(end)])
legend("Direct", "Data points", "Lagrange")
