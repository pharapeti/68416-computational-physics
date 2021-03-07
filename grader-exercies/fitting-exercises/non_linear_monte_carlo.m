clear; clc; close all;

x = linspace(0,(1+rand(1)))*5;
y = (1+rand(1))*exp(-(1+rand(1))*x) + randn([10+randi(10,1), length(x)])*0.05;

% Model function
f = @(p) p(1) .* exp(-p(2) .* x);

aVals = nan(1, length(y));
kVals = nan(1, length(y));

% loop over the rows in y
for r = 1:length(y)
    % Create error function with selected row of y
    m = @(p) norm(y(r) - f(p));

    % Minimise the error function to find A and k
    p = fminsearch(m, [1,0]);

    % store the result in an array
    A = p(1);
    k = p(2);

    aVals(1, r) = A;
    kVals(1, r) = k;
end

% average of repeats
mA= mean(aVals);
mk = mean(kVals);

% % standard deviation of repeats
sA = std(aVals);
sk = std(kVals);
