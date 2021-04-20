% test 1
disp("Test 1 - small matrix")
A1 = [4 1 0 0; 1 5 -1 0; 0 -1 6 3; 0 0 3 7];
M1 = [0 4 1; 1 5 -1; -1 6 3; 3 7 0]
b1 = [1 2 3 4]';
disp("Expected x:")
x1 = A1 \ b1
x10 = [0 0 0 0]'
omega1 = 1.1
acc1 = 1e-13;
[x1, j1, g1, s1] = iter3(M1, b1, x10, acc1, omega1)

%f1 = @(o) iter3(M1, b1, x10, acc1, o);
%[allOmega1, allS1] = measure_omega(f1, 0.1, 30, 0.05);

%figure(1)
%plot(allOmega1, allS1);
%title("Example 1");
%xlabel("Values of omega");
%ylabel("Number of iterations");

% test 2
disp("Test 2 - big Laplace matrix")
n = 10^6;
M2 = zeros(n, 3);
M2(:, 1) = -1;
M2(:, 2) = 2;
M2(:, 3) = -1;
M2(1, 1) = 0;
M2(end, end) = 0;

b2 = zeros(n, 1);
b2(1) = 1;
b2(end) = 1;

x20 = zeros(n, 1);

omega2 = 1;

tic
[x2, j2, g2, s2] = iter3(M2, b2, x20, 1e-3, 1.7);
toc
j2
g2
s2

%f2 = @(o) iter3(M2, b2, x20, 1e-3, o);
%[allOmega2, allS2] = measure_omega(f2, 0.5, 14, 0.1);

%figure(2)
%plot(allOmega2, allS2);
%title("Example 2");
%xlabel("Values of omega");
%ylabel("Number of iterations");


function [allOmega, allS] = measure_omega(f, startValue, numIterations, increment)
    omega = startValue;
    allOmega = [] * numIterations;
    allS = [] * numIterations;
    for i = 1 : numIterations
        allOmega(i) = omega;
        [~, ~, ~, s] = f(omega);
        allS(i) = s;
        omega = omega + increment;
    end
end
