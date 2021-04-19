% test 1
disp("Test 1")
M = [0 2 -1; -1 2 -1; -1 2 -1; -1 2 0]
b = [1 0 0 1]'
x0 = [0 0 0 0]'
omega = 2
acc = 1e-6;
[x, j, g, s] = iter3(M, b, x0, acc, omega)

omega = 0.1;
allOmega = [] * 100;
allS = [] * 100;
for i = 1 : 100;
    allOmega(i) = omega;
    [~, ~, ~, s] = iter3(M, b, x0, acc, omega);
    allS(i) = s;
    omega = omega + 0.05;
end

figure(1)
plot(allOmega, allS);

allJ = [] * 15;
allG = [] * 15;
allS = [] * 15;
allAcc = [] * 15;
for i = 1 : 15
    acc = 10^(-i);
    [~, j, g, s] = iter3(M, b, x0, acc, 1.5);
    allJ(i) = j;
    allG(i) = g;
    allS(i) = s;
    allAcc(i) = acc;
end

figure(2)
plot(allAcc, allJ, 'r');
hold on
plot(allAcc, allG, 'g');
plot(allAcc, allS, 'b');