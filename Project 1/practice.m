clear all; close all; clc

% Problem 1
y = @(t) pi*exp(3*(cos(t) - 1))/sqrt(2);
f = @(t, y) -3*y*sin(t);
y_true(1) = y(0);
y_num(1) = y(0);
dt = 2.^(-[2:8]);
eidx = 1;
for i = dt
    t = [0:i:5];
    lg = length(t);
    for j = 2:lg
        y_true(j) = y(t(j));
        y_num(j) = y_num(j - 1) + i * f(t(j - 1), y_num(j - 1));
    end
    E(eidx) = mean(abs(y_true - y_num));
    eidx = eidx + 1;
end
figure; plot(log(dt), log(E));
p = polyfit(log(dt), log(E),1);
% p_v = polyval(p, log(dt));
% plot(log(dt), p_v);
A1 = y_num';
A2 = E;
A3 = p(1);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% y = @(t) pi*exp(3*(cos(t) - 1))/sqrt(2);
% f = @(t, y) -3*y*sin(t);
% y_true_2(1) = y(0);
% y_num_2(1) = y(0);
% dt = 2.^(-[2:8]);
% eidx = 1;
% for i = dt
%     t = [0:i:5];
%     lg = length(t);
%     for j = 2:lg
%         y_true_2(j) = y(t(j));
%         y_num_2(j) = y_num_2(j-1)+(i/2)*(f(t(j-1), y_num_2(j-1)) + f((t(j-1)+i), (y_num_2(j-1)+i*f(t(j-1), y_num_2(j-1)))));
%     end
%     E(eidx) = mean(abs(y_true_2 - y_num_2));
%     eidx = eidx + 1;
% end
% figure; plot(log(dt), log(E));
% p = polyfit(log(dt), log(E),1);
% % p_v = polyval(p, log(dt));
% % plot(log(dt), p_v);
% A4 = y_num_2';
% A5 = E;
% A6 = p(1);

% Problem 2

% eps = [0.1, 1, 20];
% tp = [0:0.5:32];  
% y0 = [sqrt(3); 1];
% index = 1;
% for modes = eps
%     [t, y] = ode45(@(t,y) polfun(t, y, modes), tp, y0);
%      ans_1(:, index) = y(:, 1);
%      index = index + 1;
% end

% tp = [0, 32];
% eps = 1;
% y0 = [2; pi^2];
% TOL = 10.^(-[4:10]);
% index = 1;
% for i = TOL
%     options = odeset('RelTol', i, 'AbsTol', i);
%     [T45, Y45] = ode45(@(t,y) polfun(t, y, eps), tp, y0, options);
%     [T23, Y23] = ode23(@(t,y) polfun(t, y, eps), tp, y0, options);
%     [T113, Y113] = ode113(@(t,y) polfun(t, y, eps), tp, y0, options);
%     average45(index) = mean(diff(T45));
%     average23(index) = mean(diff(T23));
%     average113(index) = mean(diff(T113));
%     index = index + 1;
% end
% figure; plot(log(average45), log(TOL));
% figure; plot(log(average23), log(TOL));
% figure; plot(log(average113), log(TOL));
% p45 = polyfit(log(average45), log(TOL),1);
% p23 = polyfit(log(average23), log(TOL),1);
% p113 = polyfit(log(average113), log(TOL),1);
% A8 = p45(1);
% A9 = p23(1);
% A10 = p113(1);

% Problem 3
% 
% a1 = 0.05;
% a2 = 0.25;
% b = 0.01;
% c = 0.01;
% I = 0.1;
% t = [0:0.5:100];
% d12 = 0; d21 = 0;
% init = [0.1, 0.1, 0, 0];
% [t, y1] = ode15s(@(t, Y) fizint(t, Y, a1, a2, b, c, I, d12, d21), t, init);
% d12 = 0; d21 = 0.2;
% [t, y2] = ode15s(@(t, Y) fizint(t, Y, a1, a2, b, c, I, d12, d21), t, init);
% d12 = -0.1; d21 = 0.2;
% [t, y3] = ode15s(@(t, Y) fizint(t, Y, a1, a2, b, c, I, d12, d21), t, init);
% d12 = -0.3; d21 = 0.2;
% [t, y4] = ode15s(@(t, Y) fizint(t, Y, a1, a2, b, c, I, d12, d21), t, init);
% d12 = -0.5; d21 = 0.2;
% [t, y5] = ode15s(@(t, Y) fizint(t, Y, a1, a2, b, c, I, d12, d21), t, init);
% A11 = y1; 
% A12 = y2;
% A13 = y3;
% A14 = y4;
% A15 = y5;
















