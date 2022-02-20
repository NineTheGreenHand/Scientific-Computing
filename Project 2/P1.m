%%% main script %%%

clear all; close all; clc;
K = 1;
tol = 10^-4;                            %Setting the tolerance
L = 4;                                  %Setting the span
x_span = [-L:0.1:L];                    %Setting the span
y_init = 1;                             %Initial value of y
eps_start = 1;                          %Setting the initial shooting value of epsilon
beta_out = zeros(5,1);
j_modes = zeros(5,1000);
y_out = zeros(81, 5);
rhsfunc = @(x, y, eps, K) [y(2); (K*x^2 - eps)*y(1)]; 

for modes = 1:5
    epsilon= eps_start;
    d_epsilon = 0.1;
    
    for j=1:1000
        y_dash_init = sqrt(K*L^2-epsilon)*y_init;
        x0 = [y_init; y_dash_init];
        j_modes(modes,j) = j;
        [t,y] = ode45( @(t,y) rhsfunc(t,y,epsilon,K),x_span,x0);
        y_dash_end = y(end,2);
        y_end = y(end,1);
        check = (y_dash_end+sqrt(K*L^2-epsilon)*y_end);
        if abs(y_dash_end+sqrt(K*L^2-epsilon)*y_end)<tol
            beta_out(modes) = epsilon;
            y_out(:,modes) = y(:,1);
            break;
        elseif (-1)^(modes)*(y_dash_end+sqrt(K*L^2-epsilon)*y_end)>0
            epsilon = epsilon - d_epsilon/2;
            d_epsilon = d_epsilon/2;
        else
            epsilon = epsilon + d_epsilon;
            
        end
    end
    
    eps_start = epsilon+0.1;
    
    plot(t,y(:,1)); hold on;
   
end
close all;
for i = 1:5
    norm(i) = trapz(t, y_out(:,i).^2);
    y_out(:,i) = abs(y_out(:,i)/sqrt(norm(i)));
   % figure(1),plot(t,y_out(:,i)); hold on;
end
A1 = y_out(:,1);
A2 = y_out(:,2);
A3 = y_out(:,3);
A4 = y_out(:,4);
A5 = y_out(:,5);
A6 = beta_out;

