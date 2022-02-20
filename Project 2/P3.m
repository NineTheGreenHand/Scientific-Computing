b = zeros(99,1); b(1) = 1;%%% main script %%%

clc; clear variables; close all;
K = 1;
tol = 10^-4;                            %Setting the tolerance
L = 2;                                  %Setting the span
x_span = [-L:0.1:L];                    %Setting the span
y_init = 0.1;                           %Initial value of y
eps_start = 1;                          %Setting the initial shooting value of epsilon
gamma1 = 0.05;
gamma2 = -0.05;
beta_out = zeros(2,1);
y_out = zeros(41, 2);
nonlinearfunc = @(x, y, gamma, eps, K) [y(2); (gamma*abs(y(1))^2 + K*x^2 - eps)*y(1)];

for modes = 1:2

    epsilon= eps_start;
    d_epsilon = 0.1;
    
    for j=1:1000
        d_epsilon = 0.1;
        y_dash_init = sqrt(K*L^2-epsilon)*y_init;
        y0 = [y_init; y_dash_init];
        [x,y] = ode45( @(x,y) nonlinearfunc(x,y,gamma1,epsilon,K),x_span,y0);
        norm = trapz(x, y(:,1).^2);
        if abs(norm-1) < tol
            break;
        else
            y_init = y_init/sqrt(norm);
        end
        for k = 1:1000
            y_dash_init = sqrt(K*L^2-epsilon)*y_init;
            y0 = [y_init; y_dash_init];
            [x,y] = ode45( @(x,y) nonlinearfunc(x,y,gamma1,epsilon,K),x_span,y0);
            y_dash_end = y(end,2);
            y_end = y(end,1);
            check = (y_dash_end+sqrt(K*L^2-epsilon)*y_end);
            if abs(check)<tol
                beta_out(modes) = epsilon;
                y_out(:,modes) = abs(y(:,1));
                break;
            elseif (-1)^(modes)*(check)>0
                epsilon = epsilon - d_epsilon/2;
                d_epsilon = d_epsilon/2;
                
            else
                epsilon = epsilon + d_epsilon;
                
            end
            
        end
    end
    
    eps_start = epsilon+0.1;
end

A13 = y_out(:,1);
A14 = y_out(:,2);

A1 = beta_out;