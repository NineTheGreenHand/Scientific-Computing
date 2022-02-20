% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5] = solution()
 [consoleout, A1, A2, A3, A4, A5] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5] = student_solution(dummy_argument)

    % First part: Create Matrices A, B, C
    % =====================================================================

    n = 64;
    N = n^2;
    L = 20;
    [A, B, C] = getMatrices(n, N, L);
    d = L/n;
    A(1, 1) = 2*(1/d^2);

    % Second Part: Setup other conditions
    % =====================================================================

    % Given conditions

    nu = 0.001;
    tspan = 0:0.5:4;

    % Define our initial ocndition

    x2 = linspace(-L/2, L/2, n+1);
    x = x2(1:n);
    y = x;
    [X, Y] = meshgrid(x, y);
    w0 = reshape(exp(-(X.^2)-((Y.^2)/20)), [N, 1]);

    % Define our rhsfun for the system

    rhsfun = @(psi, w) (-B*psi).*(C*w) + (C*psi).*(B*w) + nu.*A*w;

    % Third Part: Solve the problem
    % =====================================================================

    % For A\b:
    % ---------------------------------------------------------------------

    fun_BS = @(t, w, A, rhsfun) rhsfun(A\w, w);

    tic;
    [t_BS, w_BS] = ode45(@(t, w) fun_BS(t, w, A, rhsfun), tspan, w0);
    time_BS = toc;

    % For LU:
    % ---------------------------------------------------------------------

    [L, U, P] = lu(A);
    fun_LU = @(t, w, L, U, P, rhsfun) rhsfun(U\(L\(P*w)), w);

    tic;
    [t_LU, w_LU] = ode45(@(t, w) fun_LU(t, w, L, U, P, rhsfun), tspan, w0);
    time_LU = toc;

    % For bicgstab:
    % ---------------------------------------------------------------------

    fun_BIC = @(t, w, A, rhsfun) rhsfun(bicgstab(A, w, 10^-6, 500), w);

    tic;
    [t_BIC, w_BIC] = ode45(@(t, w) fun_BIC(t, w, A, rhsfun), tspan, w0);
    time_BIC = toc;

    % For gmres:
    % ---------------------------------------------------------------------

    restart = 100;
    fun_GM = @(t, w, A, rhsfun) rhsfun(gmres(A, w, restart, 10^-6), w);

    tic;
    [t_GM, w_GM] = ode45(@(t, w) fun_GM(t, w, A, rhsfun), tspan, w0);
    time_GM = toc;

    % For FFT:
    % ---------------------------------------------------------------------
    % Domain for the FFT method

    L = 20;
    kx = (2*pi/L)*[0:(n/2-1), (-n/2):-1];
    ky = kx;
    kx(1) = 10^-6;
    ky(1) = 10^-6;
    [KX, KY] = meshgrid(kx, ky);

    tic;
    [t_FFT, w_FFT] = ode45(@(t, w) fun_FFT(t, w, KX, KY, n, rhsfun), tspan, w0);
    time_FFT = toc;

    % Part Four: Store A1 - A5
    % =====================================================================

    A1 = w_BS;
    A2 = w_LU;
    A3 = w_BIC;
    A4 = w_GM;
    A5 = w_FFT;

end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)