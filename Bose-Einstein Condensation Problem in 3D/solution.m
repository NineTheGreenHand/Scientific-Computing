% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4] = solution()
 [consoleout, A1, A2, A3, A4] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4] = student_solution(dummy_argument)
    
    % Set up
    % -------------------------------------------------------------------------
    % List the given conditions

    L = 2*pi;
    n = 16;
    tspan = 0:0.5:4;
    A = [-1, -1, -1];
    B = -A;

    % Set up Fourier Space

    kxyz = (2*pi/L)*[0:(n/2-1), (-n/2):-1];
    kxyz(1) = 10^-6;
    [KX, KY, KZ] = meshgrid(kxyz, kxyz, kxyz);

    % Set up X, Y, Z grid

    xyz2 = linspace(-L/2, L/2, n+1);
    xyz = xyz2(1:n);
    [X, Y, Z] = meshgrid(xyz, xyz, xyz);

    % Set Laplacian

    K = KX.^2 + KY.^2 + KZ.^2;
    Lap = -K;

    % Part A: Initial Condition as cos(x)cos(y)cos(z)
    % -------------------------------------------------------------------------
    % Set up initial condition 

    psi0 = cos(X).*cos(Y).*cos(Z);
    psif0 = reshape(fftn(psi0), [n^3, 1]);

    % Solve the system

    [tA, psifA] = ode45(@(t, psif) rhs(t, psif, A, B, X, Y, Z, n, Lap), tspan, psif0);

    A1 = real(psifA);
    A2 = imag(psifA);

    % Part B: Initial Condition as sin(x)sin(y)sin(z)
    % -------------------------------------------------------------------------
    % Set up initial condition 

    psi0 = sin(X).*sin(Y).*sin(Z);
    psif0 = reshape(fftn(psi0), [n^3, 1]);

    % Solve the system

    [tB, psifB] = ode45(@(t, psif) rhs(t, psif, A, B, X, Y, Z, n, Lap), tspan, psif0);

    A3 = real(psifB);
    A4 = imag(psifB);
    
    % Visualize using isosurface
    % -------------------------------------------------------------------------

    % psif = psifA;
    % 
    % for i = 1:size(psif, 1)
    %     
    %     f_cur = reshape(psif(i, :), [n, n, n]);
    %     cur = ifftn(f_cur);
    %     abscur = cur.*conj(cur);
    %     
    %     isosurface(X, Y, Z, abscur, 0.5);
    %     colormap(jet(9));
    %     axis('square');
    %     
    % end

    % psif = psifB;
    % 
    % for i = 1:size(psif, 1)
    %     
    %     f_cur = reshape(psif(i, :), [n, n, n]);
    %     cur = ifftn(f_cur);
    %     abscur = cur.*conj(cur);
    %     
    %     isosurface(X, Y, Z, abscur, 0.5);
    %     colormap(jet(9));
    %     axis('square');
    %     
    % end

end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)