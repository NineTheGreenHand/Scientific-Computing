%% Amath 481 Homework#4 AS COOL AS IT GETS: 
%% I choose to use FFT since it is the fastest.

clear all; close all; clc;

%% First part: Create Matrices A, B, C 

n = 128;
N = n^2;
L = 20;
[A, B, C] = getMatrices(n, N, L);

%% Second Part: Setup other conditions

% Set values to parameters

nu = 0.001;
tspan = 0:0.5:80;

% Set up the space

l = linspace(-L/2, L/2, n+1);
x = l(1:n);
y = x;
[X, Y] = meshgrid(x, y);

% Domain for the FFT method

kx = (2*pi/(2*L))*[0:(n/2-1), (-n/2):-1];
ky = kx;
kx(1) = 10^-6;
ky(1) = 10^-6;
[KX, KY] = meshgrid(kx, ky);

% Define our rhsfun for the system

rhsfun = @(psi, w) (-B*psi).*(C*w) + (C*psi).*(B*w) + nu.*A*w;

%% Third Part: Setup the initial condition

r = 4;
l = 5;

w0 = (exp(-1*(X + r - l).^2 - 10*(Y + l).^2) - exp(-1*(X - r - l).^2 - 10*(Y + r).^2)) + ...
    (exp(-1*(X - l).^2 - 10*(Y + r + l).^2) - exp(-1*(X - l).^2 - 10*(Y - r + l).^2)) - ...
    (exp(-1*(X + r + l).^2 - 10*(Y - l).^2) - exp(-1*(X - r + l).^2 - 10*(Y - l).^2)) - ...
    (exp(-1*(X + l).^2 - 10*(Y + r - l).^2) - exp(-1*(X + l).^2 - 10*(Y - r - l).^2));

w0 = reshape(-w0 + rot90(w0), [N, 1]);

%% Fourth Part: Solve the problem using FFT

[t_ACAIG, w_ACAIG] = ode45(@(t, w) fun_FFT(t, w, KX, KY, n, rhsfun), tspan, w0);

%% Fifth Part: make the 2-D movie of the above simulation 

% Make the movie object first

M(length(t_ACAIG)) = struct('cdata', [], 'colormap', []);
h.Visible = 'off';

for i = 1:length(t_ACAIG)
    
    w_cur = reshape(w_ACAIG(i, :), n, n);
    set(gcf, 'Position',  [500, 100, 1200, 420]);
    
    ax1 = subplot(1, 2, 1);
    surf(X, Y, w_cur); 
    colormap(prism(10));
    shading interp; axis off;
    
    ax2 = subplot(1, 2, 2);
    pcolor(X, Y, w_cur);
    colormap(prism(10));
    shading interp; axis off;
    zlim([-1, 1])
  
    drawnow
    pause(0.001)
    M(i) = getframe(gcf);
    
end

% Write the video file

v = VideoWriter('ACAIG.mp4', 'MPEG-4');
v.FrameRate = 18;
open(v);
writeVideo(v, M);
close(v);














