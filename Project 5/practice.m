clear all; close all; clc;

%% Part (a): Periodic Boundary Conditions
% List given conditions
% -------------------------------------------------------------------------

L = 20;
n = 64;
beta = 1;
D1 = 0.1;
D2 = D1;             
tspan = 0:0.5:4;

% Set up the initial conditions 
% -------------------------------------------------------------------------

m = 1;  % number of spirals
x2 = linspace(-L/2, L/2, n+1);
x = x2(1:n);
y = x;
[X, Y] = meshgrid(x, y);

u0 = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
v0 = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
uf0 = reshape(fft2(u0), [n^2, 1]);
vf0 = reshape(fft2(v0), [n^2, 1]);
uvf0 = [uf0; vf0];

kx = (2*pi/L)*[0:(n/2-1), (-n/2):-1];
ky = kx;
kx(1) = 10^-6;
ky(1) = 10^-6;
[KX, KY] = meshgrid(kx, ky);
K = KX.^2 + KY.^2;
Lap_fft = -K;   % Laplacian

% Solve for the system
% -------------------------------------------------------------------------

[t_FFT, uvf_FFT] = ode45(@(t, uvf) fftrhs(t, Lap_fft, n, D1, D2, beta, uvf), tspan, uvf0);

A1 = real(uvf_FFT);
A2 = imag(uvf_FFT);

%% Part (b): Dirichlet Boundary Conditions
% List given conditions
% -------------------------------------------------------------------------

L = 20;
n = 30;
beta = 1;
D1 = 0.1;
D2 = D1;
m = 1;
tspan = 0:0.5:4;

% Set up initial conditions
% -------------------------------------------------------------------------

n = 30;
[D, x] = cheb(n);

D = (2/L).*D;
x = (L/2).*x;

D_sq = D^2;
D_sq(1, :) = zeros(1, n+1);
D_sq(n+1, :) = zeros(1, n+1);

y = x;
[X, Y] = meshgrid(x, y);

u0 = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
v0 = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2))); 

uinitvec = reshape(u0, [(n+1)^2, 1]);
vinitvec = reshape(v0, [(n+1)^2, 1]);
uv0 = [uinitvec; vinitvec];

I = eye(length(D_sq));
Lap_cheb = kron(D_sq, I) + kron(I, D_sq);   % Laplacian

% Solve for the system
% -------------------------------------------------------------------------

[t_Cheb, uv_Cheb] = ode45(@(t, uv) chebrhs(t, Lap_cheb, n+1, D1, D2, beta, uv), tspan, uv0);

A3 = uv_Cheb;







