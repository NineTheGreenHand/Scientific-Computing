function rhs = rhs(t, psif, A, B, X, Y, Z, n, Lap)

% Define psi in fourier space, and psi. 

psif = reshape(psif, [n, n, n]);
psi = ifftn(psif);

% Linear part:

linf = (1/2).*Lap.*psif;

% Nonlinear part:

nlin_3D = (A(1).*sin(X).^2 + B(1)).*(A(2).*sin(Y).^2 + B(2)).*(A(3).*sin(Z).^2 + B(3));
nlin = (-psi.*conj(psi) + nlin_3D).*psi;
nlinf = fftn(nlin);

% Combine linear and nonlinear part to get rhs = psift

rhs = reshape((i*(linf + nlinf)), [n^3, 1]);

end