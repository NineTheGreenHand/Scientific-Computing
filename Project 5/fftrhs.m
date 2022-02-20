function rhs = fftrhs(t, Lap, n, D1, D2, beta, uvf)

uf = reshape(uvf(1:n^2), [n, n]);
vf = reshape(uvf((n^2+1):end), [n, n]);

u = reshape(real(ifft2(uf)), [n^2, 1]);
v = reshape(real(ifft2(vf)), [n^2, 1]);

A_sq = u.^2 + v.^2;
lambda = 1 - A_sq;
w = -beta * A_sq;

utf_nonlinear = fft2(reshape((lambda.*u - w.*v), [n, n]));
vtf_nonlinear = fft2(reshape((w.*u + lambda.*v), [n, n]));

utf = reshape((utf_nonlinear + D1*Lap.*uf), [n^2, 1]);
vtf = reshape((vtf_nonlinear + D2*Lap.*vf), [n^2, 1]);

rhs = [utf; vtf];

end