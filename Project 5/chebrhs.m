function rhs = chebrhs(t, Lap, n, D1, D2, beta, uv)

u = uv(1: n^2);
v = uv((n^2+1):end);

A_sq = u.^2 + v.^2;
lambda = 1 - A_sq;
w = -beta * A_sq;

ut = lambda.*u - w.*v + D1*Lap*u;
vt = w.*u + lambda.*v + D2*Lap*v;

rhs = [ut; vt];

end

