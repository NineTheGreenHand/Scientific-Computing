% FitzHugh Model

function dydt = fitzhugh(t, y, a, b, c, I)

v = y(1);
w = y(2);

dv = -v^3 + (1 + a) * v^2 - a*v - w + I;
dw = b*v - c*w;
dydt = [dv; dw];

