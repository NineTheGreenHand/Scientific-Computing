function f = polfun(t,y,eps)

f1 = y(2);
f2 = eps*(1 - y(1)^2)*y(2) - y(1);

f = [f1; f2]; 