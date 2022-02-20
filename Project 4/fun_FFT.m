function f = fun_FFT(t, w, KX, KY, n, rhsfun) 

w_f = fft2(reshape(w, [n, n]));
psi_f = -w_f./(KX.^2 + KY.^2);
psi = reshape(real(ifft2(psi_f)), [n^2, 1]);
f = rhsfun(psi, w);

end