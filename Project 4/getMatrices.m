function [A, B, C] = getMatrices(n, N, L)

del = L/n;

% For matrix A
% -------------------------------------------------------------------------

e0 = zeros(N, 1);
e1 = ones(N, 1);

% set index 8, 16, 24 as 1s, others keep to be 0s
e2 = e0;
e2(mod(1:N, n) == 0) = 1;


% set index 1, 9, 17... as 0s, others keep to be 1s. 
e3 = e1;
e3(mod(0:N-1, n) == 0) = 0;

diag = [e3 e2 e1 e1];
idx = [1 n-1 n N-n];

% construct the diag and index of corresponding diags of A
diag_A = [rot90(diag, 2) -4*e1 diag];
idx_A = [-fliplr(idx) 0 idx];

A = (1/del^2) * spdiags(diag_A, idx_A, N, N);

% For matrix B
% -------------------------------------------------------------------------

idx = [n N-n];

% construct the diag and index of corresponding diags of B
diag_B = [e1 -e1 e1 -e1];
idx_B = [-fliplr(idx) idx];

B = (1/(2*del)) * spdiags(diag_B, idx_B, N, N);

% For matrix C
% -------------------------------------------------------------------------

diag = [e3 -e2];
idx = [1, n-1];

% construct the diag and index of corresponding diags of C
diag_C = [-1*rot90(diag, 2) diag];
idx_C = [-fliplr(idx) idx];

C = (1/(2*del)) * spdiags(diag_C, idx_C, N, N);

end