clear all; close all; clc;

% construct matrix D

n = 25;
e1 = ones(n, 1);

A = spdiags([e1 -2*e1 e1], [-1 0 1], n, n);

A(1,n) = 1; A(n, 1) = 1;

A = full(A);

spy(A)