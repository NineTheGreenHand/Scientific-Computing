%Homework 2 Problem 2
clc; clear variables; close all;
%Setting parameters like delta_x, K and L
x = [-4:0.1:4];
delta_x = 0.1;
K = 1;
%Creating a matrix of 81 rows 81 colums
A = zeros(79,79);

%Creating first and last row
A(1,1) = 2/3 + (delta_x^2)*K*(x(2)^2);
A(1,2) = -2/3;
A(79,78) = -2/3;
A(79,79) = 2/3 + (delta_x^2)*K*(x(80)^2);

%Loop for all the central difference functions
for i = 2:78
    A(i,i-1) = -1;
    A(i,i) = 2+((delta_x^2)*K*(x(i+1)^2));
    A(i,i+1) = -1;
end

%Solve for eigen values and eigen vectors
[y, epsilon] = eig(A);

%Since we solved the equation Ax = (dx^2*epsilon) we have to divide by dx^2
epsilon = epsilon/(delta_x^2);

%Getting a vector of eigen values
eig_values = epsilon(epsilon>0);

%Sorting the eigen values
eig_values_sort = sort(epsilon(epsilon>0));
lowest5_eig = eig_values_sort(1:5);

for i = 1:5
    index_of_lowest_eig_values = (epsilon == lowest5_eig(i));
    [row, col(i)] = find(index_of_lowest_eig_values);
end
lowest_eig_vectors = y(:,col);
y0 = zeros(1,5);
for i = 1:5
    y0(i) = (4*lowest_eig_vectors(1,i) - lowest_eig_vectors(2,i))/(3+2*delta_x*sqrt(K*16-lowest5_eig(i)));
    yn(i) = (-4*lowest_eig_vectors(79,i) + lowest_eig_vectors(78,i))/(3+2*delta_x*sqrt(K*16-lowest5_eig(i)));
end
full_eigen_vectors = [y0; lowest_eig_vectors; yn];
for i = 1:5
    norm(i) = trapz(x, full_eigen_vectors(:,i).^2);
    final_out(:,i) = abs(full_eigen_vectors(:,i)/sqrt(norm(i)));
    figure(2),plot(x,final_out(:,i)); hold on;
end



A7 = final_out(:,1);

A8 = final_out(:,2);

A9 = final_out(:,3);

A10 = final_out(:,4);

A11 = final_out(:,5);

A12 = lowest5_eig;
