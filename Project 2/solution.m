% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = solution()
 [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = student_solution(dummy_argument)
    
    % Problem 1
    % ------------------------------------------------------------------------------------------
    tol = 10^(-4);
    L = 4;
    K = 1;
    A = 1;                  % phi(-L) = A;
    xspan = -L:0.1:L;

    eigenVal = zeros(5, 1);
    eigenFun = zeros(length(xspan), 5);

    % define function
    phiFun = @(x, y, K, eps) [y(2); (K*x^2 - eps)*y(1)];

    iters = ones(5, 1).*(-1);

    eps_start = 1;

    for modes = 1:5
    
        eps = eps_start;
        deps = eps_start/100;
    
        for j = 1:1000
        
            y0 = [A, A * sqrt(K*L^2-eps)];
            [x, y] = ode45(@(x, y) phiFun(x, y, K, eps), xspan, y0);
        
            % the value to check
        
            check = y(end, 2) + sqrt(K*L^2 - eps)*y(end, 1);
        
            if abs(check) < tol
                iters(modes) = j;
                break
            end
        
            if (-1) ^ (modes) * check > 0
                eps = eps - deps/2;
                deps = deps/2;
            else
                eps = eps + deps;
            end
        end
    
        eigenVal(modes) = eps;
        eigenFun(:, modes) = abs(y(:, 1) / sqrt(trapz(x, y(:, 1).^2)));
    
        eps_start = eps + 0.1;
    
        % plot(x,y(:,1)); hold on;
    end

    A1 = eigenFun(:, 1);
    A2 = eigenFun(:, 2);
    A3 = eigenFun(:, 3);
    A4 = eigenFun(:, 4);
    A5 = eigenFun(:, 5);
    A6 = eigenVal;
    
    %------------------------------------------------------------------------------------------S
    % Problem 2
    
    x = -4:0.1:4;
    dx = 0.1;
    K = 1;
    A = zeros(79, 79);
    A(1, 1) = 2/3 + (dx^2)*K*(x(2)^2);
    A(1, 2) = -2/3;
    A(end, end - 1) = -2/3;
    A(end, end) = 2/3 + (dx^2)*K*(x(end - 1)^2);

    % Get the inner matrix part
    for i = 2:78
        A(i,i-1) = -1;
        A(i,i+1) = -1;
        A(i,i) = 2 + ((dx^2)*K*(x(i+1)^2));
    end

    %Solve for eigenvalues and eigenvectors
    [y, eps] = eig(A);
    eps = eps/(dx^2);

    % Get a vector of eigenvalues
    eig_values = eps(eps>0);

    % Get the first five eigenvalues in order
    %Sorting the eigen values
    eigSort = sort(eig_values);
    lowestEigenval = eigSort(1:5);

    % Now working on finding eigenvectors
    for i = 1:5

        % fine the index of columns of the eigenvalue to get eigenvectors
        [row, col(i)] = find(eps == lowestEigenval(i));

    end

    % fine the first five eigenvectors 
    lowestEigenvec = y(:,col);

    % now find the first and last element of the complete version of eigenvectors
    y0 = zeros(1, 5);
    yn = zeros(1, 5);

    for i = 1:5

        y0(i) = (4*lowestEigenvec(1,i) - lowestEigenvec(2,i))/(3+2*dx*sqrt(16*K-lowestEigenval(i)));
        yn(i) = (4*lowestEigenvec(79,i) - lowestEigenvec(78,i))/(3+2*dx*sqrt(16*K-lowestEigenval(i)));

    end

    % combine the first row and last row of the eigenvectors
    eigenVec = [y0; lowestEigenvec; yn];

    % Normalize the eigenfunctions 
    for i = 1:5

        norm(i) = trapz(x, eigenVec(:,i).^2);
        finalVec(:,i) = abs(eigenVec(:,i)/sqrt(norm(i)));
        % plot(x,finalVec(:,i)); hold on;

    end

    A7 = finalVec(:, 1);
    A8 = finalVec(:, 2);
    A9 = finalVec(:, 3);
    A10 = finalVec(:, 4);
    A11 = finalVec(:, 5);
    A12 = lowestEigenval;
    
    %------------------------------------------------------------------------------------------S
    % Problem 3
    
    % Problem 3
    % Setup all values and function

    K = 1;
    L = 2;
    tol = 10^-4;
    xspan = -L:0.1:L;
    eps_start = 1;
    A = 0.1;
    gamma1 = 0.05;
    gamma2 = -0.05;
    eigenV1 = zeros(2, 1);
    eigenV2 = zeros(2, 1);
    eigenF1 = zeros(41, 2);
    eigenF2 = zeros(41, 2);
    nonlinearFun = @(x, y, gamma, eps, K) [y(2); (gamma*abs(y(1))^2 + K*x^2 - eps)*y(1)];

    % for gamma = 0.05

    for modes = 1:2

        eps = eps_start;

        for j = 1:1000

            deps = 0.1;
            y0 = [A; A*sqrt(K*L^2-eps)];
            [x, y] = ode45(@(x, y) nonlinearFun(x, y, gamma1, eps, K), xspan, y0);
            norm = trapz(x, (y(:, 1).^2));
            if abs(norm - 1) < tol
                break
            else
                A = A/sqrt(norm);
            end
            % start shooting scheme for epsilon
            for i = 1:1000

                y0 = [A; A*sqrt(K*L^2-eps)];
                [x, y] = ode45(@(x, y) nonlinearFun(x, y, gamma1, eps, K), xspan, y0);

                check = y(end, 2) + sqrt(K*L^2 - eps)*y(end, 1);
                if abs(check) < tol
                    break
                end

                if (-1) ^ (modes) * check > 0
                    eps = eps - deps/2;
                    deps = deps/2;
                else
                    eps = eps + deps;
                end
            end
        end

        eigenV1(modes) = eps;
        eigenF1(:,modes) = abs(y(:,1));
        eps_start = eps + 0.1;
    end

    A13 = eigenF1(:, 1);
    A14 = eigenF1(:, 2);
    A15 = eigenV1;

    % for gamma = -0.05

    eps_start = 1;
    A = 0.1; 

    for modes = 1:2

        eps = eps_start;

        for j = 1:1000

            deps = 0.1;
            y0 = [A; A*sqrt(K*L^2-eps)];
            [x, y] = ode45(@(x, y) nonlinearFun(x, y, gamma2, eps, K), xspan, y0);
            norm = trapz(x, (y(:, 1).^2));
            if abs(norm - 1) < tol
                break
            else
                A = A/sqrt(norm);
            end
            % start shooting scheme for epsilon
            for i = 1:1000

                y0 = [A; A*sqrt(K*L^2-eps)];
                [x, y] = ode45(@(x, y) nonlinearFun(x, y, gamma2, eps, K), xspan, y0);

                check = y(end, 2) + sqrt(K*L^2 - eps)*y(end, 1);
                if abs(check) < tol
                    break
                end

                if (-1) ^ (modes) * check > 0
                    eps = eps - deps/2;
                    deps = deps/2;
                else
                    eps = eps + deps;
                end
            end
        end

        eigenV2(modes) = eps;
        eigenF2(:,modes) = abs(y(:,1));
        eps_start = eps + 0.1;
    end

    A16 = eigenF2(:, 1);
    A17 = eigenF2(:, 2);
    A18 = eigenV2;
end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)