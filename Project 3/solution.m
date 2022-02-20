% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5] = solution()
 [consoleout, A1, A2, A3, A4, A5] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4, A5] = student_solution(dummy_argument)
    % keep these lines to load Fmat and permvec
    % DO NOT submit the mat files to Gradescope
    load Fmat.mat 
    load permvec.mat
    
    % Problem 1
    % ------------------------------------------------------------------------------------------

    n = 8;
    N = n^2;
    L = 10;
    
    del = (L + L)/n;
    
    % For matrix A
    
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
    
    idx = [n N-n];
    
    % construct the diag and index of corresponding diags of B
    diag_B = [e1 -e1 e1 -e1];
    idx_B = [-fliplr(idx) idx];
    
    B = (1/(2*del)) * spdiags(diag_B, idx_B, N, N);
    
    % For matrix C
    
    diag = [e3 -e2];
    idx = [1, n-1];
    
    % construct the diag and index of corresponding diags of C
    diag_C = [-1*rot90(diag, 2) diag];
    idx_C = [-fliplr(idx) idx];
    
    C = (1/(2*del)) * spdiags(diag_C, idx_C, N, N);
    
    A1 = full(A);
    A2 = full(B);
    A3 = full(C);
    
    % Problem 2
    % ------------------------------------------------------------------------------------------
    % This is the center matrix of the encryption
    centerEct = Fmat(161:240, 161:240);
    
    % how many blocks
    blockNum = length(permvec);
    
    % how many blocks for row/column
    num = sqrt(blockNum);
    
    % 20*20
    blockSize = size(centerEct, 1) / num;
    
    % make up the center matrix of the decryption
    centerDct = zeros(size(centerEct));
    
    % get the row and col indexes for permvec
    sz = [num, num];
    [per_row, per_col] = ind2sub(sz, permvec);
    [cen_row, cen_col] = ind2sub(sz, 1:blockNum);
    
    % now construct the decrypt center matrix
    for i = 1:blockNum
        
        % row and col index of given in permvec
        per_rowIdx = (per_row(i)*blockSize - (blockSize-1)) : per_row(i)*blockSize;
        per_colIdx = (per_col(i)*blockSize - (blockSize-1)) : per_col(i)*blockSize;
        
        % row and col index should replace in the decrypt center matrix
        cen_rowIdx = (cen_row(i)*blockSize - (blockSize-1)) : cen_row(i)*blockSize;
        cen_colIdx = (cen_col(i)*blockSize - (blockSize-1)) : cen_col(i)*blockSize;
        
        % update the matrix
        centerDct(cen_rowIdx, cen_colIdx) = centerEct(per_rowIdx, per_colIdx);
        
    end
    
    % Construct the decrypt Fourier matrix
    Dmat = Fmat;
    Dmat(161:240, 161:240) = centerDct;
    
    % shift the permuted matrix, and get the reconstructed image matrix
    Drmat = ifftshift(Dmat);
    reconst = abs(ifft2(Drmat));
    
    % Plot
    set(gcf,'colormap',gray);
    imagesc(uint8(reconst));
    
    A4 = abs(Dmat);
    A5 = reconst;
end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)