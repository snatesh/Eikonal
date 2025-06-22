function [Vm, Hm] = arnoldi_iteration(A, b, tol, m_max)
% ARNOLDI_ITERATION Computes orthonormal basis Vm and Hessenberg matrix Hm
% for the Krylov subspace K_m(A,b) using Arnoldi iteration.
%
% Inputs:
%   A     : n x n matrix (assumed sparse or large)
%   b     : n x 1 vector (starting vector for Krylov subspace)
%   tol   : scalar tolerance for termination based on residual norm
%   m_max : (optional) maximum number of iterations (default: length(b))
%
% Outputs:
%   Vm    : n x m orthonormal basis matrix for Krylov subspace
%   Hm    : m x m upper Hessenberg matrix such that A*Vm ≈ Vm*Hm

    if nargin < 4
        m_max = length(b);
    end

    n = length(b);
    Vm = zeros(n, m_max + 1);       % Preallocate with one extra column
    Hm = zeros(m_max + 1, m_max);   % One extra row

    % Normalize b to get the first basis vector
    beta = norm(b);
    if beta < eps
        error('Input vector b must be non-zero.');
    end
    Vm(:,1) = b / beta;

    for j = 1:m_max
        w = A * Vm(:,j);  % Apply matrix A

        % Orthogonalize against previous vectors
        for i = 1:j
            Hm(i,j) = Vm(:,i)' * w;
            w = w - Hm(i,j) * Vm(:,i);
        end

        % Compute norm and check convergence
        Hm(j+1,j) = norm(w);
        if Hm(j+1,j) < tol
            break;
        end

        % Normalize and store next basis vector
        Vm(:,j+1) = w / Hm(j+1,j);
    end

    % Trim output
    m = j;
    Vm = Vm(:,1:m);        % size: n x m
    Hm = Hm(1:m,1:m);      % size: m x m
end
