function [u, reactions] = solve_with_dirichlet(K, bc_idx, bc_val, f)
% K        : (ndof x ndof) global stiffness
% bc_idx   : column vector of constrained DOF indices
% bc_val   : column vector of prescribed displacement values at those DOFs
% f        : forces (ndof x 1)
%
% RETURNS
% u        : full displacement vector (nodes x 6)
% reactions: reaction forces at constrained DOFs (|bc_idx| x 1)

    ndof = size(K,1);
    if nargin < 4 || isempty(f), f = zeros(ndof,1); end
    if ~iscolumn(bc_idx), bc_idx = bc_idx(:); end
    if ~iscolumn(bc_val), bc_val = bc_val(:); end

    % Make index sets
    all = (1:ndof).';
    free_idx = setdiff(all, bc_idx);

    % Partition
    Kff = K(free_idx, free_idx);
    Kfc = K(free_idx, bc_idx);
    Kcf = K(bc_idx, free_idx);
    Kcc = K(bc_idx, bc_idx);

    ff = f(free_idx);
    fc = f(bc_idx);

    % Solve reduced system: Kff * uf = ff - Kfc * uc
    uf = Kff \ (ff - Kfc * bc_val);

    % Assemble full u
    u = zeros(ndof,1);
    u(free_idx) = uf;
    u(bc_idx)   = bc_val;
    u = reshape(u', 6, [])';

    % Reactions at constrained DOFs: rc = Kcf*uf + Kcc*uc - fc
    reactions = Kcf * uf + Kcc * bc_val - fc;
end