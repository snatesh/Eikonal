function eik_path = eikonal_extract_path(n, cu_glb, x0, goal, DT, tol)

a = 0.5; b = 0.5; c = 0.5;
meshT = DT.ConnectivityList;
meshP = DT.Points';
 
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);

% derivative matrices
Dx_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);

rhs = @(t,x) -normalize_grad(evaluate_grad_u(x(:), cu_glb, meshT, meshP, DT, ...
                             Dx_a1bc1, Dy_ab1c1, H_a1bc1, H_ab1c1, n, H_abc) );
% ODE options
opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-2, 'MaxStep', 0.1);

% Integrate
[t, X] = ode45(rhs, [0, 10], x0(:), opts);  % time horizon can be adjusted

% Stop early if close enough to goal
dists = vecnorm(X - goal(:)', 2, 2);
cut_idx = find(dists < tol, 1);
if isempty(cut_idx)
    eik_path = X;
else
    eik_path = X(1:cut_idx, :);
end


end

function v = normalize_grad(grad)

if norm(grad) < 1e-10
  v = [0; 0];  % freeze eik_path
else
  v = grad / norm(grad);
end

end
