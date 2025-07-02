function [g,iTri,tx] = evaluate_grad_u(x, cu_glb, meshT, meshP, DT, Dx, Dy, Hx, Hy, n, varargin)

plot(x(1),x(2),'r.'); hold on; drawnow;
if any(~isfinite(x))
    warning('x is non-finite: [%g, %g]', x(1), x(2));
end

a = 0.5; b = 0.5; c = 0.5;
M = n*(n+1)/2;

iTri = pointLocation(DT, x');
if isnan(iTri)
  for k = 1:size(meshT,1)
    P = meshP(:,meshT(k,:));  % 2x3: triangle vertices
    A = IncidenceMatrix(P);
    rs = A \ (x - P(:,1));
    r = rs(1); s = rs(2);
    if r >= -1e-12 && s >= -1e-12 && r + s <= 1 + 1e-12
      iTri = k;
      break;
    end
  end
end

if isnan(iTri)
  warning('Could not locate triangle for x = [%f, %f]', x(1), x(2));
  g = [0; 0];  % safest fallback
  return;
end


cu_T = cu_glb((iTri-1)*M+1 : iTri*M);
Pts_Ti = meshP(:,meshT(iTri,:));
J = IncidenceMatrix(Pts_Ti); invJ = inv(J);
rs = invJ * (x - Pts_Ti(:,1)); r = rs(1); s = rs(2);

if r < -0.01 || s < -0.01 || r + s > 1.01
    warning('Bad reference coords: r = %.4f, s = %.4f', r, s);
end

V = jPoly_tri(r, s, varargin{1}, n-1, a, b, c);
fprintf("time = %.3e\n", V*cu_T);
tx = V*cu_T;
Vx = jPoly_tri(r, s, Hx, n-1, a+1, b, c+1); 
Vy = jPoly_tri(r, s, Hy, n-1, a, b+1, c+1); 
ur = Vx * (Dx * cu_T);
us = Vy * (Dy * cu_T);
g = invJ' * [ur; us];

end
