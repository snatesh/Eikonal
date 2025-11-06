clear all; close all; clc;

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2; kap = abs(a+b+c);
m = 20;
[Jx,Jy,A1,A2,B1,B2,Hn] = jMatON_tri(m,a,b,c);
[V, D] = joint_diag_real([Jx Jy], 1e-12);
   
n = size(Jx,1);
D1 = D(:,1:n);
D2 = D(:,n+1:2*n);
offD1 = D1 - diag(diag(D1));
offD2 = D2 - diag(diag(D2));

col_res_x = sqrt(sum( offD1.^2 , 1 )); % 1-by-n
col_res_y = sqrt(sum( offD2.^2 , 1 ));
r = max([col_res_x; col_res_y], [], 1).';              % n-by-1

fprintf('r stats: min %.2e  median %.2e  max %.2e\n', min(r), median(r), max(r));

% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(m+1,a,b,c);
R = diag(D1); S = diag(D2);
% vandermonde under (a,b,c), (a+1,b,c+1), (a,b+1,c+1)
V_abc = jPoly_tri(R,S,H_abc,m-1,a,b,c);

%%

% jacobi matrices
n = size(Jx,1);
[V,D] = joint_diag_real([Jx Jy]);
D1 = D(:,1:n);
D2 = D(:,n+1:2*n);
offD1 = D1 - diag(diag(D1));
offD2 = D2 - diag(diag(D2));

col_res_x = sqrt(sum( offD1.^2 , 1 )); % 1-by-n
col_res_y = sqrt(sum( offD2.^2 , 1 ));
r = max([col_res_x; col_res_y], [], 1).';              % n-by-1

fprintf('r stats: min %.2e  median %.2e  max %.2e\n', min(r), median(r), max(r));
F2 = norm(Jx,'fro')^2 + norm(Jy,'fro')^2; 
delta_off_norm = sqrt( sum(offD1(:).^2) + sum(offD2(:).^2) ) / sqrt(F2);
fprintf('normalized delta_off = %.3e\n', delta_off_norm);

% simple greedy alignment of (diag D1, diag D2) to make the scatter “tight”
x = diag(D1); y = diag(D2);
[~,perm] = sortrows([x y],[1 2]);      % crude; Hungarian gives best
x = x(perm); y = y(perm);
xyeigs = eig(Jx + 1j*Jy);

plot(x,y,'.'); axis equal; xlim([0 1]); ylim([0 1]); hold on;
plot(real(xyeigs),imag(xyeigs),'o');

%%
% Inputs: Jx, Jy (n×n).  Run both methods:
[Ujevd, Dcat] = joint_diag_real([Jx Jy]);    % Cardoso–Souloumiac
n = size(Jx,1);
Dx = Dcat(:,1:n); Dy = Dcat(:,n+1:2*n);
Xj = diag(Dx); Yj = diag(Dy);                 % JEVD seeds

[Vc, z] = eig(Jx + 1i*Jy);                    % complexification
Xc = real(diag(Vc' * Jx * Vc));               % Ritz coords in Vc basis
Yc = real(diag(Vc' * Jy * Vc));

% Build cost matrix of pairwise distances and solve assignment
Xj = Xj(:); Yj = Yj(:); Xc = Xc(:); Yc = Yc(:);
n  = numel(Xj);
C  = hypot(Xj - Xc.', Yj - Yc.');   % n×n

% --- solve assignment exactly
assign = hungarian(C);              % returns a 1×n vector: column chosen for each row
idxJ = (1:n).'; idxC = assign(:);

% --- reorder and plot
Xj = Xj(idxJ); Yj = Yj(idxJ);
Xc = Xc(idxC); Yc = Yc(idxC);
plot(Xc,Yc,'k.', Xj,Yj,'bo'); axis equal; xlim([0 1]); ylim([0 1]);

%%
col_res_x = sqrt(sum(offD1.^2,1));
col_res_y = sqrt(sum(offD2.^2,1));
r = max(col_res_x, col_res_y).';   % n×1

fprintf('r stats: min %.2e  median %.2e  max %.2e\n', min(r), median(r), max(r));

good = r <= 1e-3;                   % pick your safety margin δ
plot(Xc(good),Yc(good),'k.', Xj(good),Yj(good),'bo'); axis equal