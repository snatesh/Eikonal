function [V_abc,dxP,dyP,Vl,Vb,Vh,...
          Vl_flip,Vb_flip,Vh_flip,...
          Rl,Sl,Rb,Sb,Rh,Sh,...
          Rl_flip,Sl_flip,Rb_flip,Sb_flip,...
          Rh_flip,Sh_flip,wleg,...
          intVlVl,intVbVb,intVhVh,...
          intVlVl_flip,intVbVb_flip,...
          intVhVh_flip,varargout] = preAssemble_poisson1(n,R,S)
a = 1/2; b = a; c = a;
M = n*(n+1)/2;

% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);

% vandermonde under (a,b,c), (a+1,b,c+1), (a,b+1,c+1)
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
V_a1bc1 = jPoly_tri(R,S,H_a1bc1,n-1,a+1,b,c+1);
V_ab1c1 = jPoly_tri(R,S,H_ab1c1,n-1,a,b+1,c+1);

% derivative amatrices
Dx_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);


Nrs = length(R);

% interior basis derivatives at abscissa
dxP = V_a1bc1 * Dx_a1bc1;
dyP = V_ab1c1 * Dy_ab1c1;
%for i = 1:M
%    dxP(:,i) = V_a1bc1 * Dx_a1bc1(:,i);
%    dyP(:,i) = V_ab1c1 * Dy_ab1c1(:,i);
%end

nleg = 50;
[xleg,wleg,~] = gjQuad(nleg,0,0);
wleg = wleg'/2; xleg = (xleg+1)/2;
Rl = 0*xleg; Sl = xleg;
Rb = xleg; Sb = 0*xleg;
Rh = xleg; Sh = 1-Rh;

Vl = jPoly_tri(Rl,Sl,H_abc,n-1,a,b,c);
Vb = jPoly_tri(Rb,Sb,H_abc,n-1,a,b,c);
Vh = jPoly_tri(Rh,Sh,H_abc,n-1,a,b,c);


intVlVl = zeros(M,M);
intVbVb = zeros(M,M);
intVhVh = zeros(M,M);

for i = 1:M
    for j = 1:M
        VlVlij = reshape(Vl(:,i).*Vl(:,j),nleg,1);
        VbVbij = reshape(Vb(:,i).*Vb(:,j),nleg,1);
        VhVhij = reshape(Vh(:,i).*Vh(:,j),nleg,1);
        intVlVl(i,j) = wleg'*VlVlij;
        intVbVb(i,j) = wleg'*VbVbij;
        intVhVh(i,j) = wleg'*VhVhij;
    end
end

% now compute flipped versions of this for shared edges
% where the parametric shared edges runs [0,1] in one triangle,
% but [1,0] in the other

% Flip quadrature points (reverse edge orientation)
xleg_flip = 1 - xleg;

% Coordinates for flipped edges
Rl_flip = 0*xleg_flip; Sl_flip = xleg_flip;
Rb_flip = xleg_flip;  Sb_flip = 0*xleg_flip;
Rh_flip = xleg_flip;  Sh_flip = 1 - Rh_flip;

% Evaluate basis on flipped edge coordinates
Vl_flip = jPoly_tri(Rl_flip,Sl_flip,H_abc,n-1,a,b,c);
Vb_flip = jPoly_tri(Rb_flip,Sb_flip,H_abc,n-1,a,b,c);
Vh_flip = jPoly_tri(Rh_flip,Sh_flip,H_abc,n-1,a,b,c);

% Initialize flipped integrals
intVlVl_flip = zeros(M,M);
intVbVb_flip = zeros(M,M);
intVhVh_flip = zeros(M,M);

% Compute flipped inner products
for i = 1:M
    for j = 1:M
        VlVlij_flip = reshape(Vl_flip(:,i) .* Vl_flip(:,j), nleg, 1);
        VbVbij_flip = reshape(Vb_flip(:,i) .* Vb_flip(:,j), nleg, 1);
        VhVhij_flip = reshape(Vh_flip(:,i) .* Vh_flip(:,j), nleg, 1);
        intVlVl_flip(i,j) = wleg' * VlVlij_flip;
        intVbVb_flip(i,j) = wleg' * VbVbij_flip;
        intVhVh_flip(i,j) = wleg' * VhVhij_flip;
    end
end

% if using neumann conditions
if nargout == 40

Vl_a1bc1 = jPoly_tri(Rl,Sl,H_a1bc1,n-1,a+1,b,c+1);
Vl_ab1c1 = jPoly_tri(Rl,Sl,H_ab1c1,n-1,a,b+1,c+1);
Vl_a1bc1_flip = jPoly_tri(Rl_flip,Sl_flip,H_a1bc1,n-1,a+1,b,c+1);
Vl_ab1c1_flip = jPoly_tri(Rl_flip,Sl_flip,H_ab1c1,n-1,a,b+1,c+1);
Vb_a1bc1 = jPoly_tri(Rb,Sb,H_a1bc1,n-1,a+1,b,c+1);
Vb_ab1c1 = jPoly_tri(Rb,Sb,H_ab1c1,n-1,a,b+1,c+1);
Vb_a1bc1_flip = jPoly_tri(Rb_flip,Sb_flip,H_a1bc1,n-1,a+1,b,c+1);
Vb_ab1c1_flip = jPoly_tri(Rb_flip,Sb_flip,H_ab1c1,n-1,a,b+1,c+1);
Vh_a1bc1 = jPoly_tri(Rh,Sh,H_a1bc1,n-1,a+1,b,c+1);
Vh_ab1c1 = jPoly_tri(Rh,Sh,H_ab1c1,n-1,a,b+1,c+1);
Vh_a1bc1_flip = jPoly_tri(Rh_flip,Sh_flip,H_a1bc1,n-1,a+1,b,c+1);
Vh_ab1c1_flip = jPoly_tri(Rh_flip,Sh_flip,H_ab1c1,n-1,a,b+1,c+1);


dxPl = Vl_a1bc1 * Dx_a1bc1; 
dyPl = Vl_ab1c1 * Dy_ab1c1; 
dxPb = Vb_a1bc1 * Dx_a1bc1; 
dyPb = Vb_ab1c1 * Dy_ab1c1; 
dxPh = Vh_a1bc1 * Dx_a1bc1; 
dyPh = Vh_ab1c1 * Dy_ab1c1; 

dxPl_flip = Vl_a1bc1_flip * Dx_a1bc1; 
dyPl_flip = Vl_ab1c1_flip * Dy_ab1c1; 
dxPb_flip = Vb_a1bc1_flip * Dx_a1bc1; 
dyPb_flip = Vb_ab1c1_flip * Dy_ab1c1; 
dxPh_flip = Vh_a1bc1_flip * Dx_a1bc1; 
dyPh_flip = Vh_ab1c1_flip * Dy_ab1c1; 

varargout{1} = dxPl;
varargout{2} = dyPl;
varargout{3} = dxPb;
varargout{4} = dyPb;
varargout{5} = dxPh;
varargout{6} = dyPh;
varargout{7} = dxPl_flip;
varargout{8} = dyPl_flip;
varargout{9} = dxPb_flip;
varargout{10} = dyPb_flip;
varargout{11} = dxPh_flip;
varargout{12} = dyPh_flip;

end

end
