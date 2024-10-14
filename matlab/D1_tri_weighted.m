function D = D1_tri_weighted(a,b,c,H,H1,mode)
n = size(H,1)-2;
a = a-1/2; b = b-1/2; c = c-1/2; 
N = (n+1)*(n+2)/2;
D = zeros(N);
disp(H)
disp(H1)
% Dx
if mode == 0
row = 1; col = 2;
for blk = 0:n
  kk = 0:blk;
  cPnp1k = H1(kk+1,blk+2)'.*(kk+c).*(blk-kk+1) ./ ...
    (-H(kk+1,blk+1)'.*(2*kk+b+c+1));
  cPnp1kp1 = H1(kk+2,blk+2)'.*(kk+1).*(blk-kk+a) ./ ...
    (-H(kk+1,blk+1)'.*(2*kk+b+c+1));
    
  for j = 0:blk
    if row+j <= N && col+j <= N
      D(row+j,col+j) = cPnp1k(j+1);
    end
    if row+j <= N && col+j+1 <= N
      D(row+j,col+j+1) = cPnp1kp1(j+1);
    end
  end
  row = row+blk+1;
  col = col+blk+2;
end
  
% Dy
elseif mode == 1
row = 1; col = 3;
for blk = 0:n
  kk = 0:blk;
  cPnp1kp1 =  -H1(kk+2,blk+2)'.*(kk+1)./H(kk+1,blk+1)';
  for j = 0:blk
    if row+j <= N && col+j <= N
      D(row+j,col+j) = cPnp1kp1(j+1);
    end
  end
  row = row+blk+1;
  col = col+blk+2;
end
else
  error('invalid mode (0,1)');
end
% we write as transpose 
D = D';
end