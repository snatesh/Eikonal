function L = lowering_mat_tri(a,b,c,H,H1,mode)
n = size(H,1)-2;
a = a-1/2; b = b-1/2; c = c-1/2; 
N = (n+1)*(n+2)/2;
L = zeros(N);
% (a,b,c) -> (a-1,b,c)
if mode == 0
row = 1; col = 1;
for blk = 0:n
  kk = 0:blk;
  cPnk = H1(kk+1,blk+1)'.*(blk-kk+a) ./ ...
    (H(kk+1,blk+1)'.*(2*blk+a+b+c+2));
  cPnp1k = H1(kk+1,blk+2)'.*(blk-kk+1) ./ ...
    (H(kk+1,blk+1)'.*(2*blk+a+b+c+2));
  for j = 0:blk
    if row+j <= N && col+j <= N
    L(row+j,col+j) = cPnk(j+1);
    end
    if row+j <= N && col+j+blk+1 <= N
    L(row+j,col+j+blk+1) = cPnp1k(j+1);
    end
  end
  row = row+blk+1;
  col = col+blk+1;
end
  
% (a,b,c) -> (a,b-1,c)
elseif mode == 1
L(1,1) = b*H1(1,1)/((a+b+c+2)*H(1,1));
L(1,2) = -b*H1(1,2)/((b+c+1)*(a+b+c+2)*H(1,1));
L(1,3) = H1(2,2)/((b+c+1)*H(1,1));
row = 2; col = 2;
for blk = 1:n
  kk = 0:blk;
  cPnk = H1(kk+1,blk+1)'.*(kk+b).*(blk+kk+b+c+1) ./ ...
    (H(kk+1,blk+1)'.*(2*kk+b+c+1).*(2*blk+a+b+c+2));
  cPnkp1 = -H1(kk+2,blk+1)'.*(kk+1).*(blk-kk+a) ./ ...
    (H(kk+1,blk+1)'.*(2*kk+b+c+1).*(2*blk+a+b+c+2));
  cPnp1k = -H1(kk+1,blk+2)'.*(kk+b).*(blk-kk+1) ./ ...
    (H(kk+1,blk+1)'.*(2*kk+b+c+1).*(2*blk+a+b+c+2));
  cPnp1kp1 = H1(kk+2,blk+2)'.*(kk+1).*(blk+kk+a+b+c+2) ./ ...
    (H(kk+1,blk+1)'.*(2*kk+b+c+1).*(2*blk+a+b+c+2)); 
  cPnkp1(blk+1) = 0;
  for j = 0:blk
    if row+j <= N
      if col+j <= N
        L(row+j,col+j) = cPnk(j+1);
      end
      if col+j+1 <= N
        L(row+j,col+j+1) = cPnkp1(j+1);
      end
      if col+j+blk+1 <= N
        L(row+j,col+j+blk+1) = cPnp1k(j+1);
      end
      if col+j+blk+2 <= N
        L(row+j,col+j+blk+2) = cPnp1kp1(j+1);
      end
    end
  end
  row = row+blk+1;
  col = col+blk+1;
end
else
  error('invalid mode (0,1)');
end
% we write as transpose 
L = L';
end