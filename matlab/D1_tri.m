function D = D1_tri(a,b,c,H,H1,mode)
n = size(H,1)-2;
a = a-1/2; b = b-1/2; c = c-1/2; 
N = (n+1)*(n+2)/2;
D = zeros(N);
% Dx
if mode == 0
D(2,1) = (1+a+b+c+2)*H1(1,1)/H(1,2);
D(3,1) = (1+b)*(1+1+b+c+1)*H1(1,1)/((2+b+c+1)*H(2,2));
row = 4; col = 2;
for blk = 2:n
  kk = 0:blk;
  cPnm1k = H1(kk+1,blk)'.*(blk+kk+a+b+c+2).*(kk+b+c+1) ./ ...
    (H(kk+1,blk+1)'.*(2*kk+b+c+1));
  cPnm1km1 = [0,H1(kk+1,blk)'.*(kk+1+b).*(blk+kk+1+b+c+1) ./ ...
    (H(kk+2,blk+1)'.*(2*kk+2+b+c+1))]; 
  cPnm1k(blk+1) = 0;
  for j = 0:blk
    D(row+j,col+j) = cPnm1k(j+1);
    D(row+j,col+j-1) = cPnm1km1(j+1);
  end
  row = row+blk+1;
  col = col+blk;
end
  
% Dy
elseif mode == 1
D(3,1) = (1+b+c+1)*H1(1,1)/H(2,2);
row = 4; col = 2;
for blk = 2:n
  kk = 0:blk;
  cPnm1km1 = [0,H1(kk+1,blk)'.*(kk+1+b+c+1)./H(kk+2,blk+1)']; 
  for j = 0:blk
    D(row+j,col+j-1) = cPnm1km1(j+1);
  end
  row = row+blk+1;
  col = col+blk;
end
else
  error('invalid mode (0,1)');
end
% we write as transpose 
D = D';
end