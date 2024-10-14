function K = promotion_mat_tri(a,b,c,H,H1,mode)
n = size(H,1)-2;
a = a-1/2; b = b-1/2; c = c-1/2; 
N = (n+1)*(n+2)/2;
K = zeros(N); 
% a+1
if mode == 0
  K(1,1) = 1; 
  row = 2; col = 1;
  for blk = 1:n  
    kk = 0:blk;
    cPnk = H1(kk+1,blk+1)'.*(blk+kk+a+b+c+2)./(H(kk+1,blk+1)'.*(2*blk+a+b+c+2));
    cPnm1k = H1(kk+1,blk)'.*(blk+kk+b+c+1)./(H(kk+1,blk+1)'.*(2*blk+a+b+c+2));
    cPnm1k(blk+1) = 0; 
    for j = 0:blk
      K(row+j,col+j) = cPnm1k(j+1);
      K(row+j,col+j+blk) = cPnk(j+1);
    end
    row = row+blk+1;
    col = col+blk;
  end
% b+1
elseif mode == 1
  K(1,1) = 1; 
  K(2,1) = -(1+a)*(b+c+1)*H1(1,1)/((2+a+b+c+2)*(b+c+1)*H(1,2)); 
  K(2,2) = (1+a+b+c+2)*(b+c+1)*H1(1,2)/((2+a+b+c+2)*(b+c+1)*H(1,2));
  K(3,1) = (1+c)*(2+b+c+1)*H1(1,1)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
  K(3,2) = -(1+c)*H1(1,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
  K(3,3) = (2+a+b+c+2)*(1+b+c+1)*H1(2,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
  row = 4; col = 2;
  for blk = 2:n  
    kk = 0:blk;
    cPnk = H1(kk+1,blk+1)'.*(blk+kk+a+b+c+2).*(kk+b+c+1) ./ ...
      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
    
    cPnm1k = -H1(kk+1,blk)'.*(blk-kk+a).*(kk+b+c+1) ./ ... 
      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
    
    cPnkm1 = [0,-H1(kk+1,blk+1)'.*(kk+1+c).*(blk-kk-1+1) ./ ...
      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
    
    cPnm1km1 = [0,H1(kk+1,blk)'.*(kk+1+c).*(blk+kk+1+b+c+1) ./ ... 
      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
    
    cPnm1k(blk+1) = 0; 
    for j = 0:blk
      K(row+j,col+j) = cPnm1k(j+1);
      K(row+j,col+j+blk) = cPnk(j+1);
      K(row+j,col+j-1) = cPnm1km1(j+1);
      K(row+j,col+j+blk-1) = cPnkm1(j+1);
    end
    row = row+blk+1;
    col = col+blk;
  end 
% c+1
elseif mode == 2
  K(1,1) = 1; 
  K(2,1) = -(1+a)*(b+c+1)*H1(1,1)/((2+a+b+c+2)*(b+c+1)*H(1,2)); 
  K(2,2) = (1+a+b+c+2)*(b+c+1)*H1(1,2)/((2+a+b+c+2)*(b+c+1)*H(1,2));
  K(3,1) = -(1+b)*(2+b+c+1)*H1(1,1)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
  K(3,2) = (1+b)*H1(1,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
  K(3,3) = (2+a+b+c+2)*(1+b+c+1)*H1(2,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
  row = 4; col = 2;
  for blk = 2:n  
    kk = 0:blk;
    cPnk = H1(kk+1,blk+1)'.*(blk+kk+a+b+c+2).*(kk+b+c+1) ./ ...
      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
    
    cPnm1k = -H1(kk+1,blk)'.*(blk-kk+a).*(kk+b+c+1) ./ ... 
      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
    
    cPnkm1 = [0,H1(kk+1,blk+1)'.*(kk+1+b).*(blk-kk-1+1) ./ ...
      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
    
    cPnm1km1 = [0,-H1(kk+1,blk)'.*(kk+1+b).*(blk+kk+1+b+c+1) ./ ... 
      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
    
    cPnm1k(blk+1) = 0; 
    for j = 0:blk
      K(row+j,col+j) = cPnm1k(j+1);
      K(row+j,col+j+blk) = cPnk(j+1);
      K(row+j,col+j-1) = cPnm1km1(j+1);
      K(row+j,col+j+blk-1) = cPnkm1(j+1);
    end
    row = row+blk+1;
    col = col+blk;
  end 
else
  error('invalid mode (0,1,2)');
end
% we write as transpose 
% so P_{a+1,b,c}^T K = P_{a,b,c}
K = K';
end  
%%
    %nrows = blk+1; ncols = blk+id;
    %K_blk = cell(n,1);
    %K_blk = zeros(nrows,ncols);
    %K_blk(j,j) = cPnm1k(j);
    %K_blk(j,j+blk) = cPnk(j);
    %K_blk(j,j) = 0;
    %K(row:row+nrows-1,col:col+ncols-1) = K_blk
    
    


% inds   = [1,1]; 
% for nn = 1:n-1
%   Jn1(inds(1):inds(2),inds(1):inds(2)) = B1{nn}; 
%   Jn2(inds(1):inds(2),inds(1):inds(2)) = B2{nn};
%   inds1 = inds + [nn,nn+1];
%   Jn1(inds(1):inds(2),inds1(1):inds1(2)) = A1{nn} ;
%   Jn1(inds1(1):inds1(2),inds(1):inds(2)) = A1{nn}';
%   Jn2(inds(1):inds(2),inds1(1):inds1(2)) = A2{nn} ;
%   Jn2(inds1(1):inds1(2),inds(1):inds(2)) = A2{nn}';
%   inds = inds + [nn,nn+1];
% end