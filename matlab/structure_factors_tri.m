function H = structure_factors_tri(N,a,b,c)
kap = abs(a+b+c);
wabc = gamma(kap+3/2)/(gamma(a+1/2)*gamma(b+1/2)*gamma(c+1/2));
H = zeros(N);
for nn = 0:(N-1)
  kk = (0:nn)';
  H(kk+1,nn+1) = sqrt(wabc./((2*nn+kap+1/2).*(2*kk+b+c)) .* ...
    gamma(nn+kk+b+c+1).*gamma(nn-kk+a+1/2).*gamma(kk+b+1/2) .* ...
    gamma(kk+c+1/2)./(factorial(nn-kk).*factorial(kk) .* ...
    gamma(nn+kk+kap+1/2).*gamma(kk+b+c)));
end
end