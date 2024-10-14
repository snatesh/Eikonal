function scale = structure_factors(N,a,b)
% Scaling:
NN = (0:(N-1))';
scale = 2^(a+b+1)*gamma(NN+a+1).*gamma(NN+b+1) ./ ...
     ((2*NN+a+b+1).*gamma(NN+a+b+1).*factorial(NN));
scale(1) = 2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
end