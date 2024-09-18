function Pfunc = Pk(k,a,b)

q = k.^2 .* (k + a + b) .* (k.^2 + a + b - 2);
r = (k.^2 + a + b - 1) .* (a.^2 - b.^2);
s = (k.^2 + a + b - 2) .* (k.^2 + a + b - 1) .* (k.^2 + a + b);
t = 2.0 * (k + a - 1) .* (k + b - 1) .* (k.^2 + a + b);


syms x;
if k == 0
  Pfunc = 1;
elseif k == 1
  Pfunc = 0.5 * ((a + b + 2) .* (x - 1) + 2 * (a + 1)); 
else
  Pfunc = ((r + s .* x) .* Pk(k-1,a,b) - t .* Pk(k-2,a,b)) ./ q; 
end


end
