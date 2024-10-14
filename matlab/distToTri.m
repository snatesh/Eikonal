function dist = distToTri(X,Y)
  
Np = length(X);
dist = zeros(Np,1);
for j = 1:Np
  x = [X(j),Y(j)];
  z = [(X(j)-Y(j)+1)/2, (Y(j)-X(j)+1)/2];
  dtoh = norm(x-z);
  dist(j) = min([X(j),Y(j),dtoh]);
end

end
