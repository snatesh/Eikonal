function phij = phij_etd(j,Letd,v,tol,Ntay_max)
phijz = 0;
v0 = v; phij = v0; 
for k = 0:Ntay_max
  vk = Letd*v0;
  term = vk/factorial(k+j);
  phij = phij + term;
  err = norm(term);
  if err < tol
    %fprintf('phij converged in %d iterations with relnorm %.3e\n',k,err)
    break;
  end
end
end
