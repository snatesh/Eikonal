% Minimal Hungarian / Munkres implementation (square C only)
function assign = hungarian(C)
  C = C - min(C,[],2);                % row reduce
  C = C - min(C,[],1);                % column reduce
  n = size(C,1);
  mask = zeros(n); rowCover = false(n,1); colCover = false(1,n);

  % Step 1: star zeros greedily
  for i=1:n
    for j=1:n
      if C(i,j)==0 && ~rowCover(i) && ~colCover(j)
        mask(i,j)=1; rowCover(i)=true; colCover(j)=true;
      end
    end
  end
  rowCover(:)=false; colCover(:)=false;

  % Step 2: cover columns with starred zeros
  colCover(any(mask==1,1)) = true;

  while sum(colCover)<n
    % Step 3: prime uncovered zero; if no starred zero in its row -> augment
    [r,c] = findZero(C,rowCover,colCover);
    while isempty(r)
      % Step 5: adjust matrix with smallest uncovered value
      m = min(C(~rowCover, ~colCover), [], 'all');
      C(rowCover,:) = C(rowCover,:) + m;
      C(:,~colCover) = C(:,~colCover) - m;
      [r,c] = findZero(C,rowCover,colCover);
    end
    mask(r,c) = 2;                            % prime
    s = find(mask(r,:)==1,1);                 % starred in row?
    if isempty(s)
      % Step 4: augment path
      path = [r c];
      done = false;
      while ~done
        rr = find(mask(:,path(end,2))==1,1);
        if isempty(rr)
          done = true;
        else
          path = [path; rr path(end,2)];
          cc = find(mask(path(end,1),:)==2,1);
          path = [path; path(end,1) cc];
        end
      end
      % flip stars/primes along path
      for k=1:size(path,1)
        if mask(path(k,1), path(k,2))==1
          mask(path(k,1), path(k,2))=0;
        else
          mask(path(k,1), path(k,2))=1;
        end
      end
      mask(mask==2)=0;
      rowCover(:)=false; colCover(:)=false;
      colCover(any(mask==1,1)) = true;       % cover columns with stars
    else
      rowCover(r)=true; colCover(s)=false;   % cover row, uncover that column
    end
  end

  % build assignment from starred zeros
  [~,assign] = max(mask,[],2);               % one star per row

  % helpers
  function [r,c] = findZero(C,rowCover,colCover)
    r=[]; c=[];
    for ii=1:n
      if ~rowCover(ii)
        jj = find(C(ii,:)==0 & ~colCover, 1);
        if ~isempty(jj), r=ii; c=jj; return; end
      end
    end
  end
end
