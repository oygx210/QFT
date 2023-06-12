function X = qlyaps(tp,A1,A2,E)
% QLYAPS Solve the general form of the Lyapunov matrix equation. (Utility)
%        X = QLYAPS(TP,A1,A2,E) solves the general form of the Lyapunov
%        matrix equation A1*X + X*A2 + E =0 in which A1 and A2 are special
%        matrices (triangular or diagonal)
%        TP : two character string defining the type of A1 and A2 matrices
%        'l': lower triangular
%        'u': upper triangular
%        'f': full
%        'd': diagonal
%        (fastest calculation for: TP='dd', then TP='du','uu','lu','fu')
%        A1 : square matrix of type tp(1), column vector for tp(1)=='d'
%        A2 : square matrix of type tp(2), column vector for tp(2)=='d'
%        E : general matrix of appropriate dimension
%        Note that          >> X=qlyaps(tp,A1,A2,E)
%        is equivalent with >> X=T1*qlyaps(tp,T1\A1*T1,T2\A2*T2,T1\E/T2)*T2

% Author: Pepijn Wortelboer
% 4/1/94
% Copyright (c) 2003, Terasoft, Inc.


if tp(2)=='l'
  [nr2,nc2]=size(A2);
  [nrE,ncE]=size(E);
  tp(2)='u';
  X=qlyaps(tp,A1,A2(nr2:-1:1,nc2:-1:1),E(:,ncE:-1:1));
  X=X(:,ncE:-1:1);
  return
end
if tp(1)=='d'
  A1=A1(:); nr1=length(A1);
else
  [nr1,nc1] = size(A1);
  if nr1~=nc1
    help qlyaps; error('matrix A1 not square');
  end
end
if tp(2)=='d'
  A2=A2(:).'; nr2=length(A2);
else
  [nr2,nc2] = size(A2);
  if nr2~=nc2
    help qlyaps; error('matrix A2 not square');
  end
end
[nrE,ncE] = size(E);
if nr1~=nrE | nr2~=ncE
  help qlyaps; error('matrix E not of compatible dimension');
end
if tp(2)=='f'
  [U2,A2]=schur(A2);
  if any(any(imag(A2)~=0))
    % already complex form
  else
    [U2,A2]=rsf2csf(U2,A2);
  end
  E=E*U2;
  tp(2)='U';
end
X(nrE,ncE) = 0;
if tp=='ld' | tp=='ud'
  A2=diag(A2); tp(2)='u';
end
if tp(1)=='l'
  A1=A1(nrE:-1:1,nrE:-1:1);
  E=E(nrE:-1:1,:);
end
if tp(2)=='l'
  A2=A2(ncE:-1:1,ncE:-1:1);
  E=E(:,ncE:-1:1);
end
if tp(1)=='d'
  if tp(2)=='d'
    X=-E./(A1(:,ones(1,ncE))+A2(ones(1,nrE),:));
  else
% solve for first column of transformed solution
    X(nrE,ncE) = 0;
    X(:,1) = -E(:,1)./(A1+ones(nrE,1)*A2(1,1));
% solve for remaining columns of transformed solution
    for k=2:ncE
      km1 = 1:(k-1);
      X(:,k) = -(E(:,k)+X(:,km1)*A2(km1,k))./...
                (A1+ones(nrE,1)*A2(k,k));
    end
  end
else
  % Solution based on lower triangularity of A2
  % solve for first column of transformed solution
  I1 = eye(nrE);
  X(:,1) = -(A1+I1*A2(1,1))\E(:,1);
  % solve for remaining columns of transformed solution
  for k=2:ncE
    km1 = 1:(k-1);
    X(:,k) = -(A1+I1*A2(k,k))\(E(:,k)+X(:,km1)*A2(km1,k));
  end
end
if tp(1)=='l'
  X=X(nrE:-1:1,:);
end
if tp(2)=='l'
  X=X(:,ncE:-1:1);
end
if tp(2)=='U'
  X=X*U2';
end
