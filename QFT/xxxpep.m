function Z = xxxpep(X,Y,tp)
% XXXPEP
%       Z=xxxpep(X,Y)      calculates Z=diax(diax(X)*diax(Y)),
%       Z=xxxpep(X,A,'x#') calculates Z=diax(X)*A, and
%       Z=xxxpep(A,Y,'#x') calculates Z=A*diax(Y) in an efficient way.
%
%         xxxpep(X,Y,[]) and
%         xxxpep(X,Y,'xx') are equivalent with xxx(X,Y)

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


if nargin==2
  tp='xx';
elseif length(tp)==0
  tp='xx';
end
[n,nc]=size(Y);
if tp=='xx'
  Z=[X(:,1).*Y(:,1)+X(:,2).*Y(n:-1:1,2),...
     X(:,1).*Y(:,2)+X(:,2).*Y(n:-1:1,1)];
else
  if tp=='#x'
    Z=X;
    X=conj([Y(:,1) Y(n:-1:1,2)]);
    Y=Z';
    [n,nc]=size(Y);
  end
%  Z=(X(:,1)*ones(1,nc)).*Y+(X(:,2)*ones(1,nc)).*Y(n:-1:1,:);
  if n==2
    Z=[X(:,1) X(:,1)].*Y+[X(:,2) X(:,2)].*Y(n:-1:1,:);
  else
    Z=X(:,ones(1,nc)).*Y+X(:,2*ones(1,nc)).*Y(n:-1:1,:);
  end
  if tp=='#x'
    Z=Z';
  end
end
