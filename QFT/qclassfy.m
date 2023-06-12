function cbdb = qclassfy(A,B)
% QCLASSFY Classifies QFT bounds into four basic sets. (Utility Function)
%          QCLASSFY takes a set of above and below bounds and consolidates
%          them into a set of only 4 bounds.  This process greatly improves
%          the performance of SECTBNDS.

% Author: Yossi Chait
% 7/15/94
% Copyright (c) 2003, Terasoft, Inc.


[ra,ca]=size(A);

myeps=1e-16;
dbeps=myeps; db1eps=1/myeps;

% pre-allocate bound matrices
Amat=ones(ra,ca)*dbeps; Bmat=ones(ra,ca)*db1eps;
A2=Amat(1,:);
B1=Bmat(1,:);

if ra==1,  % only one bound passed
 cbdb=[A';B'];
else

 locAbad1=find(A==248);      %find 'em "bad" apples for no lti
 locBbad1=find(B==-248);
 locAbad2=find(A==302);      %find 'em "bad" apples for non-connected
 locBbad2=find(B==-302);

 ia=zeros(ra,ca); ib=ia;

% find any real bounds
 loc1a=find(A~=248 & A~=302 & A~=dbeps);
 ia(loc1a)=ia(loc1a)*0+21;

 loc2b=find(B~=-248 & B~=-302 & B~=db1eps);
 ib(loc2b)=ib(loc2b)*0+15;

 loc1=find((ia-ib)==21);                % 1
 loc4=find((ia-ib)==(-15));             % 4
 loc23=find((ia-ib)==6);                % 2 and/or 3
 At=Amat; Bt=Bmat; At(loc23)=A(loc23); Bt(loc23)=B(loc23);
 loc3=find(At>=Bt);
 At=Bmat; Bt=Amat; At(loc23)=A(loc23); Bt(loc23)=B(loc23);
 loc2=find(At<Bt);

 At=Amat;
 At(loc1)=A(loc1);
 A1=max(At);

 Bt=Bmat;
 Bt(loc4)=B(loc4);
 B2=min(Bt);

 At=Amat; Bt=Bmat;
 At(loc3)=A(loc3); Bt(loc3)=B(loc3);
 A3=max(At); B3=min(Bt);

 At=Amat; Bt=Bmat;
 At(loc2)=A(loc2); Bt(loc2)=B(loc2);
 A4=max(At); B4=min(Bt);

% reassign new above and below bound set
 Anew=[A1;A2;A3;A4];  Bnew=[B1;B2;B3;B4];

% assign the "bad" ones to all cases for No LTI
 bad1colA=sort(ceil(locAbad1/ra));
 if length(bad1colA), bad1colA(find(diff(bad1colA)==0))=[]; end
 bad1colB=sort(ceil(locBbad1/ra));
 if length(bad1colB), bad1colB(find(diff(bad1colB)==0))=[]; end
 Anew(:,bad1colA)=ones(4,length(bad1colA))*248;
 Bnew(:,bad1colB)=ones(4,length(bad1colB))*(-248);

% assign the "bad" ones to all cases for non-connected
 bad2colA=sort(ceil(locAbad2/ra));
 if length(bad2colA), bad2colA(find(diff(bad2colA)==0))=[]; end
 bad2colB=sort(ceil(locBbad2/ra));
 if length(bad2colB), bad2colB(find(diff(bad2colB)==0))=[]; end
 Anew(:,bad2colA)=ones(4,length(bad2colA))*302;
 Bnew(:,bad2colB)=ones(4,length(bad2colB))*(-302);

% cbdb=[Anew';Bnew';ones(1,4);ones(1,4)];
  cbdb=[Anew';Bnew'];

end
%cbdb=[cbdb,[C;1;1]];
