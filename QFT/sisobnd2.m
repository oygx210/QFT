function bdb = sisobnd2(w,wbd,W,uP,vP,R,nom,uC,vC,ctype,ph_r,info);
%     Compute QFT bounds for the following closed-loop configuration
%
%                          |     1     |
%                          |-----------| <= WS
%                          |  1 + PGH  |
%
%     SISOBNDS(2,W,WBD,WS,P,R,NOM,C,LOC,PHS) computes bounds at frequencies
%     designated by WBD.  WS is the performance specification, P is
%     the frequency response data of the plant (complex), R is the
%     disk radius for non-parametric uncertainty, NOM designates the
%     nominal plant and controller, C.  LOC specifies location of unknown
%     controller in the loop: 1 for G, 2 for H.  PHS specifies at which
%     phases (degrees) to compute bounds.
%
%     SISOBNDS(2,W,[],WS,P,[],NOM,[],[],PHS) computes bounds using default
%     values for WBD (all frequencies in W), R (0), C (1), CTYPE (1).

% Author: Craig Borghesani
% Date: 9/6/93
% Revised: 2/16/96 10:12 AM V5 changes made
% Copyright (c) 2003, Terasoft, Inc.


m=length(ph_r); nbd=length(w);
pbds=qsubset(wbd,w);
final_state=ones(1,nbd);
state2=[];
myeps=1e-16;

% preallocating matrix for more efficient programming
bmag=[ones(m,length(pbds))*myeps;ones(m,length(pbds))*1/myeps];

[rmp,cmp]=size(uP); [rmc,cmc]=size(uC);
str='Non-existent nominal case index for nominal';
if nom(1)>rmp, error([str,' plant case']);
elseif nom(2)>rmc, error([str,' controller case']); end

% declaring sizes of replicating matricies
lo2=max(rmp,rmc);

if repltest,
    u=ones(lo2,1); v=ones(1,m);
else
    u=ones(1,lo2); v=ones(1,m);
end

if rmp~=rmc, val=min(rmp,rmc); else val=1; end

j=1; pct=1; cct=1;
while j<=length(pbds),

    if cmp>1, pct=pbds(j); end
    if cmc>1, cct=pbds(j); end

    % offset phase by phases of nominal plants
    phi=ph_r-(vP(nom(1),pct)+vC(nom(2),cct));

    cnt=1;
    if rmp>rmc, rp=1:rmp; rc=cnt;
    elseif rmc>rmp, rp=cnt; rc=1:rmc;
    else rp=1:rmp; rc=1:rmc; end

    while cnt<=val,

        %%%%%% V4.2 code
        %  rad=R(rp,pct); mP=uP(rp,pct); mC=uC(rc,cct); pP=vP(rp,pct); pC=vC(rc,cct);
        %  Ws=W(rp,pct);

        %%%%%% V5 code
        % Reason: V5 does "the right thing" with replicating matrices.  in V4.2,
        % if a vector was 10 elements long and you indexed it with a boolean
        % matrix of the same length, the identical matrix was returned.  in V5,
        % a matrix of the first element is returned.  this is consistent
        % behavior and the if-statement below makes sure that the
        % vectors to be replicated are of length 1; which they will be if
        % length(rp) == 1

        rad=R(rp,pct); mP=uP(rp,pct); pP=vP(rp,pct); Ws=W(rp,pct);
        if length(rp) == 1,
            rad=rad(u); mP=mP(u); pP=pP(u); Ws=Ws(u);
        end

        mC=uC(rc,cct); pC=vC(rc,cct);
        if length(rc) == 1,
            mC=mC(u); pC=pC(u);
        end

        psi=pP+pC; psi=psi(:,v)+phi(u,:);
        A=mP.^2 .*mC.^2 .*(1-rad.^2); A=A(:,v);
        B=2*mP(:,v).*mC(:,v).*(cos(psi)-rad(:,v)./Ws(:,v));
        C=1- 1 ./Ws(:,v).^2;

        [g1,g2]=quadrtic(A,B,C); size_g1 = size(g1);
        cbdb=[[[g1';g2'],bmag(:,j)];ones(2,size_g1(1)+1)];
        [abvblw,state(j)]=sectbnds(cbdb,0);
        bmag(:,j)=abvblw(1:2*m);
        if rmp>rmc, rc=rc+1; cnt=rc;
        elseif rmc>rmp, rp=rp+1; cnt=rp; else cnt=2; end
    end

    z=find(bmag(:,j)~=myeps & bmag(:,j)~=1/myeps & ...
        bmag(:,j)~=-248 & bmag(:,j)~=248 & ...
        bmag(:,j)~=-302 & bmag(:,j)~=302);
    if length(z)
        bmag(z,j)=bmag(z,j)*uP(nom(1),pct)*uC(nom(2),cct);
    end

    frac=j/length(pbds);
    set(info(2),'xdata',[0,frac,frac,0]);
    set(info(3),'string',[int2str(floor(100*frac)),'%']);
    drawnow;

    j=j+1;
end
close(info(1));

bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302) = ...
    20*log10(bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302));
bnd=[bmag; w(pbds); 2*ones(1,length(pbds))];

if any(R~=0),
    bndt=qrobust(w,w(pbds),1 ./R,uP,vP,0,nom,uC,vC,1,ph_r);
    [bnd,state2]=sectbnds([bnd,bndt]);
    if length(bnd), bnd(2*m+2,:)=2*ones(1,length(pbds)); end
end

mesgbnds(w,pbds,state,state2,ph_r,bnd);

[jk,t]=sort(w(pbds));
if nargout==0,
    plotbnds(bnd(:,t),[],180/pi*ph_r);
    title('SISOBND2 Bounds'),xlabel('X: Phase (degrees)  Y: Magnitude (dB)')
else
    bdb=bnd(:,t);
end
