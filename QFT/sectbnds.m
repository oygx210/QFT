function [ubds,flag] = sectbnds(bds,flagin)
% SECTBNDS Compute an algebraic intersection of QFT bounds.
%          SECTBNDS(BNDS) computes the intersection of the QFT bounds
%          contained in BNDS.
%
%          See also SISOBNDS, GENBNDS, GRPBNDS, PLOTBNDS.

% Author: Yossi Chait
% Date: 4/18/94
% Revised: 2/16/96 10:51 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.

if nargin==1, flagin=1; end

myeps=1e-16;

ubds=[];

%%%%%% V5 code
% Reason: undefined variables
iz1 = []; iz2 = [];
iv1 = []; iv2 = [];

[r,c]=size(bds);
offset=(r-2)/2;
abvrng=1:offset; blwrng=offset+1:offset*2;
if flagin,
    dbmyeps=20*log10(myeps); db1myeps=20*log10(1/myeps);
else
    dbmyeps=myeps; db1myeps=1/myeps;
    A=bds(abvrng,:); B=bds(blwrng,:);
    bdb2=qclassfy(A',B');
    wbds=ones(1,4);
    bds=[bdb2;wbds;ones(1,4)];
end

wbds=bds(r-1,:); wbds=sort(wbds); wbds(find(diff(wbds)==0))=[];

flag=ones(length(wbds),1);         % legend: 1=ok, 2=myeps, 3=248, 4=302

bds2=bds;
for wvec=1:length(wbds),
    tbds=bds2(:,find(bds2(r-1,:)==wbds(wvec)));
    [rtbds,ctbds]=size(tbds);
    A=tbds(abvrng,1);B=tbds(blwrng,1);
    if ctbds==1,
        ubds(:,wvec)=[A ; B; wbds(wvec); 13];
    else
        for k=2:ctbds,
            A2=tbds(abvrng,k);  B2=tbds(blwrng,k);
            % find any real bounds
            i2=[]; i2=find(A==248 | A2==248);      %find 'em "bad" apples
            iinf=[]; iinf=find(A==302 | A2==302);  %find 'em "bad" apples
            % if all no LTI or non-connected go to next frequency
            if(length(i2)==offset),
                A=248*ones(offset,1);  B=-248*ones(offset,1);
                ubds(:,wvec)=[A; B; wbds(wvec); 13];
                flag(wvec,1)=3;  break;
            end
            if(length(iinf)==offset),
                A=302*ones(offset,1);  B=-302*ones(offset,1);
                ubds(:,wvec)=[A; B; wbds(wvec); 13];
                flag(wvec,1)=4;  break;
            end

            ia=zeros(offset,1);    i1=find(A~=dbmyeps & A~=248 & A~=302);
            if (i1), ia(i1)=21*ones(length(i1),1); end
            ib=zeros(offset,1);    i1=find(B~=db1myeps & B~=-248 & B~=-302);
            if (i1), ib(i1)=15*ones(length(i1),1); end
            ia2=zeros(offset,1);   i1=find(A2~=dbmyeps & A2~=248 & A2~=302);
            if (i1), ia2(i1)=21*ones(length(i1),1); end
            ib2=zeros(offset,1);   i1=find(B2~=db1myeps & B2~=-248 & B2~=-302);
            if (i1), ib2(i1)=15*ones(length(i1),1); end

            % initialize
            is=zeros(offset,1); is2=zeros(offset,1); ix=zeros(offset,1);
            iy=zeros(offset,1); ix2=zeros(offset,1); iy2=zeros(offset,1);
            iw=zeros(offset,1); iw2=zeros(offset,1);
            iwx=zeros(offset,1); iw2x=zeros(offset,1); iout=zeros(offset,1);
            iza=[]; izb=[]; iz1a=[]; iz1b=[]; iz2a=[]; iz2b=[]; iz3a=[]; iz3b=[];
            i3=[]; i5=[]; i7=[]; ie=[];

            if (any(ia) | any(ib)) | (any(ia2) | any(ib2)),      %if there's a real bound

                % 1st bound
                t=find((ia-ib)==21);      ix(t)=ones(length(t),1);        % 1
                t=find((ia-ib)==(-15));   ix2(t)=ones(length(t),1);       % 4
                t=find((ia-ib)==6);       is(t)=ones(length(t),1);        % 2 and/or 3
                t1=find(A(t)>=B(t));      iw(t(t1))=ones(length(t1),1);   % 3
                t1=find(A(t)<B(t));       iw2(t(t1))=ones(length(t1),1);  % 2
                % 2nd bound
                t=find((ia2-ib2)==21);    iy(t)=ones(length(t),1);        % 1
                t=find((ia2-ib2)==(-15)); iy2(t)=ones(length(t),1);       % 4
                t=find((ia2-ib2)==6);     is2(t)=ones(length(t),1);       % 2 and or 3
                t1=find(A2(t)>=B2(t));    iwx(t(t1))=ones(length(t1),1);  % 3
                t1=find(A2(t)<B2(t));     iw2x(t(t1))=ones(length(t1),1); % 2

                % test for empty intersection
                t=[find((ix+iy2)==2);find((ix+iw2x)==2);find((iwx+iy2)==2);];  % 1-4,1-2,2-4
                if (t),  iz1=t(find(A(t)>B2(t)));   end

                t=[find((ix2+iy)==2);find((ix2+iw2)==2);find((iw2x+ix2)==2)]; % 4-1,2-1,4-2
                if (t),  iz2=t(find(A2(t)>B(t)));   end

                t=find((iw+iw2x)==2);                                         % 2-3
                if (t),  iz3a=t(find((A(t)>B2(t)) | (B(t)<A2(t)))); end

                t=find((iw2+iwx)==2);                                         % 3-2
                if (t),  iz3b=t(find((A2(t)>B(t)) | (B2(t)<A(t)))); end

                % finally, intersect bounds
                iz=sort([iz1;iz2,iz3a,iz3b]);   % no LTI cases
                if(iz),
                    A(iz)=ones(length(iz),1)*248; B(iz)=-ones(length(iz),1)*248;
                end;

                % test for non-connected intersection in 2-3 or 3-2 (if not already no LTI case)
                ie=find(A==dbmyeps | B2==db1myeps | A2==dbmyeps | B==db1myeps);
                ie=[ie;iz]; if (ie), iout(ie)=5*ones(length(ie),1); end
                t=find((iw+iw2x+iout)==2);
                if (t), iv1=t(find((A(t)<=B2(t)) & (B(t)>=A2(t)))); end
                t=find((iw2+iwx+iout)==2);
                if (t), iv2=t(find((A2(t)<=B(t)) & (B2(t)>=A(t)))); end
                iv=[iv1;iv2];

                if(iv),
                    A(iv)=A(iv)*0+302;  B(iv)=B(iv)*0-302;
                    ubds(:,wvec)=[A; B; wbds(wvec); 13];
                    flag(wvec,1)=4;
                end

                % enforce original no LTI cases and non-connected cases
                if (i2),
                    A(i2)=ones(length(i2),1)*248; B(i2)=-ones(length(i2),1)*248;
                end
                if (iinf),
                    A(iinf)=ones(length(iinf),1)*302; B(iinf)=-ones(length(iinf),1)*302;
                end

                % take simple min and max on remaining cases

                i1=ones(offset,1);
                if (iz),  i1(iz)=zeros(length(iz),1);  end
                if (i2),  i1(i2)=zeros(length(i2),1);  end
                if length(i1==1),
                    i7=sort([find(ix+iwx==2);find(iy+iw==2)]);    %1-3 and 3-1
                    i3=sort([find(iw+iy2==2);find(ix2+iwx==2)]);  %3-4 and 4-3
                    if ([i7;i3]),
                        i1([i7;i3])=zeros(length([i7;i3]),1);
                    end
                    i5=find(i1==1);
                    if (i7), A(i7)=max(A(i7),A2(i7)); B(i7)=db1myeps*ones(length(i7),1);  end
                    if (i3), A(i3)=dbmyeps*ones(length(i3),1); B(i3)=min(B(i3),B2(i3));   end
                    if (i5), A(i5)=max(A(i5),A2(i5));        B(i5)=min(B(i5),B2(i5));   end

                end
            else
                A=dbmyeps*ones(offset,1);  B=db1myeps*ones(offset,1);
                if (i2),
                    A(i2)=ones(length(i2),1)*248; B(i2)=-ones(length(i2),1)*248;
                end
                if (iinf),
                    A(iinf)=ones(length(iinf),1)*302; B(iinf)=-ones(length(iinf),1)*302;
                end
            end
            ubds(:,wvec)=[A; B; wbds(wvec); 13];
        end
    end

    i2=[]; i2=find(A==248);
    iinf=[]; iinf=find(A==302);
    imyeps=[]; imyeps=find(A==dbmyeps & B==db1myeps);
    if (length(i2)==offset), flag(wvec,1)=3;
    elseif length(iinf), flag(wvec,1)=4;
    elseif (length(imyeps)==offset), flag(wvec,1)=2; end

end

if flagin,
    if any(flag==2),
        disp(['No bound found at w = ',sprintf('%4.2g ',wbds(flag==2))]);
    end
    if any(flag==3),
        disp(['No LTI solution possible at w = ',sprintf('%4.2g ',wbds(flag==3))]);
    end
    if any(flag==4),
        disp(['Computing non-connected bound sets at w = ',sprintf('%4.2g ',wbds(flag==4))]);
        disp('Do not intersect');
    end
end
