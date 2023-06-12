function [coora,coorb] = wherebnd(bdb)
% WHEREBND Location of bound within bound vector. (Utility Function)
%          WHEREBND returns location vectors designating the location of the
%          above and below bounds within the bound vector returned by
%          SISOBNDS and GENBNDS.

% Author: Craig Borghesani
% 10/5/92
% Copyright (c) 2003, Terasoft, Inc.


[row,col]=size(bdb);
frst(1)=1; last(1)=(row-2)/2;
frst(2)=(row-2)/2+1; last(2) = row-2;
coora=[]; coorb=[];
dbrange = [-250,250];

for q=1:col
    for h=1:2
        z = find( bdb(frst(h):last(h),q) > dbrange(1) & ...
                  bdb(frst(h):last(h),q) < dbrange(2));
        if( length(z) )
            if length(z)==1
                loc(1,1)=z(1);
                loc(1,2)=z(1);
            elseif length(z)==last(1)
                loc(1,1)=z(1);
                loc(1,2)=z(last(1));
            else
                loc(1,1)=z(1);
                zbreak=find(diff(z)>1);

                if( ~length(zbreak) )
                    loc(1,2)=z(length(z));
                else
                    for k=1:length(zbreak)
                        loc(k,2)=z(zbreak(k));
                        loc(k+1,1)=z(zbreak(k)+1);
                    end
                    loc(k+1,2)=z(length(z));
                end
            end
        else
            loc(1,1)=0; loc(1,2)=0;
        end
        if( h==1 )
            loca = loc;
        else
            locb = loc;
        end
        if( h==2 )
            [ lar, lac ] = size(loca);
            [ lbr, lbc ] = size(locb);
            if( lar > lbr )
                locb = [ locb;
                         zeros(lar-lbr, 2) ];
            elseif( lbr > lar )
                loca = [ loca;
                         zeros(lbr-lar, 2) ];
            end
            [ lar , lac ] = size( loca );
            coora = [coora; q*ones(lar,1) loca]; loc=[];
            coorb = [coorb; q*ones(lar,1) locb]; loc=[];
        end
    end
end
