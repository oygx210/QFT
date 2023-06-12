function mesgbnds(w,pbds,state,state2,ph_r,bnd)
% MESGBNDS Warning messages. (Utility Function)
%          MESGBNDS handles the different messages displayed to the user
%          concerning the different possible states of the computed bounds.

% Author: Craig Borghesani
% 5/25/94
% Copyright (c) 2003, Terasoft, Inc.


myeps = 1e-16;
dbmyeps = 20*log10(myeps); dbmy1eps = 20*log10(1/myeps);
final_state=state;

if length(state2), final_state=state2; end

% notify user of possible numerical inaccuracy
offset = length(ph_r);
ph180 = find(abs(ph_r)==pi);
wbds_inac = [];
size_bnd = size(bnd);
if length(ph180),
 for k=1:size_bnd(2),
  abv180 = bnd(ph180,k)==dbmyeps;
  blw180 = bnd(ph180+offset,k)==dbmy1eps;
  abv_all = ~all(bnd(1:offset,k)==dbmyeps);
  blw_all = ~all(bnd((offset+1):(2*offset),k)==dbmy1eps);
  if abv180 & blw180 & abv_all & blw_all,
   wbds_inac = [wbds_inac,bnd(offset*2+1,k)];
  end
 end
 if length(wbds_inac),
  wbds_inac=sort(wbds_inac); wbds_inac(find(diff(wbds_inac)==0))=[];
  disp(['Bound at w = ',sprintf('%4.4g ',wbds_inac),'suspect to numerical inaccuracies.']);
 end
end

% notify user of no bound found
if any(final_state==2),
 disp(['No bound found at w = ',sprintf('%4.4g ',w(pbds(final_state==2)))]);
end

% notify user of No LTI solution
if any(final_state==3),
 disp(['No LTI solution possible at w = ',sprintf('%4.4g ',w(pbds(final_state==3)))]);
end

% notify user of non-connected bound sets
if any(final_state==4),
 disp(['Non-connected bound sets at w = ',sprintf('%4.4g ',w(pbds(final_state==4)))]);
end
