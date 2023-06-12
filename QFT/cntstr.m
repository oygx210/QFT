function [ControllerString, ListboxInfo, ListboxValue] = cntstr(f,cont,ElementLocation)
% CNTSTR Controller string for listbox. (Utility Function)

% Author: Craig Borghesani
% Date: 1/22/02 10:38PM
% Copyright (c) 2003, Terasoft, Inc.

if nargin == 2,
   ElementLocation = 1;
end

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
T=get(bthan(13),'userdata');
hint_bar = get(bthan(36),'userdata');
sortByType = findobj(f,'tag','sortByType');

ElementEnvironment = infmat(9,1);
delay = infmat(10,1);
qflag= 0;

% list of all element code numbers
% complex lead: 7
% notch: 6
% lead/lag: 5
% second-order zero: 4
% second-order pole: 3
% real zero: 2
% real pole: 1
% continuous integrator: 0.7
% discrete delay/predictor: 0.5
% discrete integrator: 0.6
% gain: 0

if strcmp(get(sortByType,'checked'),'on'),
   elem=[0 0.6 0.5 0.7 1 2 3 4 5 6 7];
end

ListboxValue = [];

ControllerString={};
ListboxInfo=[];
ListboxCounter = 1;
for k = elem,
   l=find(k==cont(:,4));
   if length(l),
      if strcmp(get(sortByType,'checked'),'on'),
         [sortvec,indexvec]=sort(cont(l,(1+any(k==[3:5])+2*(k==6))));
         l = l(indexvec); l = fliplr(l(:)');
      end
      for h=1:length(l),
         st = [];
         lh = l(h);
         v1n=cont(l(h),1); v1=num2str(v1n,4);
         v2n=cont(l(h),2); v2=num2str(v2n,4);
         v3n=cont(l(h),3); v3=num2str(v3n,4);
         val=l(h);

         if k==0,
            st = ['Gain: [',v1,']'];
%         elseif k==0.7,
%            val=0;
%            if cont(2,1)>0,
%               st='Integrator: ';
%               if all(qflag~=[2,7]),
%                  st=[st,'[',v1,']'];
%               end
%            elseif cont(2,1)<0,
%               st='Differentiator: ';
%               if all(qflag~=[2,7]),
%                  st=[st,'[',int2str(abs(v1n)),']'];
%               end
%            end
         elseif k==0.6,
            val=0;
            if any(cont(3,1:2)~=0),
               st='Delay/Pred: ';
               if all(qflag~=[2,7]),
                  st=[st,'[',int2str(cont(3,1)),'/',int2str(cont(3,2)),']'];
               end
            end
%            if cont(3,1)>0,
%               st='Delay: ';
%               if all(qflag~=[2,7]),
%                  st=[st,'[',int2str(abs(v1n)),']'];
%               end
%            elseif cont(3,1)<0,
%               st='Predictor: ';
%               if all(qflag~=[2,7]),
%                  st=[st,'[',int2str(abs(v1n)),']'];
%               end
%            end
         elseif any(k==[0.5,0.7]),
            val=0;
            if any(cont(2,1:2)~=0),
               st='Integ/Diff: ';
               if all(qflag~=[2,7]),
                  st=[st,'[',int2str(cont(2,1)),'/',int2str(cont(2,2)),']'];
               end
            end
         elseif k==1,
            if v1n<0, val=0; end
            st=['Real Pole: [',v1,']'];
            if T > 0 & all(qflag~=[4,6]),
               v1z=real(exp(-v1n*T));
               v1nz=num2str(v1z,8);
               st=[st,' [re=',v1nz,']'];
            end
         elseif k==2,
            st=['Real Zero: [',v1,'] '];
            if T > 0 & all(qflag~=[4,6]),
               v1z=real(exp(-v1n*T));
               v1nz=num2str(v1z,8);
               st=[st,' [re=',v1nz,']'];
            end
         elseif k==3,
            rt=roots([1,2*v1n*v2n,v2n^2]);
            if T > 0,
               rt=exp(rt*T);
            end
            if T > 0 & any(abs(rt)>1), val=0;
            elseif T == 0 & any(real(rt)>=0), val=0;
            end
            st=['Complex Pole: [',v1,', ',v2,']'];
            if all(qflag~=[4,6]),
               re=num2str(-real(rt(1)),4);
               im=num2str(imag(rt(1)),4);
               st=[st,' [re=',re,', im=',im,']'];
            end
         elseif k==4,
            rt=roots([1,2*v1n*v2n,v2n^2]);
            if T > 0,
               rt=exp(rt*T);
            end
            if T > 0 & any(abs(rt)>1), val=0;
            elseif T == 0 & any(real(rt)>=0), val=0;
            end
            st=['Complex Zero: [',v1,', ',v2,']'];
            if all(qflag~=[4,6]),
               re=num2str(-real(rt(1)),4);
               im=num2str(imag(rt(1)),4);
               st=[st,' [re=',re,', im=',im,']'];
            end
         elseif k==5,
            [jk,z,p]=ldlgcplx(v1n,v2n,[],T);
            st=['Lead/Lag: [',v1,', ',v2,']'];
            if all(qflag~=[4,6]),
               st=[st,' [z=',z,', p=',p,']'];
            end
         elseif k==6,
            st=['Notch: [',v1,', ',v2,', ',v3,']'];
            if all(qflag~=[4,6]),
               rt1=roots([1,2*v1n*v3n,v3n^2]);
               if T > 0,
                  rt1=exp(rt1*T);
               end
               rt2=roots([1,2*v2n*v3n,v3n^2]);
               if T > 0,
                  rt2=exp(rt2*T);
               end
               re1=num2str(-real(rt1(1)),4);
               im1=num2str(imag(rt1(1)),4);
               re2=num2str(-real(rt2(1)),4);
               im2=num2str(imag(rt2(1)),4);
               st=[st,' [re1=',re1,', im1=',im1,', re2=',re2,', im2=',im2,']'];
            end
         elseif k==7,
            st=['Complex Lead: [',v1,', ',v2,', ',v3,']'];
%            if all(qflag~=[4,6]),
%               rt1=roots([1,2*v1n*v3n,v3n^2]);
%               if T > 0,
%                  rt1=exp(rt1*T);
%               end
%               rt2=roots([1,2*v2n*v3n,v3n^2]);
%               if T > 0,
%                  rt2=exp(rt2*T);
%               end
%               re1=num2str(-real(rt1(1)),4);
%               im1=num2str(imag(rt1(1)),4);
%               re2=num2str(-real(rt2(1)),4);
%               im2=num2str(imag(rt2(1)),4);
%               st=[st,' [re1=',re1,', im1=',im1,', re2=',re2,', im2=',im2,']'];
%            end
         end % if k==0
         if ~isempty(st),
            ControllerString = [ControllerString, {st}];
            ListboxInfo = [ListboxInfo, lh];
            if l(h) == ElementLocation,
               ListboxValue = ListboxCounter;
            else
               ListboxCounter = ListboxCounter + 1;
            end
         end % if ~isempty(st)
      end % for h=1:length(l)
   end % if length(l)
end % for k = elem

if isempty(ListboxValue),
   ListboxValue = 1;
end
