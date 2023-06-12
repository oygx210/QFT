function qprint(UIOperation)
%QPRINT Print shaping window.

% Author: Craig Borghesani
% Date: 3/10/01 3:21PM
% Copyright (c) 2003, Terasoft, Inc.

f = gcf;
bthan=get(f,'userdata');
infmat = get(bthan(16),'userdata');
LoopAxes = infmat(24,1);
bnd_bt = get(bthan(12),'userdata');

switch UIOperation,
   case('Print'),

      PrintFigure = figure('vis','off');
      PrintLoopAxes = copyobj(LoopAxes, PrintFigure);
      set(PrintLoopAxes,'units','norm','pos',get(0,'defaultaxesposition'));
      if ~isempty(bnd_bt),
         bnd_ct = 1;
         for k = 1:size(bnd_bt,1),
%            if strcmp(get(bnd_bt{k,2}(1),'vis'),'on'),
               if length(get(bnd_bt{k,2}(1),'xdata')) > 1,
                  BoundHandle(bnd_ct) = bnd_bt{k,2}(1);

               else
                  BoundHandle(bnd_ct) = bnd_bt{k,2}(2);

               end;%if
               BoundLabel(bnd_ct) = bnd_bt(k,1);
               loccomma = find(BoundLabel{bnd_ct} == ',');
               BoundLabel{bnd_ct} = BoundLabel{bnd_ct}(1:loccomma-1);
               bnd_ct = bnd_ct + 1;
%            end
         end
         legend(PrintLoopAxes, BoundHandle, BoundLabel);
      end
      drawnow;
      set(PrintFigure,'vis','on');

   case('PrintSetup'),

      printdlg('-setup');

end

