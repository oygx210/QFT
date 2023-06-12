function qunitchng
%QUNITCHNG Change from rad/sec to hertz. (Utility)

% Author: Craig Borghesani <cborg@terasoft.com>
% Date: 2/28/02 11:50AM

f=gcf;
bthan = get(f,'userdata');
QFTToolData = getappdata(f,'QFTToolData');
BoundListbox = QFTToolData.Bounds(3);
wbs = get(bthan(11),'userdata');

CurrentRadioButton = gco;
if get(get(CurrentRadioButton,'userdata'),'value'),
   set(CurrentRadioButton,'value',1);
   set(get(CurrentRadioButton,'userdata'),'value',0);

   % update bound unit display
   BoundListboxString = get(BoundListbox,'string');
   FrequencyUnitw2hz = 1/(2*pi);
   FrequencyUnithz2w = 2*pi;
   FrequencyRadioButton = get(CurrentRadioButton,'string');

   if ~isempty(BoundListboxString),
      for k = 1:length(BoundListboxString),
         LocComma = findstr(BoundListboxString{k},',');
         FrequencyValue = str2num(BoundListboxString{k}(1:LocComma-1));
         ColorValue = BoundListboxString{k}(LocComma+1:end);

         if strcmp(FrequencyRadioButton,'Hertz'),
            FrequencyValue = FrequencyValue * FrequencyUnitw2hz;

         else
            FrequencyValue = wbs(k);

         end;%if

         BoundListboxString{k} = sprintf('%-0.4g,%s', FrequencyValue, ColorValue);
      end; %for

      set(BoundListbox,'string',BoundListboxString);

   end;%if

end;%if

