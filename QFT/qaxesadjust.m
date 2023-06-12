function [TickValues] = qaxesadjust(CurrentLims)
%QAXESADJUST Adjust axes tick labels to 45 degree increments. (Utility)

% Author: Craig Borghesani <cborg@terasoft.com>
% Date: 2/27/02 11:26PM

if abs(diff(CurrentLims)) >= 360,

   MinRem = rem(CurrentLims(1), 45);
   MaxRem = rem(CurrentLims(2), 45);
   if MinRem ~= 0 & MaxRem ~= 0,
      StartValue = CurrentLims(1) - MinRem;
      EndValue = CurrentLims(2) + MaxRem;
      TickValues = [StartValue:45:EndValue];

   elseif MinRem ~= 0,
      StartValue = CurrentLims(1) - MinRem;
      EndValue = CurrentLims(2);
      TickValues = [StartValue:45:EndValue];

   elseif MaxRem ~= 0,
      StartValue = CurrentLims(1);
      EndValue = CurrentLims(2) + MaxRem;
      TickValues = [StartValue:45:EndValue];

   else
      StartValue = CurrentLims(1);
      EndValue = CurrentLims(2);
      TickValues = [StartValue:45:EndValue];

   end % if MinRem ~= 0 & MaxRem ~= 0

end % if diff(CurrentLims) >= 360

