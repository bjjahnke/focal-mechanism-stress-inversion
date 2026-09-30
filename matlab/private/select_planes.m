function [normals, slips] = select_planes(planes, choice)
%SELECT_PLANES Normals and slips of the chosen plane (1 or 2) for each event.
use2 = choice(:) == 2;
normals = planes.normal1; normals(use2, :) = planes.normal2(use2, :);
slips = planes.slip1;     slips(use2, :) = planes.slip2(use2, :);
end
