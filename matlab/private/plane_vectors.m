function [normal, slip] = plane_vectors(strike, dip, rake)
%PLANE_VECTORS Unit normal and slip vectors (North-East-Down) of fault planes.
%   [normal, slip] = plane_vectors(strike, dip, rake), angles in degrees.
%   normal points up, out of the footwall; slip is the hanging-wall slip direction.
strike = strike(:); dip = dip(:); rake = rake(:);
normal = [-sind(dip).*sind(strike), sind(dip).*cosd(strike), -cosd(dip)];
slip = [cosd(rake).*cosd(strike) + cosd(dip).*sind(rake).*sind(strike), ...
        cosd(rake).*sind(strike) - cosd(dip).*sind(rake).*cosd(strike), ...
        -sind(rake).*sind(dip)];
end
