%% Latitude/Longitude to TM Coord. Converter

function [PosLon,PosLat]=lin_to_geo(localX, localY, PosRef)

% PosLat : Latitude in Degrees
% PosLon : Longitude in Degrees
% PosRef : Origin of Local coordinate
% lpos : [ localX, localY ] in meter, East-North Coordinate

% Convert Geographic coordinate into Linear coordinate with WGS84
% Ellipsoid model constants (actual values here are for WGS84

R0 = 6378137.0;
E=1/298.257223563;

Rn = R0*(1-E^2)/((1-E^2 * (sind(PosRef(1))^2))^(3/2));
Re = R0/((1-E^2 *(sind(PosRef(1))^2))^(1/2));
Ra = Re*Rn / sqrt( Re^2*sind(PosRef(1))^2 + Rn^2*cosd(PosRef(1))^2 );

deltaLon = localX * 180/pi /Ra / cosd(PosRef(1));
deltaLat = localY * 180/pi /Rn;

PosLon = deltaLon + PosRef(2);
PosLat = deltaLat + PosRef(1);

end
% display('geo_to_lin.m >> DONE')
