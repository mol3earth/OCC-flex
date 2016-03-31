function meters = ll2m(lats,lons)
% ll2m Finds distance (in meters) between 2 latitude-longitude pairs
% Explicit: Finds distance between lat(1) lons(1) and lat(2) lons(2) 
% example: meters = ll2m(lats,lons)
% Mark Oscar Larson 
%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lats:     2 latitudes
% lons:     2 longitudes 
%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meters:   distance in meters
%%% Internal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R:    Radius of the Earth

lat1=lats(1)*pi/180;
lat2=lats(2)*pi/180;
lon1=lons(1)*pi/180;
lon2=lons(2)*pi/180;
R = 6371000; % meters
dellat = abs(lat1-lat2);
dellon = abs(lon1-lon2);
a = sin(dellat/2) * sin(dellat/2) + ...
         cos(lat1) *  cos(lat2) * ...
        sin(dellon/2) * sin(dellon/2);
c = 2 * atan2( sqrt(a), sqrt(1-a) );
meters = R * c;
