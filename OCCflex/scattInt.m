function [ fctSI, mapbu, newmap, x_map , y_map, long, lat] = scattInt(map,lat,long,long1,long2,lat1,lat2)

[nul lnsidx]=min(abs(long+long1));
[nul lnbidx]=min(abs(long+long2));
long = long(lnsidx:lnbidx);
[nul ltsidx]=min(abs(lat-lat1));
[nul ltbidx]=min(abs(lat-lat2));
lat = lat(ltsidx:ltbidx);
map = map(ltsidx:ltbidx,lnsidx:lnbidx);
% create a backup
mapbu = map;
[x_map,y_map]=meshgrid(long,lat);

% do some interpolation, so we can filter
[yyy xxx] = find(isfinite(map)==1);
[yyi xxi] = find(isnan(map));
xxx = long(xxx);
yyy = lat(yyy);
xxi = long(xxi);
yyi = lat(yyi);
zzz = map(find(isfinite(map)==1));
%F = griddedInterpolant(x_bathy',y_bathy',map','spline');
% this interpolation affects the original data IN NO WAY
fctSI = scatteredInterpolant(xxx',yyy',zzz,'natural','nearest');
newmap = map;
newmap(find(isnan(map))) = fctSI(x_map(find(isnan(map))),y_map(find(isnan(map))));
