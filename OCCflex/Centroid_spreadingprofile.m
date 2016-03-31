function [gp1 gp2 AC TopProf BotProf ] = Centroid_spreadingprofile(bathy,axis,slopestruct,lat,long,n,SD)
% Centroid_spreadingprofile Multiple Profiles along spreading direction.
% Usage: Centroid_spreadingprofile(bathy,axis,lat,long,slopestruct,n,SD,LBF,RBF)
% Mark Oscar Larson 2015
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bathy:    bathymetry map, can ontain NaNs
% axis:     matrix of axis segments
% lat:      vector of latitudes, must match size(bathy,2) 
% long:     vector of longitudes, must match size(bathy,1)
% LL:       sturcture of lat long points from which profiles will be
%               collected
% n:        which lat long point in LL is currently being worked on
% SD:       azimuth of spreading direction
% RBF:      right/east bounding fault
% LBF:      left/west bounding fault
% 
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gp1 (good Profile):    the x-axis (lat or long) points of profile
% gp2 (good Profile):    the y-axis (lat or long) points of profile
% AC (axis_crossing):     point on the axis that the main profile crosses
%%%
% ACB (axis_crossing):    point on the axis that the bottom profile crosses
% ACT (axis_crossing):    point on the axis that the top profile crosses
%%%
% TopProf:          A profile parallel to the good profile, but above it
% BotProf:          A profile parallel to the good profile, but below it
%%%  unassigned -mol
% LBF_crossing:     point on the left bounding fault profile intersection
% RBF_crossing:     point on the right bounding fault profile intersection


% first, find which orientation axis is
%if abs(axis(1,1)-axis(end,1)) < abs(axis(1,2)-axis(end,2))
if sum(abs(diff(axis(isfinite(axis(:,1)),1)))) < sum(abs(diff(axis(isfinite(axis(:,2)),2))))
    MA = 1;
    MM = 2;
    l1 = long;
    l2 = lat;
else
    MA = 2;
    MM = 1;
    l1 = lat;
    l2 = long;
end

% set up the stuff
mid = (slopestruct.rps(n).Centroid);
TopMid = slopestruct.rps(n).MajorAxisTop;
BotMid = slopestruct.rps(n).MajorAxisBot;
TopMid1 = TopMid(MA);
TopMid2 = TopMid(MM);
BotMid1 = BotMid(MA);
BotMid2 = BotMid(MM);
mid1ix = round(mid(MA));
mid2ix = round(mid(MM));
mid1 = l1(mid1ix);
mid2 = l2(mid2ix);

% do axes stuff
[ax1 ax2] = find_crossing_line(axis,MA,MM,mid1,mid2ix);
% [LBF1 LBF2] = find_crossing_line(LBF,MA,MM,mid1,mid2ix);
% [RBF1 RBF2] = find_crossing_line(RBF,MA,MM,mid1,mid2ix);

% find which side of ridge starting point is on
% so we multiply our profile distances accordingly
if mid1 > ax1(mid2ix)
    % to da right!
    toda1 = 2.25;
    toda2 = 0.25;
else
    % to da left!
    toda1 = 0.25;
    toda2 = 2.25;
end
%%
sino = sind(90-SD);
coso = cosd(90-SD);
sico(MA) = sino;
sico(MM) = coso;
slope = (tand(90-SD));
% Now other side of profile properties
sino2 = sind(90-SD+180);
coso2 = cosd(90-SD+180);
sico2(MA) = sino2;
sico2(MM) = coso2;
% find y axis intersections
b = mid2 - slope*mid1;
bB = BotMid2 -slope*BotMid1;
bT = TopMid2 -slope*TopMid1;
testy = slope*ax1 + b;
testyB = slope*ax1 + bB;
testyT = slope*ax1 + bT;
% find the index for our axis point
% as well as left and right bounding faults
[nul axidx] = min( abs( abs(testy) - abs(ax2) ) );
[nul axidxB] = min( abs( abs(testyB) - abs(ax2) ) );
[nul axidxT] = min( abs( abs(testyT) - abs(ax2) ) );
% [nul LBFidx] = min( abs( abs(slope*LBF1 + b) - abs(LBF2) ) );
% [nul RBFidx] = min( abs( abs(slope*RBF1 + b) - abs(RBF2) ) );
% find distance to this point
length2axis = ll2m([ax1(axidx) mid1],[ax2(axidx) mid2]);
hypo = sqrt((ax1(axidx) - mid1)^2 + (ax2(axidx) - mid2)^2);

%%
% figure(10)
% clf
%  
% imshade(long,lat, bathy);
% %surf(long,lat,bathybu)
% hold on
% plot(ax1,ax2,'-','Color',[.8 .8 .8],'LineWidth',1.3)
% shading interp
% view([0 90])
% plot(ax1,testy,'--k')
% plot(mid1,mid2,'b*')
% plot(ax1(axidx),ax2(axidx),'r*')
% %%
% plot(gp1,gp2,'-r')

%%
% make sure profile is sufficient length
% if length2axis*2<minproflength
%     hypo = hypo*(minproflength/length2axis);
% end   
extentA = [sign(slope) -1].*hypo.*toda1.*sico + [mid2 mid1];
extBotA = [sign(slope) -1].*hypo.*toda1.*sico + [BotMid2 BotMid1];
extTopA = [sign(slope) -1].*hypo.*toda1.*sico + [TopMid2 TopMid1];
% and other side of profile 
extentB = [sign(slope) -1].*hypo.*toda2.*sico2 + [mid2 mid1];  
extBotB = [sign(slope) -1].*hypo.*toda2.*sico2 + [BotMid2 BotMid1];     
extTopB = [sign(slope) -1].*hypo.*toda2.*sico2 + [TopMid2 TopMid1];


% make the profiles contain 'sufficiently' many points
% find out how big?
% maybe later
% [ nul idx1 ] = min( abs( abs(l1) - abs(extentA(2)) ) );
% [ nul idx2 ] = min( abs( abs(l1) - abs(extentA(2)) ) );
% [ nul idx3 ] = min( abs( abs(l1) - abs(extentA(2)) ) );
% [ nul idx4 ] = min( abs( abs(l1) - abs(extentA(2)) ) );
%
clear good_profile
% since ll2m([13 13],[-44 -44.001845955]) == 200.0000
% use spacing increment of 0.001845955 
%%% Old way
% good_profile(:,2) = linspace(extentA(1),extentB(1),1000);
% good_profile(:,1) = linspace(extentA(2),extentB(2),1000);
% BotProf(:,2) = linspace(extBotA(1),extBotB(1),1000);
% BotProf(:,1) = linspace(extBotA(2),extBotB(2),1000);
% TopProf(:,2) = linspace(extTopA(1),extTopB(1),1000);
% TopProf(:,1) = linspace(extTopA(2),extTopB(2),1000);
try good_profile(:,2) = extentA(1):0.00001845955:extentB(1) ;
catch good_profile(:,2) = fliplr(extentB(1):0.00001845955:extentA(1) );
end
good_profile(:,1) = linspace(extentA(2),extentB(2),length(good_profile(:,2)));

try BotProf(:,2) = extBotA(1):0.00001845955:extBotB(1);
catch BotProf(:,2) = fliplr(extBotB(1):0.00001845955:extBotA(1));
BotProf(:,1) = linspace(extBotA(2),extBotB(2),length(BotProf(:,2)));
end
try TopProf(:,2) = extTopA(1):0.00001845955:extTopB(1);
catch TopProf(:,2) = fliplr(extTopB(1):0.00001845955:extTopA(1));
TopProf(:,1) = linspace(extTopA(2),extTopB(2),length(TopProf(:,2)));
end

% check if all inside for profile
l2lim = [l2(1) l2(end)];
l1lim = [l1(1) l1(end)];
good_profile = lineinside(l1lim,l2lim,good_profile);
BotProf = lineinside(l1lim,l2lim,BotProf);
TopProf = lineinside(l1lim,l2lim,TopProf);


% set output
AC = [ax1(axidx) ax2(axidx)];
ACB = [ax1(axidxB) ax2(axidxB)];
ACT = [ax1(axidxT) ax2(axidxT)];
% RBFcrossing = [RBF1(axidx) RBF2(axidx)];
% LBFcrossing = [LBF1(axidx) LBF2(axidx)];
gp1 = good_profile(:,1);
gp2 = good_profile(:,2);
% test 
% findbestprofilelength(axis_interp,slopeE_Wcc,lat,long,469)
