%CALCULATE SLOPE DISTRIBUTION FOR A REGION

% give the name of the grid file
gridfile = 'Bathy12_15N.grd';
% an axis file if you have on, otherwise, you will pick the axis points
Axisfile ='Axis_12_15N_Mallow&Searle2012_supp_Fujiwara2003.txt'; 
% East and west bounding fault files if you have them,
WestBoundingFaults = 'WBF_Mallow&Searle2012_Aug12.txt';
EastBoundingFaults = 'EBF_Mallow&Searle2012_Aug12.txt';
% a gravity map of the region if you have it.
do_grav = 1; % binary switch
gravgridfile = 'cutrmba.grd';
gravtype ='RMBA'; % RMBA or MBA, etc
% if there is EQ data
EQfname = 'EQ_GMT_Aug18_2015.m'; % if no data, set to ''

% set the subregion, if you want
% otherwise, region will be the entire grid file
long1 =  45.6;
long2 = 44.4; 
lat1 = 12.63;
lat2 =  15.25;

% Spreading direction, I got this from Okino's plate calculator
SD = 102.9;
% these were the coordinates for my spreading direction
% (lat,lon) = (30.10, 317.90) is 23.64 mm/year at an azimuth of N102.9E

% inputs for fft filtering of the region
% set the low cut/highpass filter here
% highpass1 will assign all wavelengths greater than itself a '1' in fft
% domain
% highpass2 will grade sinusoidally from 1 (at highpass1) to 0, at
% highpass2
highpass1x = 100; % these are in units of km roughly
highpass2x = 50;
% Set y direction to zero to uniformly filter 
highpass1y = 100;
highpass2y = 50;
% do the filter swith
dohpf =1;
% do the filter twice switch? this removes the large features
removebig=1;
if dohpf 
    if highpass2y > 0
        HPFstr = sprintf('HPFxy_%s_%s_%s_%s',num2str(highpass1x),...
                                             num2str(highpass2x),...
                                             num2str(highpass1y),...
                                             num2str(highpass2y));
    else
        HPFstr = sprintf('HPFx%s_%s',num2str(highpass1x),num2str(highpass2x));
    end
else
   HPFstr = 'noHPF';
end

% set cutoff area (in pixels)
% or cutoff length of major axis
% set the other value to 0
cutoff=0;
cutofflength = 30;
cutofflength2 = 5;
connectivity=8;
if cutoff > 0
    CLstring = sprintf('CL%d',cutoff);
else
    CLstring = sprintf('CL%d_%d',cutofflength2,cutofflength);
end
% set slope filter amounts in degrees
lowslope = 15;
highslope = 70;

% set profile length in km
minimum_prof_length = 80;
% set to 0, to let variable length be decided
orthog_to_ridge = 0;

%%% Save Figures? 
print_y_n = 0;

% This is a cuttoff for direction of identified features facing
% give the average spreading direction, 
% or whatever you want for west azimuth spreading
% direction
% east will be -180 of that
westaz = 270; 
eastaz = westaz-180;
% give two values bounding the east direction
eastlowaz=40;
easthighaz=180;
% and two values bounding the west direction
westlowaz =180;
westhighaz=350;

% make title for figures from provided parameters
disp(sprintf('\n\n:::: Parameters are understood to be ::::'))
app_title = sprintf('%s_%s_SL%.1f_%.1f_Waz%d_%d_Eaz%d_%d_%s',...
            HPFstr,...
            CLstring,...
            lowslope,highslope,...
            westlowaz,westhighaz,...
            eastlowaz,easthighaz,...
            datestr(now,'mm_dd_yyyy'));
disp(app_title)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up interpolation for Bathy and Grav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for bathymetry
disp(sprintf(':::: Loading :::: \n \t\t %s ',gridfile))
[long,lat, bathy]=grdread2(gridfile);
disp(':::: File Grid Spacing :::: :::: :::: :::: :::: :::: :::: :::: :::: :::: ::::')
disp(sprintf('\t :: x dimension:\n\t\tMax: %.2f Deg \n\t\tMin: %.2f Deg \n\t\tSpacing:\t%.4f Deg \n\t\tLat %.2f:\t%.4f m\n\t\tLat %.2f:\t%.4f m',   max(long), min(long), abs(long(1)-long(2)), min(lat),  ll2m([min(lat) min(lat)],[long(1) long(2)]), max(lat),  ll2m([max(lat) max(lat)],[long(1) long(2)]) ))
disp(sprintf('\t :: y dimension:\n\t\tMax: %.2f Deg \n\t\tMin: %.2f Deg \n\t\tSpacing:\t%.4f Deg \n\t\t\t\t%.4f m', max(lat),  min(lat),  abs(lat(1)-lat(2)),ll2m([lat(1) lat(2)],[long(1) long(1)]) ))
disp(sprintf('\t :: z dimension:\n\t\tMax: %.2f m \n\t\tMin: %.2f m ',max(max(bathy)),min(min(bathy))))
disp(':::: :::: :::: :::: :::: :::: :::: :::: :::: :::: :::: :::: :::: :::: :::::::')
bathy = double(bathy);
% bathy = bathy(1:5:end,1:5:end);
% long = long(1:5:end);
% lat=lat(1:5:end);
%% do interpolation, this is because some functions do not deal with NaNs
[ F , bathybu , bathy , x_bathy, y_bathy, long, lat ] = scattInt(bathy,lat,long,long1,long2,lat1,lat2);

% plot it up, for visual comparison
figure(99)
clf
subplot(2,4,1)
imshade(long,lat, bathy);colorbar('SouthOutside')
title('Interpolated Bathy')
subplot(2,4,2)
imshade(long,lat, bathybu);colorbar('SouthOutside')
title('Raw Bathy')
%% Now do for gravity

if do_grav
    % import grav file
    [Glong,Glat, gravy]=grdread2(gravgridfile);
    gravy = double(gravy);
    % interpolate it
    [ G , gravybu , gravy, ~ , ~ , Glong, Glat ] = scattInt(gravy,Glat,Glong,long1,long2,lat1,lat2);
    % now interpolate to make gravity map same size as bathymetry
    gravy = G(x_bathy,y_bathy);
    % plot the interpolated and the non
    figure(99)
    subplot(2,4,3)
    imshade(long,lat, gravy);colorbar('SouthOutside')
    title(sprintf('Interpolated %s',gravtype))
    subplot(2,4,4)
    imshade(Glong,Glat, gravybu);colorbar('SouthOutside')
    title(sprintf('Raw %s',gravtype))
    % now fft filter it
    HPF_gravy = fftfiltermap(gravy,lat,long,.4,.2);
    figure(99)
    subplot(2,4,8) 
    imshade(Glong,Glat, HPF_gravy);colorbar('SouthOutside')
    title(sprintf('High Pass Filtered %s',gravtype))
end

%%
if dohpf
    if highpass2y > 0
        if removebig == 1
                HPF_bathy_forrem = directional_fftfiltermap(bathy,lat,long,...
                    highpass1x,highpass2x,highpass1y,highpass2y);
                bathy2 = bathy - HPF_bathy_forrem;
                HPF_bathy = directional_fftfiltermap(bathy2,lat,long,...
                    lowpass1x,lowpass2x,lowpass1y,lowpass2y);
                
        else
                HPF_bathy = directional_fftfiltermap(bathy,lat,long,...
                    highpass1x,highpass2x,highpass1y,highpass2y);
            
        end
                
                
    else
        HPF_bathy = fftfiltermap(bathy,lat,long,highpass1x,highpass2x);
    end
    %bathy = HPF_bathy; 
    Fbathy = bathy - HPF_bathy;
    figure(99)
    subplot(2,4,5)
    imshade(long,lat, Fbathy);colorbar('SouthOutside')
    title('Band Pass Filtered Bathy')    
end
%% Now plot the results and compare to original
figure(999)
subplot(1,2,1)
imshade(long,lat, bathy);
title('Original')
subplot(1,2,2)
imshade(long,lat, Fbathy);
title(sprintf('Band Pass Filtered Bathy \n %.0f km - %.0f km EW & %.0f km - %.0f km NS',highpass2x,lowpass1x,highpass2y,lowpass1y)  )
%% do the gradient
% using matlab mapping toolbox function 'gradientm'
% this needs the lat and long, and expect z data in meters
[aspect, slope, gradN, gradE] = gradientm(y_bathy,x_bathy, bathy);

%dirslope = gradN*(1-90/SD) + gradE*(90/SD);
dirslope = atand(gradN*cosd(SD) + gradE*sind(SD));
slopebu = slope;
slope = dirslope;
[ S ] = scattInt(dirslope,lat,long,long1,long2,lat1,lat2);

[Gaspect, Gslope, GgradN, GgradE] = gradientm(y_bathy,x_bathy, gravy);

% visual check of min max
%min(min(slope))
%max(max(slope))
%min(min(gradN))
%max(max(gradN))
%% plot pretty
dirslope = atand(gradN*cosd(SD)+gradE*sind(SD));
fig104=figure(104);
clf
subplot(1,2,1)
surf(long,lat,slopebu.*isfinite(bathybu));
view(0,90)
axis equal
axis tight
shading interp
% xlim([-44.98 -44.8])
% ylim([13.43 13.58 ])
hold on
title('Slope Map')
%title('Azimuth of Gradient','Fontsize',28,'Fontname','Ubuntu')
subplot(1,2,2)
surf(long,lat,abs(dirslope).*isfinite(bathybu));
view(0,90)
axis equal
axis tight
shading interp
% xlim([-44.98 -44.8])
% ylim([13.43 13.58 ])
hold on
title('Directional Slope Map')


fig104.PaperPositionMode = 'auto';
if print_y_n
print(fig104,sprintf('Azimuth_Prefilter_%s.png',sprintf('%s_%s',...
            HPFstr,...
            datestr(now,'mm_dd_yyyy'))),'-dpng','-r0')
end

%% Now take care of the axis
try % try loading an axis file, and plot it
    axis1=load(Axisfile);
    
    Ax = axis1(:,1);
    Ay = axis1(:,2);
    disp(sprintf('<<<<<<<<< Using Axis data from %s',Axisfile))
    disp('<<<<<<<<<< <<<<<<<<< >>>>>>>>> >>>>>>>>>>')
    fig104=figure(104);
    clf
    surf(long,lat,bathybu);
    shading interp
    view(0 ,90)
    axis equal
    axis tight
    hold on
    plot(Ax,Ay,'r*')
    plot(Ax,Ay,'r--')
catch % otherwise, get user input
    fig104=figure(104);
    clf
    surf(long,lat,bathybu);
    shading interp
    view(0 ,90)
    lightangle(SD,1e-5)
    lightangle(SD+180,1e-5)
    axis equal
    axis tight
    hold on
    disp('<<<<<<<<< Select Axis by clicking along it')
    disp('<<<<<<<<< Press the Return key when done selecting points >>>>>>>>')
    [Ax,Ay] = ginput;
    plot(Ax,Ay,'r*')
    disp('<<<<<<<<<< Axis Accepted, Thank you >>>>>>>>>>')
    disp('<<<<<<<<<< <<<<<<<<<<<< >>>>>>>>>>> >>>>>>>>>>')
    savethisaxis = input(sprintf('Want to save this as ''%s'' ?\n(1 for yes, 0 for no)\n  ',Axisfile))
    if savethisaxis 
        file = [Ax,Ay];
        save(Axisfile,'file', '-ascii')
    end        
end


%% plot EQs, 
if strcmp(EQfname,'') ~= 1
    plotEQs(EQfname,lat1,lat2,long1,long2,SD)
end

%% filter out slopes with 'extreme slopes' & filter out slopes that face certain directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create new copies, so we don't mess with original data
slope2=slope;
gradE2=gradE; 
gradN2 =gradN;
aspect2=aspect;

disp(sprintf('<<<<<<<<<< Slope Filter: Removing %.0f elements with slope greater than %.1f',numel(find(slope>highslope)),highslope))
disp(sprintf('<<<<<<<<<<\tand %.0f elements with slope less than %.1f\n',numel(find(slope<lowslope)),lowslope))
yes=find(abs(slope)>highslope | abs(slope) < lowslope);
slope2(yes)=NaN;
aspect2(yes)=NaN;
figure(101)
clf
subplot(1,2,1)
surf( long, lat, slope2);
hold on
plot(Ax,Ay,'r-')
colorbar('SouthOutside')
title('''Filtered by Aspect'' Slope Map')
shading interp
view(0 ,90)
axis equal
%%
figure(101)
clf
subplot(1,2,1)
surf(slope2);
colorbar('SouthOutside')
title('''Filtered by Aspect'' Slope Map')
shading interp
view(0 ,90)
axis equal
xlim([0 size(slope2,2)])
ylim([0 size(slope2,1)])

%% remove slopes that face the wrong ways
%%% 'Crude' Catch
if dohpf
disp('<<<<<<<<<< Crude Boundary Filter: Getting rid of boundary effected slopes:')
    disp(sprintf('<<<<<<<<<<\tWhich Account for %.1f %% of total Identified Slopes\n',100*sum(   [numel(isfinite(slope2(1:5,1:end)))    numel(isfinite(slope2((end-5):end,1:end)))    numel(isfinite(slope2(1:end,1:5)))   numel(isfinite(slope2(1:end,(end-5):end)))]   )/numel(isfinite(slope2))))
slope2(1:5,1:end) =NaN; slope2((end-5):end,1:end) = NaN; slope2(1:end,1:5) = NaN; slope2(1:end,(end-5):end) = NaN;
aspect2(1:5,1:end) =NaN; aspect2((end-5):end,1:end) = NaN; aspect2(1:end,1:5) = NaN; aspect2(1:end,(end-5):end) = NaN;
gradN2(1:5,1:end) =NaN; gradN2((end-5):end,1:end) = NaN; gradN2(1:end,1:5) = NaN; gradN2(1:end,(end-5):end) = NaN;
gradE2(1:5,1:end) =NaN; gradE2((end-5):end,1:end) = NaN; gradE2(1:end,1:5) = NaN; gradE2(1:end,(end-5):end) = NaN;
end
slope2 = slope2.*isfinite(bathybu);
aspect2 = aspect2.*isfinite(bathybu);
gradN2 = gradN2.*isfinite(bathybu);
gradE2 = gradE2.*isfinite(bathybu);

%%%
% figure(3);clf;surf((aspect2));view([0 90]);shading interp
% colorbar
%SPREADING DIRECTION IS 87 DEGREES
if exist('slopebu')==0
    % find east facers
    notE = find(aspect2<eastlowaz | aspect2>easthighaz);
    disp(sprintf('<<<<<<<<<< Azimuth East Filter: Removing %.0f elements which face not between Azimuth %.1f & %.1f',numel(notE),easthighaz,eastlowaz))
    disp(sprintf('<<<<<<<<<<\tWhich are %.1f & %.1f degrees different from the spreading direction of %.1f\n',eastaz-eastlowaz,eastaz-easthighaz,eastaz))
    % now for west facers
    notW = find(aspect2<westlowaz | aspect2>westhighaz);
    disp(sprintf('<<<<<<<<<< Azimuth West Filter: Removing %.0f elements which face not between Azimuth %.1f & %.1f',numel(notW),westhighaz,westlowaz))
    disp(sprintf('<<<<<<<<<<\tWhich are %.1f & %.1f degrees different from the spreading direction of %.1f\n',westaz-westlowaz,westaz-westhighaz,westaz))
else
    % find east facers
    notW = find(slope2<0);
    disp(sprintf('<<<<<<<<<< Azimuth East Filter: Removing %.0f elements which face West',numel(notW)))
    % now for west facers
    notE = find(slope2>0);
    disp(sprintf('<<<<<<<<<< Azimuth West Filter: Removing %.0f elements which face East',numel(notE)))
end
    
% allocate the new matrices
slopeE = abs(slope2);
slopeW = abs(slope2);
aspectE=aspect2;
aspectW=aspect2;
% assign values to them
slopeE(notE) = NaN;
slopeW(notW) = NaN;
aspectE(notE)=NaN;
aspectW(notW)=NaN;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(1);
clf
%Debbie: imagesc(long,lat,slopeE);axis xy; title('East facing slopes')
surf(long,lat,slopeE)
hold on
plot(Ax,Ay,'k*','LineWidth',3)
view([0 90])
shading interp
axis equal
xlim([long(1) long(end)])
ylim([lat(1) lat(end)])
title('East facing slopes')
% xlabel('Longitude (Degrees)')
% ylabel('Latitude (Degrees)')
cb=colorbar('Location','South','Fontsize',10);
cbPos = cb.Position;cbPos(4) = .5*cbPos(4);cb.Position = cbPos;
set(gca,'Fontsize',18,'Fontname','Ubuntu')
hold off
fig1.PaperPositionMode = 'auto';
if print_y_n
print(fig1,sprintf('E_Slopes__%s.png',sprintf('%s_SL%.1f_%.1f_Waz%d_%d_Eaz%d_%d_%s',...
            HPFstr,...
            lowslope,highslope,...
            westlowaz,westhighaz,...
            eastlowaz,easthighaz,...
            datestr(now,'mm_dd_yyyy'))),'-dpng','-r0')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2=figure(2);
clf
%Debbie: imagesc(long,lat,slopeW); axis xy; title ('West facing slopes');
surf(long,lat,slopeW)
hold on
plot(Ax,Ay,'k--','LineWidth',4)
view([0 90])
shading interp
axis equal
xlim([long(1) long(end)])
ylim([lat(1) lat(end)])
title('West facing slopes')
% xlabel('Longitude (Degrees)')
% ylabel('Latitude (Degrees)')
cb=colorbar('Location','South','Fontsize',10);
cbPos = cb.Position;cbPos(4) = .5*cbPos(4);cb.Position = cbPos;
set(gca,'Fontsize',18,'Fontname','Ubuntu')
hold off
% write these to file, if want
% grdwrite2(long,lat, slopeE,'AllEast20_60.grd');
% grdwrite2(long,lat, slopeW, 'AllWest20_60.grd');
fig2.PaperPositionMode = 'auto';
if print_y_n
print(fig2,sprintf('W_Slopes__%s.png',sprintf('%s_SL%.1f_%.1f_Waz%d_%d_Eaz%d_%d_%s',...
            HPFstr,...
            lowslope,highslope,...
            westlowaz,westhighaz,...
            eastlowaz,easthighaz,...
            datestr(now,'mm_dd_yyyy'))),'-dpng','-r0')
end


%% seperate axis into overlapping parts and create new axis matrix
axis1 = axisbreaker(axis1);
% intitialize these
slopeE_Wside=slopeE;
slopeE_Eside=slopeE;
slopeW_Wside=slopeW;
slopeW_Eside=slopeW;

clear axis_interp
[axis_interp,slopeE_Eside,slopeW_Wside,slopeE_Wside,slopeW_Eside] = axisinterpolater(axis1,long,lat,slopeE_Eside,slopeW_Wside,slopeE_Wside,slopeW_Eside);

%% plot new axis up
fig3 = figure(3);
clf
hold on
surf(long,lat,slopeW_Wside)
for m = 1:2:size(axis_interp,2)
    plot(axis_interp(:,m),axis_interp(:,m+1),'-','Color',[.5 .5 .5],'LineWidth',1.3)
end
view([0 90])
shading interp
axis equal
xlim([long(1) long(end)])
ylim([lat(1) lat(end)])
title('West Facing Slopes West of Ridge Axis')
% xlabel('Longitude (Degrees)')
% ylabel('Latitude (Degrees)')
cb=colorbar('Location','South','Fontsize',10);
cbPos = cb.Position;cbPos(4) = .5*cbPos(4);cb.Position = cbPos;
set(gca,'Fontsize',18,'Fontname','Ubuntu')
fig3.PaperPositionMode = 'auto';
if print_y_n
print(fig3,sprintf('W_W_Slopes__%s.png',sprintf('%s_SL%.1f_%.1f_Waz%d_%d_Eaz%d_%d_%s',...
            HPFstr,...
            lowslope,highslope,...
            westlowaz,westhighaz,...
            eastlowaz,easthighaz,...
            datestr(now,'mm_dd_yyyy'))),'-dpng','-r0')
end
%% plot up slopes 
fig4 = figure(4);
clf
surf(long,lat,slopeE_Eside)
hold on
for m = 1:2:size(axis_interp,2)
    plot(axis_interp(:,m),axis_interp(:,m+1),'-','Color',[.5 .5 .5],'LineWidth',1.3)
end
    view([0 90])
shading interp
axis equal
xlim([long(1) long(end)])
ylim([lat(1) lat(end)])
title('East Facing Slopes East of Ridge Axis')
cb=colorbar('Location','South','Fontsize',10);
cbPos = cb.Position;cbPos(4) = .5*cbPos(4);cb.Position = cbPos;
set(gca,'Fontsize',18,'Fontname','Ubuntu')
fig4.PaperPositionMode = 'auto';
if print_y_n
print(fig4,sprintf('E_E_Slopes__%s.png',sprintf('%s_SL%.1f_%.1f_Waz%d_%d_Eaz%d_%d_%s',...
            HPFstr,...
            lowslope,highslope,...
            westlowaz,westhighaz,...
            eastlowaz,easthighaz,...
            datestr(now,'mm_dd_yyyy'))),'-dpng','-r0')
end

clear slopeW_Wcc slopeE_Wcc slopeW_Ecc slopeE_Ecc
slopeW_Wcc = bwconncomp(slopeW_Wside>0,connectivity);
slopeW_Ecc = bwconncomp(slopeW_Eside>0,connectivity);
slopeE_Ecc = bwconncomp(slopeE_Eside>0,connectivity);
slopeE_Wcc = bwconncomp(slopeE_Wside>0,connectivity);
% which properties do you want?
propsofinterest = {  'Centroid','Area','Orientation',...
                     'MajorAxisLength','MinorAxisLength',...
                     'Eccentricity','Orientation'};
slopeE_Wcc.rps = regionprops(slopeE_Wcc,propsofinterest);
slopeW_Wcc.rps = regionprops(slopeW_Wcc,propsofinterest);
slopeE_Ecc.rps = regionprops(slopeE_Ecc,propsofinterest);
slopeW_Ecc.rps = regionprops(slopeW_Ecc,propsofinterest);
if cutoff>0
    disp('<<<<<<<<<< Using Cutoff >>>>>>>>>>')
    disp(sprintf('\tCutoff is %.0f\n\tCutting:',cutoff))
    disp(sprintf('\t\t%.0f West Facing Slopes, East of Axis',numel(find([slopeW_Ecc.rps.Area]>cutoff))))
    disp(sprintf('\t\t%.0f East Facing Slopes, West of Axis',numel(find([slopeE_Ecc.rps.Area]>cutoff))))
    disp(sprintf('\t\t%.0f West Facing Slopes, West of Axis',numel(find([slopeW_Wcc.rps.Area]>cutoff))))
    disp(sprintf('\t\t%.0f East Facing Slopes, West of Axis',numel(find([slopeE_Wcc.rps.Area]>cutoff))))
    disp('<<<<<<<<<< <<<<<<<<< >>>>>>>>> >>>>>>>>>>')
elseif cutofflength>0
    disp('<<<<<<<<<< Using CutoffLength >>>>>>>>>>')
    disp(sprintf('\tCutoffLength is %.0f\n\tCutting:',cutofflength))
    disp(sprintf('\t\t%.0f West Facing Slopes, East of Axis',numel(find([slopeW_Ecc.rps.MajorAxisLength]<cutofflength))))
    disp(sprintf('\t\t%.0f East Facing Slopes, West of Axis',numel(find([slopeE_Ecc.rps.MajorAxisLength]<cutofflength))))
    disp(sprintf('\t\t%.0f West Facing Slopes, West of Axis',numel(find([slopeW_Wcc.rps.MajorAxisLength]<cutofflength))))
    disp(sprintf('\t\t%.0f East Facing Slopes, West of Axis',numel(find([slopeE_Wcc.rps.MajorAxisLength]<cutofflength))))
    disp('<<<<<<<<<< <<<<<<<<< >>>>>>>>> >>>>>>>>>>')
end
%%
fig6 = figure(6);
clf
hold on
large_slopeW_Wside = slopeW_Wside*NaN;
for n=1:slopeW_Wcc.NumObjects
    if cutoff>0
    if slopeW_Wcc.rps(n).Area>cutoff
        large_slopeW_Wside(slopeW_Wcc.PixelIdxList{n})=slopeW_Wside(slopeW_Wcc.PixelIdxList{n});
    end
    elseif cutofflength>0 
    if slopeW_Wcc.rps(n).MajorAxisLength>cutofflength & slopeW_Wcc.rps(n).MinorAxisLength>cutofflength2
        large_slopeW_Wside(slopeW_Wcc.PixelIdxList{n})=slopeW_Wside(slopeW_Wcc.PixelIdxList{n});
    end 
    end
end
surf(long,lat,large_slopeW_Wside)
hold on
for m = 1:2:size(axis_interp,2)
    plot(axis_interp(:,m),axis_interp(:,m+1),'--','Color',[.2 .2 .2],'LineWidth',2)
end
view([0 90])
shading interp
axis equal
xlim([long(1) long(end)])
ylim([lat(1) lat(end)])
title('West Facing Slopes West of Ridge Axis')
cb=colorbar('Location','South','Fontsize',10);
cbPos = cb.Position;cbPos(4) = .5*cbPos(4);cb.Position = cbPos;
set(gca,'Fontsize',18,'Fontname','Ubuntu')
hold off
fig6.PaperPositionMode = 'auto';
if print_y_n
print(fig6,sprintf('W_W_Slopes__%s.png',app_title),'-dpng','-r0')
end
%%
fig7 = figure(7);
clf
hold on
large_slopeE_Eside = slopeE_Eside*NaN;
for n=1:slopeE_Ecc.NumObjects
    if cutoff>0
    if slopeE_Ecc.rps(n).Area>cutoff
        large_slopeE_Eside(slopeE_Ecc.PixelIdxList{n})=slopeE_Eside(slopeE_Ecc.PixelIdxList{n});
    end
    elseif cutofflength>0
    if slopeE_Ecc.rps(n).MajorAxisLength>cutofflength & slopeE_Ecc.rps(n).MinorAxisLength>cutofflength2
        large_slopeE_Eside(slopeE_Ecc.PixelIdxList{n})=slopeE_Eside(slopeE_Ecc.PixelIdxList{n});
    end  
    end
end
surf(long,lat,large_slopeE_Eside)
hold on
for m = 1:2:size(axis_interp,2)
    plot(axis_interp(:,m),axis_interp(:,m+1),'-','Color',[.5 .5 .5],'LineWidth',1.3)
end
view([0 90])
shading interp
axis equal
xlim([long(1) long(end)])
ylim([lat(1) lat(end)])
title('East Facing Slopes East of Ridge Axis')
cb=colorbar('Location','South','Fontsize',10);
cbPos = cb.Position;cbPos(4) = .5*cbPos(4);cb.Position = cbPos;
set(gca,'Fontsize',18,'Fontname','Ubuntu')
hold off
fig7.PaperPositionMode = 'auto';
if print_y_n
print(fig7,sprintf('E_E_Slopes__%s.png',app_title),'-dpng','-r0')
end
%%  Get user picks
% first plot up the map
cmin = min(min(bathybu));
cmax = max(max(bathybu));

clear s s1
cmaplim=64;
figure(9)
clf
s=surf(long,lat,bathybu,'CData',round((-1+cmaplim)*(bathybu-cmin)/(cmax-cmin))+1,'cdatamapping','direct');
hold on
s1(1)=surf(long,lat,large_slopeE_Eside*0,'CData',large_slopeE_Eside*0+(cmaplim+1),'cdatamapping','direct','Facealpha',0.5);%,'cdatamapping','direct',
s1(2)=surf(long,lat,large_slopeW_Wside*0,'CData',large_slopeW_Wside*0+(cmaplim+1),'cdatamapping','direct','Facealpha',0.5);
hold on

set(gcf,'Colormap',[parula(cmaplim);[0 0 0];[0 1 0];[1 0 0]])
view(0,90)
axis equal
axis tight
shading interp
for m = 1:2:size(axis_interp,2)
    plot(axis_interp(:,m),axis_interp(:,m+1),'-r','LineWidth',1.3)
end
lightangle(SD+90,1e-100)
lightangle(SD-90,1e-100)
s.FaceLighting = 'gouraud';
s.AmbientStrength = 0.2;
s.DiffuseStrength = .3;
s.SpecularStrength = .1;
s.SpecularExponent = .1;
s.BackFaceLighting = 'lit'; %'reverselit' | 'unlit' | 'lit'
set(gca,'Fontsize',20,'Fontname','Ubuntu')
title('Bathymetry Overlain by Filtered Slopes')

%% now ask for points

% these OCC points are mine, set to 0 to pick your own
if 1
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Using Mark Larson Picks')
    disp('Go to lime 647 and set to ''0'' to pick your own')
    CCs = [ -45.0262    15.08
            -45.1       14.875
            -44.93      14.843
            -44.9100    14.6650
            -44.8396    13.8339
            -45.0800    13.6980
            -44.9500    13.5150
            -44.9500    13.3200
            -44.9000    13.1040
            -44.7053    12.8254];
else
    [Cx,Cy] = ginput;
    CCs = [Cx Cy];
end

%% get bounding faults
% because the bounding faults are not that important for manual picks, 
% you don't need to always 
if input('Use axis as EBF and WBF for now? ( ''1'' for YES else NO) :\n       ')
    Ebf_interp = axis_interp;
    Wbf_interp = axis_interp;
else
try 
    Wbf=load(WestBoundingFaults);
    Ebf=load(EastBoundingFaults);
    disp(sprintf('<<<<<<<<< Using Bounding Fault data from %s & %s',WestBoundingFaults,EastBoundingFaults))
    disp('<<<<<<<<<< <<<<<<<<< >>>>>>>>> >>>>>>>>>>')
    fig104=figure(104)
    clf
    imshade(long,lat, bathybu);
    hold on
    plot(Ebf(:,1),Ebf(:,2),'r*')    
    plot(Wbf(:,1),Wbf(:,2),'b*')
    Ebf = axisbreaker(Ebf);
    Ebf_interp = interpolateaxis(lat,long,Ebf);
    Wbf = axisbreaker(Wbf);
    Wbf_interp = interpolateaxis(lat,long,Wbf);
catch
    GetBoundingFault
    if input('Want to save these bounding faults?')
        disp('!!!!!!!!! Don''t forget the quotes !!!!!!!!!!!')
        WBFfilename = input(sprintf('Name for West Bounding Fault?\n:'));
        save(WBFfilename,'Wbf_interp','-ascii')
        disp('!!!!!!!!! Don''t forget the quotes !!!!!!!!!!!')
        EBFfilename = input(sprintf('Name for East Bounding Fault?\n:'));
        save(EBFfilename,'Ebf_interp','-ascii')
        clear EBFfilename WBFfilename
    end
end    
end

clear slopes

fig8=figure(8);
clf 
%imshade(long,lat, bathybu);
s=surf(long,lat,bathybu);
% lightangle(100,1e-5)
% lightangle(-80,1e-5)
% light
shading interp
view([0 90])
lightangle(90,1e-2)
lightangle(-60,1e-2)
s.FaceLighting = 'flat';
s.AmbientStrength = 0.3;
s.DiffuseStrength = .6;
s.SpecularStrength = 01;
s.SpecularExponent = 1;
s.BackFaceLighting = 'lit'; %'reverselit' | 'unlit' | 'lit'
hold on
for m = 1:2:size(axis_interp,2)
    plot(axis_interp(:,m),axis_interp(:,m+1),'-m','LineWidth',1.3)
end
shading interp
view([0 90])
axis equal
axis tight
%title(sprintf('Major and Minor Axis of Filtered Slopes \nw/Profile Numbers and profiles drawn'))
%  xlabel('Longitude (Degrees)')
%  ylabel('Latitude (Degrees)')
title('Profile Locations')
 set(gca,'Fontsize',18,'Fontname','Ubuntu')
%%%
H = gca;
for hs = 1:length(H.XTickLabel)
    dotidx = strfind(H.XTickLabel{hs},'.');
    if length(dotidx) == 1
        H.XTickLabel{hs} = sprintf('%s%c %.0f''',H.XTickLabel{hs}(1:dotidx-1),char(176),str2num(H.XTickLabel{hs}(dotidx:end))*60);
    else
        H.XTickLabel{hs} = sprintf('%s%c 00''',H.XTickLabel{hs},char(176));
    end
end
for hs = 1:length(H.YTickLabel)
    dotidx = strfind(H.YTickLabel{hs},'.');
    if length(dotidx) == 1
        H.YTickLabel{hs} = sprintf('%s%c %.0f''',H.YTickLabel{hs}(1:dotidx-1),char(176),str2num(H.YTickLabel{hs}(dotidx:end))*60);
    else
        H.YTickLabel{hs} = sprintf('%s%c 00''',H.YTickLabel{hs},char(176));
    end
end

%%
% Now Plot major axis, minor axis, profile number and profiles
% on the bathy map
% Reminder, we get data from original map not from any filter 
% 
ridgedistance = 30;
np=1;
p2m = length(lat)/ll2m([lat(1) lat(end)],[long(1) long(1)]);
% divide by 2 because this becomes half the length of the axis
p2l = abs(long(end) - long(1))/length(long)/2;
numprofiles = 0;
for n=1:size(CCs,1)
    if CCs(n,2) < 14
        SD = 93.4;
    else
        SD = 102.9;
    end
    disp(sprintf('Spreading Direction used: %.1f',SD))
   [val C1] = min(abs(long - CCs(n,1)));
   [val C2] = min(abs(lat - CCs(n,2)));
   
    slopes.rps(n).Centroid = [C1 C2];
    
    slopes.rps(n).DistFromAxis = ll2m( [ lat(round(slopeE_Ecc.rps(n).Centroid(2)))...
               axis_interp(round(slopeE_Ecc.rps(n).Centroid(2)),2) ],...
             [ long(round(slopeE_Ecc.rps(n).Centroid(1))) ...
               axis_interp(round(slopeE_Ecc.rps(n).Centroid(2)),1) ])*1e-3 ;


    slopes.rps(n).MajorAxisLength = 150;
    slopes.rps(n).MinorAxisLength = 20;
    slopes.rps(n).Area = 20;

    x = slopes.rps(n).MajorAxisLength*cosd(SD);
    y = slopes.rps(n).MajorAxisLength*sind(SD);
    xs = p2l*x*[-1 1] + long(round(slopes.rps(n).Centroid(1)));
    ys = p2l*y*[1 -1] + lat(round(slopes.rps(n).Centroid(2)));
    slopes.rps(n).MajorAxisTop = [xs(1) ys(1)];
    slopes.rps(n).MajorAxisBot = [xs(2) ys(2)];
    plot(   xs,...
                ys,...
                '-r','LineWidth',2)
        text(-p2l*x + long(round(slopes.rps(n).Centroid(1))),...
             p2l*y + lat(round(slopes.rps(n).Centroid(2))),...
            sprintf('%d',n),'Color','k','FontSize',15)
        x = slopes.rps(n).MinorAxisLength*cosd(SD+90);%slopeW_Wcc.rps(n).Orientation+90);
        y = slopes.rps(n).MinorAxisLength*sind(SD+90);%slopeW_Wcc.rps(n).Orientation+90);
        plot(   p2l*x*[-1 1] + long(round(slopes.rps(n).Centroid(1))),...
                p2l*y*[1 -1] + lat(round(slopes.rps(n).Centroid(2))),...
                '-r','LineWidth',2)         
    [slopes.profile(n).long slopes.profile(n).lat slopes.profile(n).Axis corner1 corner2 ] = Centroid_spreadingprofile(bathy,axis_interp,slopes,lat,long,n,SD);
    plot(slopes.profile(n).long, slopes.profile(n).lat,...
             '--','Color',[1 0 0],'LineWidth',1.5);
        slopes.profile(n).corner1 = [ corner1(1,1) corner1(end,1) corner2(end,1) corner2(1,1) corner1(1,1) ];
        slopes.profile(n).corner2 = [ corner1(1,2) corner1(end,2) corner2(end,2) corner2(1,2) corner1(1,2) ];
        slopes.profile(n).zs = F(slopes.profile(n).long, slopes.profile(n).lat);
        slopes.gravprof(n).mgals = G(slopes.profile(n).long, slopes.profile(n).lat);
end

%%
% intialize
Profslopes = slopes;

%% Now we make the profile structure, 
% this is done manaully with OnePlot_byHand
% or semi automated with other .m files
for n =1:size(CCs,1) 
    axis_to_use = whichaxis(axis_interp,long(round(slopes.rps(n).Centroid(1)) ),lat(round(slopeW_Wcc.rps(n).Centroid(2))),lat); 
    % try disp('Using Already Picked!')
     %  [figi tempprof] = OnePlot_byHand(1,'orthographic',90,Profslopes,slope,n,n,gravtype,lat,long,bathybu,gravy,axis_to_use,Wbf_interp,Ebf_interp,F,G,S,5,0,200);
     %catch
        [figi1, figi2, figi3, figi4, tempprof] = OnePlot_byHand(1,'orthographic',90,slopes,slope,n,n,gravtype,lat,long,bathybu,gravy,axis_to_use,Wbf_interp,Ebf_interp,F,G,S,3,0,200);
        disp('Finished with the Profile!') 
        for js = length(tempprof.profile(n).pick)  
            Profslopes.profile(n).pick(js) = tempprof.profile(n).pick(js);
        end
    % end
    if print_y_n      
        cd proposal' sec4'/
        disp('Saving, please patience')
        disOCC = num2str(CCs(n,2));      
        figi1.PaperPositionMode = 'auto'
        figi2.PaperPositionMode = 'auto'
        figi3.PaperPositionMode = 'auto'
        figi4.PaperPositionMode = 'auto'
%          figi1.PaperUnits = 'inches';
%          figi1.PaperSize = [11 8.5];
%          figi1.PaperPosition = [.0 .0 [11 8.5]-0.5];
        filename1 = sprintf('%d_1_%s_%sN_Nov16_proposal',n,disOCC(1:2),disOCC(4:5)); %%'%s_%sN_%d_DistVSOutRot_Aug17_wShift_wM_wResSqr_wInversion.png',disOCC(1:2),disOCC(4:5),n);
        filename2 = sprintf('%d_2_%s_%sN_Nov16_proposal',n,disOCC(1:2),disOCC(4:5)); %%'%s_%sN_%d_DistVSOutRot_Aug17_wShift_wM_wResSqr_wInversion.png',disOCC(1:2),disOCC(4:5),n);
        filename3 = sprintf('%d_3_%s_%sN_Nov16_proposal',n,disOCC(1:2),disOCC(4:5)); %%'%s_%sN_%d_DistVSOutRot_Aug17_wShift_wM_wResSqr_wInversion.png',disOCC(1:2),disOCC(4:5),n);
        filename4 = sprintf('%d_4_%s_%sN_Nov16_proposal',n,disOCC(1:2),disOCC(4:5)); %%'%s_%sN_%d_DistVSOutRot_Aug17_wShift_wM_wResSqr_wInversion.png',disOCC(1:2),disOCC(4:5),n);
%        filenum = num2str(CCs(n,2));
%        filename = sprintf('%s_%sN_w3Dmisfit_Oct5.png',filenum(1:2),filenum(4:5));
        print(figi1,filename1,'-depsc')
        close(figi1)
        print(figi2,filename2,'-depsc')
        close(figi2)
        print(figi3,filename3,'-depsc')
        close(figi3)
        print(figi4,filename4,'-dpng','-r300')
        close(figi4)
        cd ..
    end
    disp(' *******  Done  ********')
    close all
    pause(1)
end
disp('%%%%%%%%%%%%%%%% All Done %%%%%%%%%%%%%%%%')

%% Now do the inversion

InversionRoutine

disp('More All Done')
