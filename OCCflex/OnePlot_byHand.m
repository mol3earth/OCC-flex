function [fighandle figi2 figi3 figi4 newslopestruct] = OnePlot_byHand(LA,projection,downangle,slopestruct,slopemap,n,profilenumber,gravtype,lat,long,bathy,gravy,axis_interp,WBF_interp,EBF_interp,F,G,S,VE,BForFF,desirespacing)
       %% first, find region
            yeslat=find(lat<max(slopestruct.profile(n).corner2) & lat>min(slopestruct.profile(n).corner2));
            yeslong=find(long<max(slopestruct.profile(n).corner1) & long>min(slopestruct.profile(n).corner1));
            tlat = lat(yeslat);
            tlong = long(yeslong);
            tBathy = bathy(yeslat,yeslong);
            tgravy = gravy(yeslat,yeslong);
            tslope = slopemap(yeslat,yeslong);
            isfiniteaxis1 = axis_interp(find( isfinite(axis_interp(:,1))),1) ;
            isfiniteaxis2 = axis_interp(find( isfinite(axis_interp(:,2))),2) ;
            isnotfiniteaxis1 = axis_interp(find( isfinite(axis_interp(:,1))==0 ),1);
            isnotfiniteaxis2 = axis_interp(find( isfinite(axis_interp(:,2))==0 ),2);
            isnotfiniteaxis1 = isfiniteaxis1(end);
            isnotfiniteaxis2 = isfiniteaxis2(end);
            axis_interp(find( isfinite(axis_interp(:,2))==0 ),2) = isnotfiniteaxis2;
            axis_interp(find( isfinite(axis_interp(:,1))==0 ),1) = isnotfiniteaxis1;            
        % yay, found region
        % Now find colors for profiles
        % grab the middlemost or whatever
        colorss = autumn(1);
        gravcolorss = winter(1);
        [val idxs]=max(tBathy);
        [val idx1]=max(val);
        indexms = idxs(idx1);
        indexms = floor(length(tlat)/2);
        %indexms = 41;  %Aug10 13_20N
        %indexms = 65; %Aug10 13_47N
        %indexms = 115;%89; % 13.1N
        
        disp(sprintf('%d: %d %d %d , and %d',n,indexms,length(tlat)))
        sp1 = 2;
        sp2 = 2;
        spp1 = 1;
        spp2 = 3;

         %%  First, go grab the pertinent information    
       if isfield(slopestruct.profile(n),'pick') == 0
            for hk = 2:length(slopestruct.profile(n).lat)
                fidminc(hk) = ll2m([slopestruct.profile(n).lat(1) slopestruct.profile(n).lat(hk)],[slopestruct.profile(n).long(1) slopestruct.profile(n).long(hk)]);
            end
            [val idxhk] = min(abs(fidminc-desirespacing)); 
            % note use of idxhk % 
            slopestruct.profile(n).pick.long = slopestruct.profile(n).long(1:idxhk:end);
            % note use of idxhk % 
            slopestruct.profile(n).pick.lats = slopestruct.profile(n).lat(1:idxhk:end);
            % note use of idxhk % 
            % assign bathymetry values to this profile
            slopestruct.profile(n).pick.zs = slopestruct.profile(n).zs(1:idxhk:end);
            % note use of idxhk % 
            % assign gravity values to this profile
            slopestruct.profile(n).pick.gs = G(...
                slopestruct.profile(n).long(1:idxhk:end) ,...
                slopestruct.profile(n).lat(1:idxhk:end));
            % assign axis
            slopestruct.profile(n).pick.Ax = slopestruct.profile(n).Axis;
            % assign bounding fault locations
            % these can be overwritten later
            slopestruct.profile(n).pick.EBF = [EBF_interp(yeslat(indexms),1) EBF_interp(yeslat(indexms),2)];    
            slopestruct.profile(n).pick.WBF = [WBF_interp(yeslat(indexms),1) WBF_interp(yeslat(indexms),2)];
            % assign slope values
            slopestruct.profile(n).pick.dzdx = S( ...
                slopestruct.profile(n).long(1:idxhk:end) ,...
                slopestruct.profile(n).lat(1:idxhk:end));
            % find entire profile length
            profdist = ll2m([slopestruct.profile(n).pick.lats(1) slopestruct.profile(n).pick.lats(end)],[slopestruct.profile(n).pick.long(1) slopestruct.profile(n).pick.long(end)]);
            % find distance from LL point to axis
            axdist = ll2m([slopestruct.profile(n).pick.lats(1) slopestruct.profile(n).pick.Ax(2)],[slopestruct.profile(n).pick.long(1) slopestruct.profile(n).pick.Ax(1)]);
            slopestruct.profile(n).pick.axdist = axdist;
            % make a vector for our profile, and offset to make origin the axis
            profdist = linspace(0,profdist,length(slopestruct.profile(n).pick.zs))  - axdist;
            slopestruct.profile(n).pick.profdist = profdist;
            [val Axidx]=min(abs(slopestruct.profile(n).pick.profdist));
            slopestruct.profile(n).pick.dzdx(Axidx) = NaN;
            % find sign of the slopes relative to facing toward or away
            % from the axis
            try slopestruct.profile(n).pick.Sdzdx = slopestruct.profile(n).pick.dzdx'.*sign(slopestruct.profile(n).pick.profdist);
            catch slopestruct.profile(n).pick.Sdzdx = slopestruct.profile(n).pick.dzdx.*sign(slopestruct.profile(n).pick.profdist);
            end
            % plot it up, and grab some user picks
            figure(10010)
            clf
            % first plot the profile, and then replot after finding new
            % axis
            plot(slopestruct.profile(n).pick.profdist,slopestruct.profile(n).pick.zs)
            disp(' Pick the axis ')                
            [ axdist Zval ] = ginput(1);
            slopestruct.profile(n).pick.profdist = profdist - axdist;
            clf
            plot(slopestruct.profile(n).pick.profdist,slopestruct.profile(n).pick.zs)
                disp(' <----   Pick the left side of the back slope')
                [ leftdist Zval ] = ginput(1);
                [val idx1]=min(abs(slopestruct.profile(n).pick.profdist-leftdist));
                disp('    Pick the right side of the back slope  ---->')
                [ rightdist Zval ] = ginput(1) ;
                [val idx2]=min(abs(slopestruct.profile(n).pick.profdist-rightdist));
                disp(' Choose the Bounding Fault, or Toe of the OCC')
                [ slopestruct.profile(n).pick.EBFdist val] = ginput(1);
                % only need one BF for this case
                slopestruct.profile(n).pick.WBFdist = slopestruct.profile(n).pick.EBFdist;
% %                 disp('Choose the WBF')
% %                 [slopestruct.profile(n).pick.WBFdist val] = ginput(1);
                % now ask if want to confine the plot to this side only
                if input('Press 1 if want to confine profile to OCC side only:  ') == 1
                if idx1 > Axidx
                    slopestruct.profile(n).pick.profdist = slopestruct.profile(n).pick.profdist(Axidx:end);
                    slopestruct.profile(n).pick.Sdzdx = slopestruct.profile(n).pick.Sdzdx(Axidx:end);
                    slopestruct.profile(n).pick.dzdx = slopestruct.profile(n).pick.dzdx(Axidx:end);
                    slopestruct.profile(n).pick.long = slopestruct.profile(n).pick.long(Axidx:end);
                    slopestruct.profile(n).pick.lats = slopestruct.profile(n).pick.lats(Axidx:end);
                    slopestruct.profile(n).pick.zs = slopestruct.profile(n).pick.zs(Axidx:end);
                    slopestruct.profile(n).pick.gs = slopestruct.profile(n).pick.gs(Axidx:end);
                    idx1 = idx1-Axidx;
                    idx2 = idx2-Axidx;
                    Axidx = 1;
                else
                    slopestruct.profile(n).pick.profdist = slopestruct.profile(n).pick.profdist(1:Axidx);
                    slopestruct.profile(n).pick.Sdzdx = slopestruct.profile(n).pick.Sdzdx(1:Axidx);
                    slopestruct.profile(n).pick.dzdx = slopestruct.profile(n).pick.dzdx(1:Axidx);
                    slopestruct.profile(n).pick.long = slopestruct.profile(n).pick.long(1:Axidx);
                    slopestruct.profile(n).pick.lats = slopestruct.profile(n).pick.lats(1:Axidx);
                    slopestruct.profile(n).pick.zs = slopestruct.profile(n).pick.zs(1:Axidx);
                    slopestruct.profile(n).pick.gs = slopestruct.profile(n).pick.gs(1:Axidx);
                end                
                % set the ids
                end
                
                slopestruct.profile(n).pick.bsidx1=idx1;
                slopestruct.profile(n).pick.bsidx2=idx2;
                    
       else
           disp('Using old shit')
            idx1 = slopestruct.profile(n).pick.bsidx1;
            idx2 = slopestruct.profile(n).pick.bsidx2;
       end
      
%%         % and plot it
        % find colormaps for our back facing slopes
        cm.colors = lines(1);
         %%   
         
         % go edit figure and subplot
    fighandle = figure(n*100000+1);
    clf
%%        subplot(sp1,sp2,spp1)
 %       subplot(2,1,1)
        hold on
            olzs = slopestruct.profile(n).pick.zs;
            olgs = slopestruct.profile(n).pick.gs;
            olss = slopestruct.profile(n).pick.Sdzdx;
            thisbathy = slopestruct.profile(n).pick.zs;
            mingravy = min(olgs);
            minslope = min(olss);
            thisgravy = olgs - mingravy ;
            % do the stuff for gravyity 
            [maxgravy gidx] = max(abs(thisgravy));
            maxgravy = maxgravy*sign(thisgravy(gidx));
            thisgravy = thisgravy*(abs(diff([max(olzs) min(olzs)]))/maxgravy);
            thisgravy = thisgravy + min(olzs);
            % do the same for slope profile            
            thisslope = olss - minslope ;
            gridspaces = [round(min(olss),-1):20:round(max(olss),-1)];
            gridslope = gridspaces -minslope;
            zeroslope = -minslope;
            [maxslope sidx] = max(abs(thisslope));
            maxslope = maxslope*sign(thisslope(sidx));
            thisslope = thisslope*(abs(diff([max(olzs) min(olzs)]))/maxslope);
            gridslope = gridslope*(abs(diff([max(olzs) min(olzs)]))/maxslope);
            zeroslope = zeroslope*(abs(diff([max(olzs) min(olzs)]))/maxslope);
            thisslope = thisslope + min(olzs);
            gridslope = gridslope + min(olzs);
            zeroslope = zeroslope + min(olzs);
            % plot a zero line for the slope
% % %             plot([ones(1,length(gridslope))*slopestruct.profile(n).pick.profdist(1) ;ones(1,length(gridslope))*slopestruct.profile(n).pick.profdist(end)],...
% % %                         [gridslope ; gridslope]*VE,'--','Color',[0 .7 0])
% % %             for mnm = 1:length(gridslope)
% % %                text(slopestruct.profile(n).pick.profdist(end),...
% % %                               gridslope(mnm)*VE,sprintf('  %s%c',num2str(gridspaces(mnm)),char(176)),'Fontname','Cambria','FontSize',15)
% % %                     
% % %             end
% % %             plot([slopestruct.profile(n).pick.profdist(1) slopestruct.profile(n).pick.profdist(end)],...
% % %                         [1 1]*zeroslope*VE,'-','color',[0 .5 0])
% % %             text((slopestruct.profile(n).pick.profdist(end)),...
% % %                         zeroslope*VE,sprintf('   0%c',char(176)),'Fontname','Cambria','FontSize',15)
% % %                     
% % %             plot(slopestruct.profile(n).pick.profdist,thisslope*VE,'--g','LineWidth',1.5)
           % plot(slopestruct.profile(n).pick.profdist*1e3/VE,thisgravy,'--c','LineWidth',1.5)
            plot(slopestruct.profile(n).pick.profdist,thisbathy*VE,'Color',colorss,'LineWidth',1.5)
            plot(slopestruct.profile(n).pick.WBFdist*[1 1],[min(thisbathy*VE) mean(thisbathy*VE)],'b-')
            % plot OCC
            slopestruct.profile(n).pick.Bfcc.OutwardRotation = [];
            slopestruct.profile(n).pick.Bfcc.DistFromAxis = [];
            slopestruct.profile(n).pick.Bfcc.Heave = [];
            
            [ slopestruct.profile(n).pick Ffcc] = bfccstuff(slopestruct.profile(n).pick,0,VE,cm.colors,idx1,idx2);
        axis tight
        axis equal
        fname = 'Cambria';

        h=get(gca);
        hold on
        plot([0 0],...
                h.YLim,...
                '--','Color',[.2 .2 .2])
        % print axis text
 %       text(sign(slopestruct.profile(n).pick.profdist(2))*500,mean(h.YLim), 'Axis','rotation',90,'Fontname','Cambria','Fontsize',17)
    %%        text(h.XLim(1)*.85,h.YLim(1),sprintf('VE: %g',VE),'Fontname',fname,'Fontsize',12)
        text(mean(h.XLim)-...
                sign(h.XLim(find(max(abs(h.XLim)))))*mean(h.XLim)*.7...
                ,h.YLim(1)*.95,sprintf('VE: %g',VE),'Fontname',fname,'Fontsize',17)

        title(sprintf('Profile #: %d',profilenumber))

%%%        ylabel('''Relative'' Depth ','Fontsize',12,'Fontname',fname)
%%%        xlabel('Distance from Axis (km)','Fontsize',12,'Fontname',fname)
        text(slopestruct.profile(n).pick.profdist(1)*.97+slopestruct.profile(n).pick.profdist(end)*.03,...
                h.YLim(1)*.9,' Depth (km)','rotation',90,'Fontname',fname,'Fontsize',18)
        text(mean(h.XLim),h.YLim(1)*.95,'Distance from Axis (km)','Fontsize',18,'Fontname',fname)

%        set(get(AX(2),'Ylabel'),'String',sprintf('%s (mGal)',gravtype),'Fontsize',10,'Fontname',fname);
        % not sure why, but this weird method is the only way I could get
        % the tick label to change
        set(gca,'Fontsize',18,'Fontname',fname)
        h=get(gca);
        for hnum = 1:length(h.XTickLabel)
            h.XTickLabel{hnum} = num2str(h.XTick(hnum)*1e-3);
        end
        set(gca,'XTickLabel',h.XTickLabel')
  
       %%
       
        SPI.sp1 = sp1;
        SPI.sp2 = sp2;
        SPI.spp1 = spp1;
        SPI.spp2 = spp2;
        if BForFF 
          %  figure(n*100000+2)
     %%      subplot(2,2,spp2) 
           plotrotationvsTe(1,0,0)
           hold on
           for m = 1:min(length(indexms),length(slopestruct.profile(n).pick))
            if isfield(slopestruct.profile(n).pick.Bfcc,'LargeObjects')
                TePlot(slopestruct.profile(n).pick,m,cm.colors(1,:))
            end
           end 
        else
           for m =  1:min(length(indexms),length(slopestruct.profile(n).pick))
               if isfield(slopestruct.profile(n).pick,'Ffcc')
                   display('Using old Ffcc')
               else
                    slopestruct.profile(n).pick.Ffcc = Ffcc;
               end
                display('**************');display('*******************');display('************************')
                display('**************');display('*******************');display('************************')
                display('**************');display('*******************');display('************************')
                onEast = 0;%input('Is OCC on East (1) or West (0)?     ');
                [slopestruct.profile(n).pick , onEast , figi2 ]= inversion_plotupFF(onEast,slopestruct.profile(n).pick,1,colorss(m,:),desirespacing,SPI,VE,0,n);                                                                              %cm.colors(1,:)            
           end
        end
       if onEast == -1
           onEast=0;
       end
       figi3 = figi2;
       newslopestruct = slopestruct;
         %% plot mapview
figi4=figure(n*100000+4)
%    figure(n*100000+1);
    clf
%    subplot(2,1,2)
%%         subplot(2,2,2)
        s=surf(tlong,tlat, tBathy);
%        l=light('Position',[axis_interp(yeslat(end-5),1)  axis_interp(yeslat(end-5),2) 1e-1],'Style','infinite');
        lightangle(onEast*180+93.4,.0001)
        zlim([min(min(tBathy)) max(max(tBathy))])
        xlim([min(tlong) max(tlong)])
        ylim([min(tlat) max(tlat)])
        s.FaceLighting = 'flat';
        s.AmbientStrength = 0.7;%.4
        s.DiffuseStrength = .9;%.5
        s.SpecularStrength = .8;%.6
        s.SpecularExponent = .9;%.5
        s.BackFaceLighting = 'reverselit'; %'reverselit' | 'unlit' | 'lit'
        %view( 93.4-90,...
        %        downangle )
        box on
        view(0, 90)
        shading interp
        axis equal
        colormap(jet)
        hold on         
        [c h]=contour3(tlong,tlat, tBathy,-7000:100:-1000,'k');
%          fill3([ slopestruct.profile(n).long(1) slopestruct.profile(n).long(end) slopestruct.profile(n).long(end) slopestruct.profile(n).long(1) slopestruct.profile(n).long(1)],...
%                [ slopestruct.profile(n).lat(1)  slopestruct.profile(n).lat(end) slopestruct.profile(n).lat(end) slopestruct.profile(n).lat(1) slopestruct.profile(n).lat(1)],...
%                [ max(max(tBathy)) max(max(tBathy)) min(min(tBathy)) min(min(tBathy)) max(max(tBathy)) ],...
%                  [1 0 0 ],'FaceAlpha',.4);
        colid = 1;
        plot3(  slopestruct.profile(n).pick.long,...
                slopestruct.profile(n).pick.lats,...
                    slopestruct.profile(n).pick.zs,...
                    '-','Color',[1 1 1],... colorss(colid,:),...
                    'LineWidth',2);
        plot3(axis_interp(yeslat,1),axis_interp(yeslat,2),...
              F(axis_interp(yeslat,1),axis_interp(yeslat,2))*1.01,...
             '-','Color',[.2 .2 .2],'LineWidth',2.5 )
        set(gca,'Visible','on','Projection',projection);
        set(gca,'GridLineStyle','none')
        set(gca,'Ylim',[min(tlat) max(tlat)])        
        set(gca,'Xlim',[min(tlong) max(tlong)])
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
function [ProPick Ffcc2] = bfccstuff(ProPick,offset,VE,colors,idx1,idx2)
%%
yolo=1;
RLid = 1; 

% find the ids for these
[val Axidx]=min(abs(ProPick.profdist));
    
    % find max rotation/slope of the Back facer 
    PPs = ProPick.dzdx(idx1:idx2);
    PPs = sort(PPs);
    PPs = PPs(1:round(length(PPs)/1.65));
    ProPick.Bfcc.OutwardRotation = max(abs(PPs));%abs(mean(PPs));
    % if it is large enough, keep it
    if sign(ProPick.profdist(round((idx1+idx2)/2))) == -1
        % on west
        pidx1 = idx2;
        [val pidx2]=min(abs(ProPick.profdist - ProPick.WBFdist));
        ProPick.Bfcc.DistFromAxis = abs(ProPick.profdist(idx2));
        ProPick.Bfcc.Heave = abs(diff([ProPick.profdist(idx2) ProPick.WBFdist]));
    else
        % on east %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [val pidx1]=min(abs(ProPick.profdist - ProPick.EBFdist));
        pidx2 = idx1;
        ProPick.Bfcc.DistFromAxis = abs(ProPick.profdist(idx1));
        ProPick.Bfcc.Heave = abs(diff([ProPick.profdist(idx1) ProPick.EBFdist]));
    end
    % plot the profile, and the slope angle
    Ffcc2.PixelIdxList = sort([pidx2:1:pidx1]);
    plot(ProPick.profdist(idx1:idx2),...
        (ProPick.zs(idx1:idx2)+offset)*VE,...
        '-','Color',colors(yolo,:),'LineWidth',2);
    text(   ProPick.profdist(idx1),...
            (mean(ProPick.zs(idx1:idx2))+offset)*VE,...
            sprintf('%.1f^o', ProPick.Bfcc.OutwardRotation),...
            'Fontsize',15,'Fontname','Cambria')
    % plot the offset bar, and the amount
    plot([ProPick.profdist(pidx1) ProPick.profdist(pidx2)],...
            ((max(ProPick.zs(idx1:idx2))+offset)*[.99 .99])*VE,...
            '-','Color',colors(yolo,:),'LineWidth',1.5);
    text(  ProPick.profdist(pidx1),...mean(ProPick.profdist(pidx1:pidx2))
               (( -100+  max(ProPick.zs(idx1:idx2)))+offset)*VE,...
                sprintf('Heave %.1f km',ProPick.Bfcc.Heave*1e-3),...
                'Fontsize',15,'Fontname','Cambria')
    yolo=yolo+1;

function TePlot(ProPick,m,colors)
%% plot slope vs offset for back face

yolo=1;        
% instead of the mean do a slope from higheset to lowest value. 
% maybe try to use different spacings, and then grab the highest
% value? 
% try averaging across 500 meters, whatever that is in pixels
% that is, taking the slope every n number of pixels instead of
% taking mean

for bcid = ProPick.Bfcc.LargeObjects
    plot(  abs( ProPick.Bfcc.Heave(bcid)),...
            abs(ProPick.Bfcc.OutwardRotation(bcid)),...
            'o',    'MarkerFaceColor',colors(yolo,:),...
                    'MarkerEdgecolor',[ 1 0 0 ])
    text(   max(abs(ProPick.Bfcc.Heave(bcid))+...
                    1),...
            abs(ProPick.Bfcc.OutwardRotation(bcid)),...
            sprintf('%d.%d',m,yolo),'Color',[.75 0 0])
    yolo=yolo+1;
end  
xlim([0 20])
ylim([0 47])
ylabel('Outward Rotation (Degrees)')
xlabel('Offset/Heave (km)') 
set(gca,'Fontsize',12,'Fontname','Ubuntu')
disp('subplot 2 fin')

