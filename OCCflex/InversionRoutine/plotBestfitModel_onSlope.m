% plotBestfitModel_onSlope
fig = figure(fetnum*100000+2);
clf
hold on
%% plot topo

% first assign some constants
minZ = min(Depths);

sp1 = subplot(2,1,1);
hold on
title(sprintf('Feature %g',fet))
VE = 2;
plot(ProPick.profdist,(ProPick.zs- minZ)*VE,'Linewidth',3)
tk2=1;
for tk = 1:length(Distances)
    plot( Distances(tk),...
            (Depths(tk) - minZ)*VE,'o',...
            'MarkerFacecolor',colorN(tk,:),...
            'MarkerEdgecolor',[1 1 1],...
            'MarkerSize',7)
    try
    if round(abs(Distances(tk))) == abs(round(ProPick.fittingdist(tk2)*1e3))
        plot( Distances(tk),...
            (Depths(tk) - minZ)*VE,'o',...
            'MarkerFacecolor',colorN(tk,:),...
            'MarkerEdgecolor',[0 0 0],...
            'MarkerSize',7)
        tk2=tk2+1;
    end
    catch
        
    end
end
axis equal
%%

side = sign(Distances(2))*1e3;
if side > 0
    side2 = 1;
else
    side2 = length(Depths);
end
%%% for posterity
% get a distance vector % dis = TES(HEidx,ANidx,IFidx,CTidx).dist(:,teidx); % tdis = TES(HEidx,ANidx,IFidx,CTidx).topodist(:,teidx);
% get an id for start of footwall
%stID = length(TES(HEidx,ANidx,IFidx,CTidx).topo(:,teidx))-length(TES(HEidx,ANidx,IFidx,CTidx).slope(:,teidx));
%SlopeBFdist = linspace(0,TES(HEidx,ANidx,IFidx,CTidx).dist(end,teidx),length(TES(HEidx,ANidx,IFidx,CTidx).topo)) + offdistX(offdistidx);
%%%

% First plot the weighted model
WeightedModelDistances = side*(TES(HEidxW,ANidxW,IFidxW,CTidxW).topodist(:,teidxW) + offdistX(offdistidxW));

[val idx] = min(abs( WeightedModelDistances - ...
                   side*min(ProPick.fittingdist)  ));
               
%%%%%%%%%%%%%%%%%%%% IF FIXED TO FIRST FOOTWALL POINT %%%%%%%%%%%%%%%%%%%%
% WeightedModelDepths = (TES(HEidxW,ANidxW,IFidxW,CTidxW).topo(:,teidxW) + ...
%     (diff( [ min(TES(HEidxW,ANidxW,IFidxW,CTidxW).topo(:,teidxW)) min(Depths) ] ) - minZ));
%%%%%%%%%%%%%%%%%%%%% IF FIXED TO FIRST FITTED POINT %%%%%%%%%%%%%%%%%%%%%
WeightedModelDepths = TES(HEidxW,ANidxW,IFidxW,CTidxW).topo(:,teidxW) - ...
     TES(HEidxW,ANidxW,IFidxW,CTidxW).topo(idx,teidxW) - minZ - abs(min(ProPick.fittingtopo));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(   WeightedModelDistances,...
        WeightedModelDepths*VE,...
        '-','Color',[.5 .5 .5 .5],'Linewidth',5)
    %%
if doAll
% get a distances vector for the slope model, we must offset it, and
% multiply by the side (-1 if on West, 1 if on East)
SlopeModelDistances = side*(TES(HEidx,ANidx,IFidx,CTidx).dist(:,teidx)  + offdistX(offdistidx));
[val idx] = min(abs( SlopeModelDistances - ...
                   side*min(ProPick.fittingdist)  ));
% find the elevation to offset out Slope model by
% then construct elevation vector
%%%%%%%%%%%%%%%%%%%% IF FIXED TO FIRST FOOTWALL POINT %%%%%%%%%%%%%%%%%%%%
% SlopeModelDepths = (  TES(HEidx,ANidx,IFidx,CTidx).topo(:,teidx) + ...
%     (diff([ min(TES(HEidx,ANidx,IFidx,CTidx).topo(:,teidx)) min(Depths) ]) - minZ));
%%%%%%%%%%%%%%%%%%%%% IF FIXED TO FIRST FITTED POINT %%%%%%%%%%%%%%%%%%%%%
SlopeModelDepths = TES(HEidx,ANidx,IFidx,CTidx).topo(:,teidx) - ...
     TES(HEidx,ANidx,IFidx,CTidx).topo(idx,teidx) - minZ - abs(min(ProPick.fittingtopo));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% do it for Topo model too
TopoModelDistances = side*(TES(HEidxT,ANidxT,IFidxT,CTidxT).topodist(:,teidxT) + offdistX(offdistidxT));
[val idx] = min(abs( TopoModelDistances - ...
                   side*min(ProPick.fittingdist)  ));               
%%%%%%%%%%%%%%%%%%%% IF FIXED TO FIRST FOOTWALL POINT %%%%%%%%%%%%%%%%%%%%
% TopoModelDepths = (TES(HEidxT,ANidxT,IFidxT,CTidxT).topo(:,teidxT) + ...
%     (diff([min(TES(HEidxT,ANidxT,IFidxT,CTidxT).topo(:,teidxT)) min(Depths) ]) - minZ));
%%%%%%%%%%%%%%%%%%%%% IF FIXED TO FIRST FITTED POINT %%%%%%%%%%%%%%%%%%%%%
TopoModelDepths = TES(HEidxT,ANidxT,IFidxT,CTidxT).topo(:,teidxT) - ...
    TES(HEidxT,ANidxT,IFidxT,CTidxT).topo(idx,teidxT) - minZ - abs(min(ProPick.fittingtopo));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(   SlopeModelDistances,...
        SlopeModelDepths*VE,...
        '-','Color',[1 0 1 .7],'Linewidth',2)

plot(   TopoModelDistances,...
        TopoModelDepths*VE,...
        '-','Color',[0 1 0 .7],'Linewidth',1.5)
end
%% plot depths to fault root for models

% find axis
[v Axid]=min(abs(ProPick.profdist));
% find some made up error for axis elevation
FaultRootdepth = abs(ProPick.zs(Axid)) - abs(mean(ProPick.zs(Axid-5:Axid+5))) - 6500;
if doAll
[vT idT] = min(TopoModelDepths);

TopoM_FR =  vT - abs(TopoModelDistances(idT))*tand(angles(ANidxT)) ;

[vS idS] = min(SlopeModelDepths);

SlopeM_FR = vS - abs(SlopeModelDistances(idS))*tand(angles(ANidx)) ;

plot([ SlopeModelDistances(idS) 0],[vS SlopeM_FR]*VE,'--','Color',[1 0 1 .7])
plot([ TopoModelDistances(idT) 0],[vT TopoM_FR]*VE,'--','Color',[0 1 0 .7])
text(-1*side,mean([vT TopoM_FR]*VE),sprintf('Topo Model Depth to Fault Root: %.1f km',abs(TopoM_FR*1e-3)))
if round(mean([vT TopoM_FR]*VE)*1e-3) == round(mean([vS SlopeM_FR]*VE)*1e-3)
    text(-1*side,mean([vS SlopeM_FR]*VE)*.75,sprintf('Slope Model Depth to Fault Root: %.1f km',abs(SlopeM_FR*1e-3)))
else
    text(-1*side,mean([vS SlopeM_FR]*VE),sprintf('Slope Model Depth to Fault Root: %.1f km',abs(SlopeM_FR*1e-3)))
end
end

[vW idW] = min(WeightedModelDepths);

WeightedM_FR = vW - abs(WeightedModelDistances(idW))*tand(angles(ANidxW)) ;

plot([ WeightedModelDistances(idW) 0],[vW WeightedM_FR]*VE,...
        '--','Color',[.5 .5 .5 .5],'Linewidth',4.5)

yval = mean([vW WeightedM_FR]*VE);

if doAll
    if yval == mean([vT TopoM_FR]*VE) | yval == mean([vS SlopeM_FR]*VE)
        yval = yval/2;
    end
end
    

text(-1*side,yval,sprintf('Weighted Model Depth to Fault Root: %.1f km',abs(WeightedM_FR*1e-3)))

% plot axis 
plot([ 0 0 ], [ -1e4 2e3],'--k','Linewidth',1.5 )

%% finally, set xlims & ylims

if doAll
[v MSid]=max(abs(TES(HEidx,ANidx,IFidx,CTidx).topo(:,teidx)));
[v MTid]=max(abs(TES(HEidxT,ANidxT,IFidxT,CTidxT).topo(:,teidxT)));

xlim([  min([ min(ProPick.profdist) TopoModelDistances(MTid) SlopeModelDistances(MSid)  ]) ...
        max([max(ProPick.profdist) TopoModelDistances(MTid) SlopeModelDistances(MSid) ]) ])

ylim(VE*[ min([min((Depths(tk) - minZ )) min(TopoM_FR) min(SlopeM_FR)]) ...
    max([max((ProPick.zs - minZ )) max(TopoModelDepths) max(SlopeModelDepths)]) ...
        ])
else
    
[v MWid]=max(abs(TES(HEidxW,ANidxW,IFidxW,CTidxW).topodist(:,teidxW)));

xlim([  min([ min(ProPick.profdist) WeightedModelDistances(MWid)  ]) ...
        max([max(ProPick.profdist) WeightedModelDistances(MWid)  ]) ])

ylim(VE*[ min([min((Depths(tk) - minZ )) min(WeightedM_FR) ]) ...
    max([max((ProPick.zs - minZ )) max(WeightedModelDepths) ]) ...
        ])
end
ylabel(sprintf('Pseudodepths (normalized, and exagerated) (m*%.0f)',VE))
xlabel('Distance from Axis (m)')
%% plot slope
subplot(2,1,2) 
% offdist = 0;
%    [te(1).dist, te(1).slope] =
%     CCslopevsTe(60000, offdist,'-k',12,tes);
hold on
colorN = jet(length(Distances));
tk2=1;
for tk = 1:length(Distances)
    plot(  abs(Distances(tk)*1e-3),... 
        -(OutwardRotations(tk)),... 
        'o',    'MarkerFacecolor',colorN(tk,:),... colors(fid,:),...
                'MarkerEdgecolor',[1 1 1],...
                'MarkerSize',8)
    try        
    if abs(round(Distances(tk))) == abs(round(ProPick.fittingdist(tk2)*1e3))
        plot(  abs(Distances(tk)*1e-3),... 
            -(OutwardRotations(tk)),... 
            'o',    'MarkerFacecolor',colorN(tk,:),... colors(fid,:),...
            'MarkerEdgecolor',[0 0 0 ],...
            'MarkerSize',8)
        tk2=tk2+1;
    end
    end
end

% plot weighted model slope
plot(   TES(HEidxW,ANidxW,IFidxW,CTidxW).dist(:,teidxW)+offdistX(offdistidxW),...
        TES(HEidxW,ANidxW,IFidxW,CTidxW).slope(:,teidxW),...
        '-','Color',[.5 .5 .5 .5],'Linewidth',5)
    
if doAll
% plot slope model slope
plot(   TES(HEidx,ANidx,IFidx,CTidx).dist(:,teidx)+offdistX(offdistidx),...
        TES(HEidx,ANidx,IFidx,CTidx).slope(:,teidx),...
        '-','Color',[1 0 1 .7],'Linewidth',1.5)

% plot topo model slope
% get a distance vector
TopoBFdist = linspace(0,TES(HEidxT,ANidxT,IFidxT,CTidxT).dist(end,teidxT),length(TES(HEidxT,ANidxT,IFidxT,CTidxT).topo)) + offdistX(offdistidxT);
plot(   TES(HEidxT,ANidxT,IFidxT,CTidxT).topodist(:,teidxT)+offdistX(offdistidxT),...TopoBFdist(stID+1:end),...
        TES(HEidxT,ANidxT,IFidxT,CTidxT).slope(:,teidxT),...
        '-','Color',[0 1 0 .7],'Linewidth',1.5)

text(11,-10,sprintf('Slope Fit (magenta):\nHeave = %g km \nAngle = %g^o \nInfill = %g km\nCrust = %g km\nOffset = %g km\nTe = %g m',...
    heaves(HEidx)*1e-3,...
    angles(ANidx),...
    infill(IFidx)*1e-3,...
    crusts(CTidx)*1e-3,...
    offdistX(offdistidx),...
    tes(teidx)),...
    'Fontsize',12 ...
    )
text(11,-30,sprintf('Topo Fit (green):\nHeave = %g km \nAngle = %g^o \nInfill = %g km\nCrust = %g km\nOffset = %g km\nTe = %g m',...
    heaves(HEidxT)*1e-3,...
    angles(ANidxT),...
    infill(IFidxT)*1e-3,...
    crusts(CTidxT)*1e-3,...
    offdistX(offdistidxT),...
    tes(teidxT)),...
    'Fontsize',12 ...
    )
end
   
    
text(15,-20,sprintf('Weighted Slope Fit (grey):\nHeave = %g km \nAngle = %g^o \nInfill = %g km\nCrust = %g km\nOffset = %g km\nTe = %g m',...
    heaves(HEidxW)*1e-3,...
    angles(ANidxW),...
    infill(IFidxW)*1e-3,...
    crusts(CTidxW)*1e-3,...
    offdistX(offdistidxW),...
    tes(teidxW)),...
    'Fontsize',15, ...
    'Fontweight','bold'...
    )

text(1,-30,sprintf('Feature: %s\nMeasured:\nHeave = %.1f km\nWeight = %.2f ',...
    fet,...
    mHeave,...
    meanTw),...ProPick.Bfcc.Heave*1e-3),...
    'Fontsize',12 ...
    )
xlim([0 20])
ylim([-41 10])
xlabel('Distance from Axis (km)')
ylabel('Slope (Degrees)')
if printme
cd mfiles
save2pdf(sprintf('%s_Feature_%s_Bestfit_Slopes',figstring,fetnum),fig)
cd ..
end