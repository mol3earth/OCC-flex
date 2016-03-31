function Modeled = errorcalcandplot(fetnum,TES,tes,heaves,crusts,infill,angles,offdistX,ProPick,figstring)

%% upack the structure
% and assign values
[ val axidx ] =min(abs(ProPick.profdist));
axDepth = ProPick.zs(axidx);
mHeave = ProPick.Bfcc.Heave*1e-3;
idx1 = ProPick.Ffcc.idx1;
idx2 = ProPick.Ffcc.idx2;
try
    fittingrot = ProPick.fittingrot;% [  ProPick.fittingrot 0 ];
catch
    fittingrot =  ProPick.fittingrot';% [  ProPick.fittingrot' 0 ];
end

try    
    fittingdist = abs(ProPick.fittingdist)';%[   abs(ProPick.fittingdist)' 0 ];
catch        
    fittingdist =abs(ProPick.fittingdist);% [   abs(ProPick.fittingdist) 0 ];
end
    

try                                                         % attempt to make axis depth matter
    fittingtopodist = abs(ProPick.fittingdist)' ; % [ abs(ProPick.fittingdist)' 0 ] ;
catch
    fittingtopodist = abs(ProPick.fittingdist); % [ abs(ProPick.fittingdist) 0 ] ;
end

try
    fittingtopo= ProPick.fittingtopo';%[ ProPick.fittingtopo' axDepth] ;% ProPick.zs(idx1:idx2);
catch
    fittingtopo=  ProPick.fittingtopo;%[ ProPick.fittingtopo axDepth] ;% ProPick.zs(idx1:idx2);
end

fulldist=abs(ProPick.profdist(idx1:idx2)*1e-3);
fulltopo=ProPick.zs(idx1:idx2);


%% keep heave and offset constant
% do for heave
% if only considering measured heave
Acceptable_heaves =  round(mHeave)*1e3; 
% for a small subset of heaves use this
%Acceptable_heaves = round(mHeave)*1e3-500:500:round(mHeave)*1e3+500
for he = 1:length(Acceptable_heaves)
    [v t] = (min(abs(Acceptable_heaves(he) - heaves)));
    accept_he(he) = t;
end
heaves = Acceptable_heaves;
% do for Offset
% measured offset
 Acceptable_offdistX =   round((min(fulldist) - 3.6)/.25)*.25;
% small subset of offsets
%Acceptable_offdistX = Acceptable_offdistX-1:.25:Acceptable_offdistX+1;
for of = 1:length(Acceptable_offdistX)    
    [v t] = (min(abs(Acceptable_offdistX(of) - offdistX)));
    accept_of(of) = t;
end
offdistX = Acceptable_offdistX;
% assign these to the model, cause now they unique
Modeled.offdistX = Acceptable_offdistX;
Modeled.heaves = Acceptable_heaves;
Modeled.crusts = crusts;
Modeled.tes = tes;
Modeled.infill = infill;
Modeled.angles = angles;
%% reset the TES
TES = TES(accept_he,:,:,:);
%% for synthetic
if 0
tesidx = 6
if tesidx == 1
    fittingdist = [TES(12,5,1,1).dist(1:50:500,tesidx)' TES(12,5,1,1).dist(500:200:1500,tesidx)'    TES(12,5,1,1).dist(2000:2000:24000,tesidx)'  TES(12,5,1,1).dist(24500:200:26000,tesidx)'                    TES(12,5,1,1).dist(26000:50:26700,tesidx)' TES(12,5,1,1).dist(26700:200:27500,tesidx)'                         TES(12,5,1,1).dist(28000:1000:end,tesidx)'];
    fittingrot = [TES(12,5,1,1).slope(1:50:500,tesidx)' TES(12,5,1,1).slope(500:200:1500,tesidx)'   TES(12,5,1,1).slope(2000:2000:24000,tesidx)' TES(12,5,1,1).slope(24500:200:26000,tesidx)'                  TES(12,5,1,1).slope(26000:50:26700,tesidx)' TES(12,5,1,1).slope(26700:200:27500,tesidx)'                  TES(12,5,1,1).slope(28000:1000:end,tesidx)'];
elseif tesidx == 6
    fittingtopoP = [TES(12,5,1,1).topo(4542:500:4541+5000,tesidx)'  TES(12,5,1,1).topo(4541+5500:1600:4541+22000,tesidx)' TES(12,5,1,1).topo(4541+23000:500:end,tesidx)' ];

    fittingdistP =  [TES(12,5,1,1).dist(1:500:5000,tesidx)'  TES(12,5,1,1).dist(5500:1600:22000,tesidx)' TES(12,5,1,1).dist(23000:500:end,tesidx)' ];
    fittingrotP =   [TES(12,5,1,1).slope(1:500:5000,tesidx)' TES(12,5,1,1).slope(5500:1600:22000,tesidx)' TES(12,5,1,1).slope(23000:500:end,tesidx)' ];
end   
    fittingdist = fittingdist(1:10);
    fittingrot = fittingrot(1:10) + randn(10,1)'*5;
    fittingtopo = fittingtopo(1:10);
end

%%   Big if
if 1
%% misfit calculation
% Now we have our datas distances, and rotations
% loop through a bunch of offsets, and Te curves, find minimum R^2

if 1
[R_ChiSqX, R_ChiSqTopo] = errorCalc(fittingtopo,fittingtopodist,fittingrot,fittingdist,fulldist,fulltopo,heaves,TES,tes,angles,infill,crusts,offdistX,accept_he,accept_of);
% , R_sqrX, R_sqrY, SumResSqr
end

%%
if 0
    
R_ChiSqX = ProPick.Model.R_ChiSqx;
R_ChiSqTopo = ProPick.Model.R_ChiSqTopo;
end
%% find minimum
%misfitfct = SumResSqr;
%msftstr = 'S/n'; % S/n (Sum Squared Residuals / n)
% % longprofmisfitfct;
% % shortprofmisfitfct;
% % shortprofmisfitfctwNoise;
% % R_ChiSqX = R_ChiSqTopo;
Modeled.R_ChiSqx = R_ChiSqX;
Modeled.R_ChiSqTopo = R_ChiSqTopo;

[bestfitmistfit, idx] = min(R_ChiSqTopo(:)); % max for some misfits
[HEidxT,ANidxT,IFidxT,CTidxT,offdistidxT,teidxT] = ind2sub(size(R_ChiSqTopo), idx);
Modeled.TopoRawminfitidx = [HEidxT,ANidxT,IFidxT,CTidxT,offdistidxT,teidxT];
minfitidxT=[HEidxT,ANidxT,IFidxT,CTidxT,offdistidxT,teidxT];


[bestfitmistfit, idx] = min(R_ChiSqX(:)); % max for some misfits
[HEidx,ANidx,IFidx,CTidx,offdistidx,teidx] = ind2sub(size(R_ChiSqX), idx);
Modeled.SlopeRawminfitidx = [HEidx,ANidx,IFidx,CTidx,offdistidx,teidx];
minfitidx=[HEidx,ANidx,IFidx,CTidx,offdistidx,teidx];


[HEsize ANsize IFsize CTsize ODsize TEsize] = size(R_ChiSqX);

%% Now do the wieght function
[Modeled.maxTw , Modeled.minTw ,Modeled.Weightedminfitidx] = misfitweights(Modeled);
meanTw = (Modeled.maxTw + Modeled.minTw)/2 ; 
minfitidxW = Modeled.Weightedminfitidx;
HEidxW = Modeled.Weightedminfitidx(1);
ANidxW = Modeled.Weightedminfitidx(2);
IFidxW = Modeled.Weightedminfitidx(3);
CTidxW = Modeled.Weightedminfitidx(4);
offdistidxW = Modeled.Weightedminfitidx(5);
teidxW = Modeled.Weightedminfitidx(6);

% keep the ids of the variables that are not constant
if TEsize == 1; minfitidxT(6) = [];minfitidxW(6) = [];minfitidx(6) = [];end
if ODsize == 1; minfitidxT(5) = [];minfitidxW(5) = [];minfitidx(5) = [];end
if CTsize == 1; minfitidxT(4) = [];minfitidxW(4) = [];minfitidx(4) = [];end
if IFsize == 1; minfitidxT(3) = [];minfitidxW(3) = [];minfitidx(3) = [];end
if ANsize == 1; minfitidxT(2) = [];minfitidxW(2) = [];minfitidx(2) = [];end
if HEsize == 1; minfitidxT(1) = [];minfitidxW(1) = [];minfitidx(1) = [];end

%% this is where figi2 also comes from
[numDF ys modelfitsbysigma] = plotErrorSurface(R_ChiSqX,    minfitidx, tes,heaves,angles,infill,crusts,offdistX,HEidx,ANidx,IFidx,CTidx,offdistidx,teidx,sprintf('%g Slope',fetnum),figstring);
%% 
[numDF ys modelfitsbysigmaT] = plotErrorSurface(R_ChiSqTopo,minfitidxT,tes,heaves,angles,infill,crusts,offdistX,HEidxT,ANidxT,IFidxT,CTidxT,offdistidxT,teidxT,sprintf('%g Topo',fetnum),figstring);

%%
[numDF ys modelfitsbysigmaW] = plotErrorSurface((R_ChiSqTopo*meanTw + R_ChiSqX*(1-meanTw)),minfitidxW,tes,heaves,angles,infill,crusts,offdistX,HEidxW,ANidxW,IFidxW,CTidxW,offdistidxW,teidxW,sprintf('%g Weighted by %g',fetnum,meanTw),figstring);
    end

%% find plotting stuff
desirespacing=200;
axisdist = 15;
Nkminc = 1;
% these for the bf dist vs outward rotation
howfar = 3;
BFNkminc = 1;
if isfield(ProPick.Ffcc,'idx1')
    idx1 = ProPick.Ffcc.idx1;
    idx2 = ProPick.Ffcc.idx2;
else
    disp('Pick the left side of the Footwall')
    [idx1 val]=ginput(1);
    [val idx1]=min(abs(ProPick.profdist-idx1));
    disp('Pick the right side of the Footwall')
    [idx2 val]=ginput(1);
    [val idx2]=min(abs(ProPick.profdist-idx2));
    ProPick.Ffcc.idx1 = idx1;
    ProPick.Ffcc.idx2 = idx2;
end
    
% find new spacing for this profile
for hk = 1:abs(idx1-idx2)
    fidminc(hk) = diff([ProPick.profdist(idx1) ProPick.profdist(idx1+hk)]);
end
[val idxhk] = min(abs(fidminc-desirespacing));
Nkminc =  idxhk;
fidminc = fidminc(idxhk);
%fidkminc = diff([ProPick.profdist(idx1) ProPick.profdist(idx1+Nkminc)]);

% subplot(2,2,spp2) 
% CCslopevsTe(ProPick.Bfcc.Heave(fid1)*1e3,abs(ProPick.WBFdist)-3)
% hold on

ProPick.Ffcc.kminc = fidminc; 
% we have to recalculate dzdx, using this new spacing,
% and a resampling as given in input
% multiply by 100 to get toi degrees from slope
onEast = sign(ProPick.profdist(ProPick.bsidx1));
%%
if isfield(ProPick,'Sdzdx')
if onEast
    disp('It'' on the east')
    ProPick.Ffcc.OutwardRotation = ProPick.Sdzdx(idx1:Nkminc:idx2);%gradient(ProPick.zs(idx1:Nkminc:idx2),fidminc)*100;
    OutwardRotations = ProPick.Ffcc.OutwardRotation;
    Distances = ProPick.profdist([idx1:Nkminc:idx2]);
    Depths = ProPick.zs([idx1:Nkminc:idx2]);
    xlimit1 = ProPick.profdist(idx1);
    xlimit2 = 1;
else
    disp('It'' on the west')    
    ProPick.Bfcc.OutwardRotation = ProPick.Sdzdx(idx1:Nkminc:idx2);%gradient(ProPick.zs(idx1:Nkminc:idx2),fidminc)*100;
    OutwardRotations = fliplr(ProPick.Bfcc.OutwardRotation);
    Distances = fliplr(ProPick.profdist([idx1:Nkminc:idx2]));
    Depths = fliplr(ProPick.zs([idx1:Nkminc:idx2])');
    xlimit1 = -1;
    xlimit2 = ProPick.profdist(idx2);
end
ProPick.Distances = Distances;
ProPick.OutwardRotations = OutwardRotations;
%%
else
Distances = ProPick.Distances;
OutwardRotations = ProPick.OutwardRotations*-1;
Depths = fulltopo;
end
countme=1;
PrevRot = 45;
colorN = jet(length(Distances));
PrevDep = -1e100;


fet=fetnum;

printme =1;
doAll = 1;
plotBestfitModel_onSlope;
