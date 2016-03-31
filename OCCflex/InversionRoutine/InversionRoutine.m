%% set up parameter vectors, and calculate our model structure
%heaves = [ProPick.Bfcc.Heave:5000:60000 ];
tes =[ 100:100:2000 ]%[10,50,100,200,250,500,750,1000,1250,1500,1750,2000];% [10,25,50,100,200,250,500,750,1000,1500,2000];

for heavers = 1:length(Profslopes.profile)
    heaves(heavers) = Profslopes.profile(heavers).pick.Bfcc.Heave;
end

% heaves = [1000:500:12000] ;%[1000:500:12000];%[ [1000:1000:10000] [15000:5000:30000]];
crusts =[0 1e3 3e3 6e3]; [0:3:10];
infill = crusts; crusts;[0:1:6]*1e3;
angles = [45:5:75];
% <45 is impossible because 
% critical_angle = ( 90 + internal_friction_angle )/ 2
% >75 impossible because


% set up offset vector
offdistX =  -3:.25:3;

if 1
TEmodels = TESmaker(heaves,angles,infill,crusts,6,tes)
end
%%
figstring = 'Mar30_absheave_fixoff_vIF_firstfitpt';

for n = [ 1:length(Profslopes.profile)]
    ErrorStructwCrustThick = errorcalcandplot(n,TEmodels,tes,heaves,crusts,infill,angles,offdistX,...
            Profslopes.profile(n).pick,figstring);
    Profslopes.profile(n).pick.ModelwCTMar30 = ErrorStructwCrustThick;
end

%% for Dantes, Kane and Tag.

DANTESpick.figstring = 'Dantes_pick_Mar30_fixheave_fixoff_firstfitpt';
KANEpick.figstring = 'Kane_pick_Mar30_fixheave_fixoff_firstfitpt';
TAGpick.figstring = 'TAG_pick_Mar30_fixheave_fixoff_firstfitpt';
time = 1;
for pickpick = [ DANTESpick KANEpick TAGpick]
    ErrorStructwCrustThick = errorcalcandplot(n+time,TEmodels,tes,heaves,crusts,infill,angles,offdistX,...
        pickpick,pickpick.figstring);
    pickpick.ModelwCT = ErrorStructwCrustThick;
    time = time+1;
end


