function [maxTw , minTw , Weightedminfitidx] = misfitweights(Misfits)

% first, find the misfits which have a value. Some models were not run, and
% have a 1e6 misfit.
TopoMF = (Misfits.R_ChiSqTopo(find( Misfits.R_ChiSqTopo(:) ~= max(Misfits.R_ChiSqTopo(:)))));
SlopMF = (Misfits.R_ChiSqx(find( Misfits.R_ChiSqx(:) ~= max(Misfits.R_ChiSqx(:)))));

% ScaledSlope = Misfits.R_ChiSqx/var(SlopMF);
% ScaledTopo = Misfits.R_ChiSqTopo/var(TopoMF);


testweights = [.0001:.0001:.9999];

idxs = zeros(length(testweights),1);
bestfitmistfits = idxs;
for m = 1:length(testweights)
    [bestfitmistfits(m) , idxs(m) ] = min(Misfits.R_ChiSqx(:)*testweights(m) ...
                                          + Misfits.R_ChiSqTopo(:)*(1-testweights(m)));            
end
%%
clear distances
distances = sqrt( ([Misfits.R_ChiSqx(:).^ 2 + Misfits.R_ChiSqTopo(:).^ 2]) );
[mindist distID ] =min(distances);
figure(2);clf
subplot(3,1,3)
plot(Misfits.R_ChiSqx(:),Misfits.R_ChiSqTopo(:)','*')
hold on
plot(Misfits.R_ChiSqx(distID),Misfits.R_ChiSqTopo(distID)','ro')
ylabel('Topography \chi^2')
xlabel('Slope \chi^2')
title('Misfit space: Slope \chi^2 vs Topography \chi^2 space')
%%
subplot(3,1,1)
h = histogram(TopoMF);
hold on
histogram(SlopMF,h.BinEdges)
xlabel('Misfit value')
ylabel('Frequency')
title('Histogram of the misfits in each model space')
legend('Topography Misfits','Slope Misfits')
%%

% this always gives Tw = max and Sw = min
% idx = find(bestfitmistfits == min(bestfitmistfits));
% 
% Sw = testweights(idx);
% Tw = 1 - Tw;
% 
% [HEidx,ANidx,IFidx,CTidx,offdistidx,teidx] = ind2sub(size(Misfits.R_ChiSqx), idxs(idx));
% Weightedminfitidx = [HEidx,ANidx,IFidx,CTidx,offdistidx,teidx];

% lets try to find the most numerous bestfit
uids = unique(idxs);

for n = 1:length(uids)
    numofmodels(n) = numel(find(idxs == uids(n)));
end

subplot(3,1,2)
bar((uids),log(numofmodels))
title('Histogram of best fit models')
xlabel('Matrix ID')
ylabel('Log_1_0 Frequency')
for n = 1:length(uids)
    text(uids(n)+100,log(numofmodels(n)),sprintf('%g',numofmodels(n)))
end

bfidx = find(bestfitmistfits == min(bestfitmistfits));
Tw = 1 - testweights(bfidx);
text(idxs(bfidx),4,sprintf('Minimun, with Topo_{weight}: %g',Tw))

%%
acceptable = find(idxs == uids(find(numofmodels == max(numofmodels))));
acceptableID = uids(find(numofmodels == max(numofmodels)));
maxTw = min(testweights(acceptable));
minTw = max(testweights(acceptable));

[HEidx,ANidx,IFidx,CTidx,offdistidx,teidx] = ind2sub(size(Misfits.R_ChiSqx), acceptableID);
if acceptableID ~= distID
    disp('Most common and Shortest Distance misfits are NOT EQUAL')
    disp('Going with distance model')
    [HEidx,ANidx,IFidx,CTidx,offdistidx,teidx] = ind2sub(size(Misfits.R_ChiSqx), distID);
end

Weightedminfitidx = [HEidx,ANidx,IFidx,CTidx,offdistidx,teidx];

text(acceptableID+100,8,sprintf('Tw_m_e_a_n = %g \nTw_r_a_n_g_e = %g - %g',(maxTw+minTw)/2,maxTw,minTw))
