function TES = TESmaker(heaves,angles,infill,crusts,depths,tes)

% initialize big mat
% % % BigMat.Slope = zeros(29999,...
% % %                          (length(angles)* ...
% % %                                 length(infill)* ...
% % %                                 length(crusts)* ...
% % %                                 length(heaves)* ...
% % %                                 length(tes)* ...
% % %                                 length(offdistX)) ...
% % %                                 );
% % % BigMat.Dist = BigMat.Slope;
% % % BigMat.Topo = zeros(34540,...
% % %                          (length(angles)* ...
% % %                                 length(infill)* ...
% % %                                 length(crusts)* ...
% % %                                 length(heaves)* ...
% % %                                 length(tes)* ...
% % %                                 length(offdistX)) ...
% % %                                 );
% % initialkize models structure
% TES(1:length(heaves),1:length(angles),1:length(infill),1:length(crusts)).topo =0;
% TES(1:length(heaves),1:length(angles),1:length(infill),1:length(crusts)).slope=0;


tic
BMcount = 1;
    for he = 1:length(heaves)
        he
        %[te(he).dist, te(he).slope] = CCslopevsTe(heaves(he), 0,'none');
        for ANid = 1:length(angles)
            ANid
            for IFid=1:length(infill)
                for CTid=1:length(crusts)
                    [  TES(he,ANid,IFid,CTid).slope TES(he,ANid,IFid,CTid).topo TES(he,ANid,IFid,CTid).dist TES(he,ANid,IFid,CTid).topodist] = ...
                        CCslopesLooper(heaves(he),angles(ANid),crusts(CTid),infill(IFid),depths,tes);
                    TES(he,ANid,IFid,CTid).crust = crusts(CTid);
                    TES(he,ANid,IFid,CTid).heave = heaves(he);
                    TES(he,ANid,IFid,CTid).angle = angles(ANid);
                    TES(he,ANid,IFid,CTid).infill = infill(IFid);     
% % %                     BigMat.Slope(:, BMcount:length(tes)+BMcount-1) = TES(he,ANid,IFid,CTid).slope;
% % %                     BigMat.Topo(:, BMcount:length(tes)+BMcount-1) = TES(he,ANid,IFid,CTid).topo;
% % %                     BigMat.Dist(:, BMcount:length(tes)+BMcount-1) = TES(he,ANid,IFid,CTid).dist;
% % %                     BMcount = length(tes)+BMcount;
                 end
            end
        end
    end
toc