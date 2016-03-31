function [numDF ys modelfitsbysigma] = plotErrorSurface(R_ChiSqX,minfitidx,tes,heaves,angles,infill,crusts,offdistX,HEidx,ANidx,IFidx,CTidx,offdistidx,teidx,fet,figstring)
% first make tes in Km
tes=tes*1e-3;
%%
R_ChiSqX = squeeze(R_ChiSqX);

% next, make a structure for looping
numDF =0;
numHe = length(heaves);
if numHe > 1
    numDF = numDF+1;
    ys(numDF).ylab = 'Heave (km)';
    ys(numDF).ys = heaves;
    ys(numDF).miny = heaves(HEidx);
    ys(numDF).minidx = HEidx;
end
numAn = length(angles);
if numAn > 1
    numDF = numDF+1;
    ys(numDF).ylab = 'Angle (degrees)';
    ys(numDF).ys = angles;
    ys(numDF).miny = angles(ANidx);
    ys(numDF).minidx = ANidx;
end
numIF = length(infill);
if numIF > 1
    numDF = numDF+1;
    ys(numDF).ylab = 'Infill (km)';
    ys(numDF).ys = infill;
    ys(numDF).miny = infill(IFidx);
    ys(numDF).minidx = IFidx;
end
numCr = length(crusts);
if numCr > 1
    numDF = numDF+1;
    ys(numDF).ylab = 'Crust (km)';
    ys(numDF).ys = crusts;
    ys(numDF).miny = crusts(CTidx);
    ys(numDF).minidx = CTidx;
end
numOS = length(offdistX);
if numOS > 1
    numDF = numDF+1;
    ys(numDF).ylab =  'Offset (km)';
    ys(numDF).ys = offdistX;
    ys(numDF).miny = offdistX(offdistidx);
    ys(numDF).minidx = offdistidx;
end
numTe = length(tes);
if numTe > 1
    numDF = numDF+1;
    ys(numDF).ylab = 'Te (km)';
    ys(numDF).ys = tes;
    ys(numDF).miny = tes(teidx);
    ys(numDF).minidx = teidx;
end

%% Sigma stuff
% do sigma vectors
sigmaMat = [        
    1.00    4.00    9.00
    2.30    6.18    11.83
    3.53    8.02    14.16
    4.72    9.72    16.25
    5.89    11.31   18.21
    7.04    12.85   20.06 
    8.18    14.34   21.85 ];
foursig=[16.00
19.33
22.06
24.50
26.77
28.91
30.96];
sigmas = [  sigmaMat(numDF,:) ];
sigmalabels = {'\sigma_1' '\sigma_2' '\sigma_3'};
% http://www.reid.ai/2012/09/chi-squared-distribution-table-with.html
%Sigma      1?      1.28	1.64	1.96	2?      2.58	3?      3.29	4?
%CI %       68.3%	80%     90%     95%     95.45%	99%     99.73%	99.9%	99.99%
%P-value	0.317	0.20	0.10	0.05	0.0455	0.01	0.0027	0.001	0.00006
%chi2(k=1)	1.00	1.64	2.71	3.84	4.00	6.63	9.00	10.83	16.00
%chi2(k=2)	2.30	3.22	4.61	5.99	6.18	9.21	11.83	13.82	19.33
%chi2(k=3)	3.53	4.64	6.25	7.81	8.02	11.34	14.16	16.27	22.06
%chi2(k=4)	4.72	5.99	7.78	9.49	9.72	13.28	16.25	18.47	24.50
%chi2(k=5)	5.89	7.29	9.24	11.07	11.31	15.09	18.21	20.52	26.77
%chi2(k=6)	7.04	8.56	10.64	12.59	12.85	16.81	20.06	22.46	28.91
%chi2(k=7)	8.18	9.80	12.02	14.07	14.34	18.48	21.85	24.32	30.96
%chi2(k=8)	9.30	11.03	13.36	15.51	15.79	20.09	23.57	26.12	32.93
%chi2(k=9)	10.42	12.24	14.68	16.92	17.21	21.67	25.26	27.88	34.85
%chi2(k=10)	11.54	13.44	15.99	18.31	18.61	23.21	26.90	29.59	36.72
disp(sprintf('%g models or %.3f%% are within 1 sigma',numel(find(R_ChiSqX<sigmas(1))),100*numel(find(R_ChiSqX<sigmas(1)))/numel(R_ChiSqX)))
disp(sprintf('%g models or %.3f%% are within 2 sigma',numel(find(R_ChiSqX<sigmas(2))),100*numel(find(R_ChiSqX<sigmas(2)))/numel(R_ChiSqX)))
disp(sprintf('%g models or %.3f%% are within 3 sigma',numel(find(R_ChiSqX<sigmas(3))),100*numel(find(R_ChiSqX<sigmas(3)))/numel(R_ChiSqX)))



%% now lets histogram the variables, 



szRsq = size(R_ChiSqX);
for n = 1:length(sigmas)
    fig=figure(11110*n);clf;
    idxs = find(R_ChiSqX<sigmas(n));    
    if numDF == 1
        [ id1 ] = ind2sub(szRsq , idxs);
        modelfitsbysigma(n).ids = id1;    
    elseif numDF == 2
        [ id1, id2 ] = ind2sub(szRsq , idxs);
        modelfitsbysigma(n).ids = [id1  id2];    
    elseif numDF == 3
        [ id1, id2, id3 ] = ind2sub(szRsq , idxs);
        modelfitsbysigma(n).ids = [id1  id2 id3 ];   
    elseif numDF == 4
        [ id1, id2, id3 , id4] = ind2sub(szRsq , idxs);
        modelfitsbysigma(n).ids = [id1  id2 id3  id4]; 
    elseif numDF == 5        
        [ id1, id2, id3, id4, id5 ] = ind2sub(szRsq , idxs);
        modelfitsbysigma(n).ids = [id1  id2 id3  id4  id5]; 
    elseif numDF == 6
        [ id1, id2, id3, id4, id5 , id6] = ind2sub(szRsq , idxs);
        modelfitsbysigma(n).ids = [id1  id2 id3  id4  id5 id6];         
    elseif numDF == 7
        [ id1, id2, id3, id4, id5 , id6 ,id7] = ind2sub(szRsq , idxs);
        modelfitsbysigma(n).ids = [id1  id2 id3  id4  id5 id6 id7]; 
    end
    hold on
    for m = 1:numDF
        sp = subplot(numDF,1,m);
        counts = hist(ys(m).ys(modelfitsbysigma(n).ids(:,m)),ys(m).ys);
        brs = bar(ys(m).ys,counts);
        brs.FaceColor = [ .05 .65 .05 ];
        for j = 1:length(counts)
            if counts(j) > 0
                t=text(ys(m).ys(j),counts(j),sprintf('%g',counts(j)),'Fontweight','bold');
                t.VerticalAlignment = 'bottom';
                t.HorizontalAlignment = 'center';
            end
        end
        xlabel(ys(m).ylab)
        ylabel('Counts')
        sp.XTick = ys(m).ys;
        sp.XTickLabel = ys(m).ys;
        if m == 1
            title(sprintf('Feature %s: %s Histogram',fet,sigmalabels{n}))
        end
    end
    
    cd mfiles
    save2pdf(sprintf('%s_Feature_%s_Histogram_%g',figstring,fet,n),fig)
    cd ..

end


%%
cmap = [ .88 .99 .66;[.88 .99 .66]*.75; [.88 .99 .66]*.5  ];
for n=1:numDF
spi=1;
fig=figure(1111110*n);clf;
% %[HEidx,ANidx,IFidx,CTidx,offdistidx,teidx]
%  1  2  3  4  5  6
%  he an if ct os te
% heaves % angles % infill % crusts % offdistX  % tes
nys = ys(n).ys ; miny = ys(n).miny ; ylab = ys(n).ylab;
for m = 1:(numDF)
    % check that we are not going to plot the same variable against itself
    if m~=n
    clear C h cl cb
    idxs = setdiff(1:numDF,[n m]);
    pR_ChiSqX=permute(R_ChiSqX,[m,n,idxs]);
    if numDF == 7
    thispR = pR_ChiSqX(:,:,minfitidx(idxs(1)),minfitidx(idxs(2)),minfitidx(idxs(3)),minfitidx(idxs(4)),minfitidx(idxs(5)));
    elseif numDF == 6
    thispR = pR_ChiSqX(:,:,minfitidx(idxs(1)),minfitidx(idxs(2)),minfitidx(idxs(3)),minfitidx(idxs(4)));
    elseif numDF == 5
    thispR = pR_ChiSqX(:,:,minfitidx(idxs(1)),minfitidx(idxs(2)),minfitidx(idxs(3)));
    elseif numDF == 4        
    thispR = pR_ChiSqX(:,:,minfitidx(idxs(1)),minfitidx(idxs(2)));
    elseif numDF == 3        
    thispR = pR_ChiSqX(:,:,minfitidx(idxs(1)));
    end
    
    subplot(numDF-1,1,spi)
    [C h] = contourf(nys,ys(m).ys, thispR ,...
        [ sigmas] );
    cl = clabel(C);
    for cln = 1:2:length(cl)
        cl(cln).Marker = '.';
    end
    for cln = 2:2:length(cl)
        for sigstr = 1:length(sigmas)
            if strfind(num2str(sigmas(sigstr)),cl(cln).String)
                cl(cln).String = sigmalabels{sigstr};
            end           
            cl(cln).HorizontalAlignment ='center'; 
            cl(cln).FontSize = 18;
        end
    end
    hold on
    plot(miny,ys(m).miny,'+r')
    colormap(cmap);
    view(0, 90)
    if m == 1
        title(sprintf('Feature %s  %s vs %s',fet,ylab,ys(m).ylab))
    else
        title(sprintf(' %s vs %s',ylab,ys(m).ylab))
    end
    ylabel(ys(m).ylab)
    xlabel(ylab)
    spi=spi+1;
    end
end

cd mfiles
save2pdf(sprintf('%s_Feature_%s_Slope_ErrorSurface_%g',figstring,fet,n),fig)
cd ..

end