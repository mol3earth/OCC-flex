function [R_ChiSqX, R_ChiSqTopo] = errorCalc(fittingtopo,fittingtopodist,fittingrot,fittingdist,fulldist,fulltopo,heaves,TES,tes,angles,infill,crusts,offdistX,accept_he,accept_of)
%,  R_sqrX, R_sqrY, SumResSqr
% find deviations from averages, which won't change for each R^2 calculation
y_barX = (1/length(fittingrot))*sum(fittingrot);
y_barY = (1/length(fittingdist))*sum(fittingdist);
TotSqX= (fittingrot - y_barX).^2;   
TotSqY= (fittingdist - y_barY).^2; 
%% initialize residual matrices;
clear R_sqrX
R_sqrX(1:length(heaves),1:length(angles),1:length(infill),1:length(crusts),1:length(offdistX), 1:length(tes)) = 1e8; 
SumResSqr = R_sqrX;
R_sqrY = R_sqrX;
R_ChiSqX = R_sqrX;
R_ChiSqTopo = R_sqrX;
clear ResSqX ChiSqX res_sqr ResSqY ChiSqTopo
ChiSqTopo = zeros(1,length(fittingtopodist));
ChiSqX = zeros(1,length(fittingdist));
% display error value, this is for chi^2 test
% error = 5; % degrees for slope. this is likely very conservative
errorS = 5;
errorT = 250;
disp(sprintf('Slope Error is %g^o',errorS))
disp(sprintf('Topo Error is %g m',errorT))
tic
aheid = 1;
aofid = 1;
for he = 1:length(heaves) % prev 'te'
    %if aheid <= length(accept_he)
    %if accept_he(aheid) == he
   % aheid = aheid+1;   
    disp(sprintf('he = %g',heaves(he)))    
% %      figure(100002);clf
    for teidx = 1:length(tes) %(size(te(he).slope,2))

%        disp(sprintf('/tte = %g',tes(teidx)))
        for ANid = 1:length(angles)
%            disp(sprintf('/tte = %g',angles(ANid)))
            for IFid=1:length(infill)
                for CTid=1:length(crusts)
                    for offdistidx = 1:length(offdistX)
                        % if aofid <= length(accept_of)
                        % if accept_of(aofid) == offdistidx
                        % aofid = aofid+1;   
% % %                         disp(sprintf('Offset = %g',offdistX(offdistidx)))   
                        % offset our distance vector
                        % testtedist = te(he).dist + offdistX(offdistidx);
                        testtedistTOPO = TES(he,ANid,IFid,CTid).topodist(:,teidx) + offdistX(offdistidx);%linspace(0,TES(he,ANid,IFid,CTid).dist(end,teidx),length(TES(he,ANid,IFid,CTid).topo)) + offdistX(offdistidx);
                        testtedist =  TES(he,ANid,IFid,CTid).dist(:,teidx) + offdistX(offdistidx);
                        %theseslopes = te(he).slope(:,teidx);
                        theseslopes =  TES(he,ANid,IFid,CTid).slope(:,teidx);
                        % we need to fix our topo to the elevation of the
                        % first point, which is either the fulltopo 1 or
                        % end, depending on if on left or right side.
                      %  thesetopos =  TES(he,ANid,IFid,CTid).topo(:,teidx) - ...
                      %         abs(min(TES(he,ANid,IFid,CTid).topo(:,teidx)) +  abs(min([ fulltopo(1) fulltopo(end)])));
                        % to fix according to the chosen pts minimizatoin...
                             [val, idx] = min(abs( testtedistTOPO - min(fittingdist) )); %(2:end-1)
                        thesetopos =  TES(he,ANid,IFid,CTid).topo(:,teidx) - ...                    % 2             end-1
                            abs(min(TES(he,ANid,IFid,CTid).topo(idx,teidx)) + abs(min([ fittingtopo(1) fittingtopo(end)])));

                        %% now find residual for each point in slope space
                        for jkkk = 1:length(fittingdist)
                            % first find the nearest point to the point, 
                            [val, idx] = min(abs( testtedist - fittingdist(jkkk)));
                            
                            %minslope = te(he).slope(idx,teidx);
                            minslope = TES(he,ANid,IFid,CTid).slope(idx,teidx);
                            
                            
% % %                             ResSqX(jkkk) = (fittingrot(jkkk) - minslope)^2;
                            ChiSqX(jkkk) = ( ( fittingrot(jkkk) - minslope ) / errorS)^2;           
% % %                             res_sqr(jkkk) = (fittingrot(jkkk) - minslope)^2;
% % %                             [val, idx] = min(abs(theseslopes  - fittingrot(jkkk)));
% % %                             ResSqY(jkkk) = (fittingdist(jkkk) - testtedist(idx))^2;     
                        end    
                        %sum the residuals
                        %R_sqrX(teidx,offdistidx,he) = 1 - sum(ResSqX)/sum(TotSqX);   
                        %SumResSqr(teidx,offdistidx,he) = sum(res_sqr)/length(fittingdist);
                        %R_sqrY(teidx,offdistidx,he) = 1 - sum(ResSqY)/sum(TotSqY);    
                        %R_ChiSqX(teidx,offdistidx,he) = sum(ChiSqX);
                       % R_sqrX(he,ANid,IFid,CTid,offdistidx,teidx) = 1 - sum(ResSqX)/sum(TotSqX);   
                       % SumResSqr(he,ANid,IFid,CTid,offdistidx,teidx) = sum(res_sqr)/length(fittingdist);
% % %                         R_sqrY(he,ANid,IFid,CTid,offdistidx,teidx) = 1 - sum(ResSqY)/sum(TotSqY);    
                        R_ChiSqX(he,ANid,IFid,CTid,offdistidx,teidx) = sum(ChiSqX);
                        ChiSqX = ChiSqX*0;
%plot for checking
%     plot(fulldist,fulltopo,'*m')
%     plot(testtedistTOPO,...
%        thesetopos,'b')
%     hold on
%     plot(fittingtopodist,fittingtopo,'r')
                        
                       
                        %% find residual for each point in topography space
                        for jkkk = 1:length(fittingtopodist)
                            [valT, idxT] = min(abs( testtedistTOPO - fittingtopodist(jkkk)));
                            ChiSqTopo(jkkk) = ( ( fittingtopo(jkkk) - thesetopos(idxT) ) / errorT)^2;
%   plot(testtedistTOPO(idxT),thesetopos(idxT),'b*')
%   plot(fittingtopodist(jkkk),fittingtopo(jkkk),'r*')
                            
                        end
%    drawnow
%   disp(sprintf('For offset: %g , misfit is: %g',offdistX(offdistidx),sum(ChiSqTopo)))
%   pause(1)
                        R_ChiSqTopo(he,ANid,IFid,CTid,offdistidx,teidx) = sum(ChiSqTopo);
                        ChiSqTopo = ChiSqTopo*0;
                        
                      %  end
                      %  end
                    end
                end
            end
        end
    end
    %end
    %end
end 
toc

