function [axis_interp,slopeE_Eside,slopeW_Wside,slopeE_Wside,slopeW_Eside] = axisinterpolater(axis,ll1,ll2,slopeE_Eside,slopeW_Wside,slopeE_Wside,slopeW_Eside)
% ll1 is lat or long, ll2 is tother
dothetwist = abs(axis(1,1) - axis(end,1)) > abs(axis(1,2) - axis(end,2));
if dothetwist
    axis = fliplr(axis);
    slopeE_Eside = slopeE_Eside';
    slopeW_Wside = slopeW_Wside';
    slopeW_Eside = slopeW_Eside';
    slopeE_Wside = slopeE_Wside';
end
    
for n = 1:length(ll2)   
    tempEE = zeros(size(ll1));
    tempWW = zeros(size(ll1));
    tempWE = zeros(size(ll1));
    tempEW = zeros(size(ll1)); 
  for m = 1:2:size(axis,2)
    ay=axis(find(isfinite(axis(:,m+1))),m+1);
    ax=axis(find(isfinite(axis(:,m))),m);
        axis_interp(n,m) = interp1(ay,ax,ll2(n),'linear');
        if isfinite(axis_interp(n,m))
            axis_interp(n,m+1) = ll2(n);
        else
            axis_interp(n,m+1) = NaN;
        end 
        [nul intind] = min(abs(ll1 - axis_interp(n,m)));
        if dothetwist
            plot(axis_interp(n,m+1),axis_interp(n,m),'.k')
        else
            plot(axis_interp(n,m),axis_interp(n,m+1),'.k')
        end
        if isfinite(nul)
            % EE
            temp = ll1*0 + 1;
            temp(1:intind) = 0;
            tempEE = tempEE + temp;
            % WW
            temp = ll1*0 + 1;
            temp(intind:end) = 0;
            tempWW = tempWW + temp;
            % WE
            temp= ll1*0 + 1;
            temp(1:intind) = 0;
            tempWE = tempWE + temp;    
            % EW
            temp = ll1*0 + 1;
            temp(intind:end) = 0;  
            tempEW = tempEW + temp;
        end 
  end
    slopeE_Eside(n,find(tempEE == 0)) = NaN;
    slopeW_Wside(n,find(tempWW == 0)) = NaN;
    slopeW_Eside(n,find(tempWE == 0)) = NaN;
    slopeE_Wside(n,find(tempEW == 0)) = NaN; 
    %%%
%     slopeE_Eside(n,1:intind) = NaN;
%     slopeW_Wside(n,intind:end) = NaN;
%     slopeW_Eside(n,1:intind) = NaN;
%     slopeE_Wside(n,intind:end) = NaN;   
    %%%
end

if dothetwist
    for m = 1:2:size(axis_interp,2)
        axis_interp(:,m:m+1) = fliplr(axis_interp(:,m:m+1));
    end    
    slopeE_Eside = slopeE_Eside';
    slopeW_Wside = slopeW_Wside';
    slopeW_Eside = slopeW_Eside';
    slopeE_Wside = slopeE_Wside';
end
hold off
disp('*******************************************************************');
disp('*******************************************************************');
disp('*******************************************************************');
disp('*******************************************************************');
disp(sprintf('Axis has %g segment(s), total of %g points that means %.2f %% overlap',m/2+.5,numel(find(isfinite(axis_interp)))/2,(1-((numel(find(isfinite(axis_interp)))/2)/length(ll2)))*100))
