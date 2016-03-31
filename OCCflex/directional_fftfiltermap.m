function [FILTMAP]= directional_fftfiltermap(MAP,LAT,LONG,CUTOFFX,GRADATIONX,CUTOFFY,GRADATIONY)
[numrows numcolumns]=size(MAP);
% if numrows or numcolumns not divisible by 2, fix it
oldrows = numrows;
oldcolumns = numcolumns;
if mod(numrows,2)
    numrows = numrows+1;
end
if mod(numcolumns,2)
    numcolumns = numcolumns+1;
end

% find distance covered in each dimension
longx=ll2m([LAT(1) LAT(1)],[LONG(1) LONG(end)])*1e-3;
longy=ll2m([LAT(1) LAT(end)],[LONG(1) LONG(1)])*1e-3;

% create vectors for plotting spectrum
X = -[-longy:longy/numrows*2:0];
Y = -[-longx:longx/numcolumns*2:0];

% convert cutoffs (in wavelengths duh) into frequencies
WHY=1/CUTOFFY;           %smaller cut-off frequency (1/wavelength in Km) e. g. 0.01, i.e. 1/100 Km
SHY=1/GRADATIONY;            %greater cut-off frequency (1/wavelength in Km) e. g. 0.012, i.e. 1/83.3 Km
WHX=1/CUTOFFX;
SHX=1/GRADATIONX;

% throw map into fft in Y 
fftMAPY = fft(MAP,numrows,1);
FILTMAPY = fft_me(fftMAPY,numrows,longy,WHY,SHY,oldrows,oldcolumns);
% Throw this result into fft for X 
fftMAPYX = fft(FILTMAPY,numcolumns,2);
FILTMAPYX = fft_me(fftMAPYX,numcolumns,longx,WHX,SHX,oldrows,oldcolumns);

% Now do the reverse
fftMAPX = fft(MAP,numcolumns,2);
FILTMAPX = fft_me(fftMAPX,numcolumns,longx,WHX,SHX,oldrows,oldcolumns);
% Throw this result into fft for X 
fftMAPXY = fft(FILTMAPX,numrows,1);
FILTMAPXY = fft_me(fftMAPXY,numrows,longy,WHY,SHY,oldrows,oldcolumns);

sum(sum(abs(FILTMAPXY-FILTMAPYX)));
sum(sum(abs(FILTMAPXY)));
sum(sum(abs(FILTMAPYX)));
sum(sum(abs(FILTMAPX)));
sum(sum(abs(FILTMAPY)));
FILTMAP = FILTMAPXY;



% Plot The Things
figure
clf
subplot(2,2,1)
surf(LONG,LAT,FILTMAPYX)   %and the spectrum is drawn only for visualization
view([0 90]);shading interp;colorbar;axis equal
lightangle(-90,1e-3)
lightangle(0,1e-3)
xlim([LONG(1) LONG(end)])
ylim([LAT(1) LAT(end)])
title('Y then X Combine Filtered Map') ...Bathymetry map') %this is the title of the new graph

subplot(2,2,2)
surf(LONG,LAT,FILTMAPY)   
view([0 90]);shading interp;colorbar;axis equal
lightangle(-90,1e-3)
lightangle(0,1e-3)
xlim([LONG(1) LONG(end)])
ylim([LAT(1) LAT(end)])
title(sprintf('Map filtered in Y- direction only\n cosine filtered Cut: %.1f %.1f',CUTOFFY,GRADATIONY))

subplot(2,2,4)
surf(LONG,LAT,FILTMAPX)   
view([0 90]);shading interp;colorbar;axis equal
lightangle(-90,1e-3)
lightangle(0,1e-3)
xlim([LONG(1) LONG(end)])
ylim([LAT(1) LAT(end)])
title(sprintf('Map filtered in X-direction only\n cosine filtered Cut: %.1f %.1f',CUTOFFX,GRADATIONX)) %this is the title of the new graph

subplot(2,2,3)
surf(LONG,LAT,FILTMAPXY)   %and the spectrum is drawn only for visualization
view([0 90]);shading interp;colorbar;axis equal
lightangle(-90,1e-3)
lightangle(0,1e-3)
xlim([LONG(1) LONG(end)])
ylim([LAT(1) LAT(end)])
title('X then Y Combine Filtered Map')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION 
%%%% The real fft filtering goes on here. 
%%
function FILTeMAP = fft_me(FFTMAP,NUM_CR,Tlength,HighCut,LowCut,oldrows,oldcolumns)

for fg=1:((NUM_CR/2)+1);
   frequency(fg)=sqrt(((fg-1)/Tlength)^2);
end

%the matrix of the negative frequencies is also computed
frequency2=fliplr(frequency);
frequencytotal=[frequency frequency2];
frequencytotal(1)=[];
frequencytotal(end)=[];

frequencytotalplot=frequencytotal.*(2*pi);  %the frequency (1/wavelength) matrix is transformed to wavenumber (2*pi/wavelength) matrix
% and now filtered
FILTER=zeros(size(FFTMAP));      %the filter matrix is set to zero
for f=1:size(FFTMAP,1);
   for g=1:size(FFTMAP,2);
        if  length(frequencytotal) == size(FFTMAP,2)
            fg = g;
            DIRECTION = 2;
        else
            fg = f;
            DIRECTION = 1;
        end
      if frequencytotal(fg)<HighCut
      FILTER(f,g)=1;  
      elseif frequencytotal(fg)<LowCut
      FILTER(f,g)=.5.*(1+cos((((2*pi)*frequencytotal(fg))-(2*pi*HighCut))/(2*(LowCut-HighCut))));
      else
      FILTER(f,g)=0;
      end
   end;
end;

% finally apply the filter, and calculate the ifft
FILTeMAP = -abs(ifft(FILTER.*FFTMAP,NUM_CR,DIRECTION));
FILTeMAP = FILTeMAP(1:oldrows,1:oldcolumns);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% OLD WAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now fileiter in Xx 
% %the matrix with the frequencies of every harmonic is computed
% clear frequencyX frequencyX2 frequencytotalplotX filterX
% for g=1:((numcolumns/2)+1);
%    frequencyX(g)=sqrt(((g-1)/longx)^2);
% end
% 
% 
% %the matrix of the negative frequencies is also computed
% frequencyX2=fliplr(frequencyX);
% frequencytotalX=[frequencyX frequencyX2];
% frequencytotalX(1)=[];
% frequencytotalX(end)=[];
% 
% frequencytotalplotX=frequencytotalX.*(2*pi);  %the frequency (1/wavelength) matrix is transformed to wavenumber (2*pi/wavelength) matrix
% % and now filtered
% filterX=zeros(size(fftMAPX));      %the filter matrix is set to zero
% for f=1:size(fftMAPX,1);
%    for g=1:size(fftMAPX,2);
%       if frequencytotalX(g)<WHX
%       filterX(f,g)=1;  
%       elseif frequencytotalX(g)<SHX
%       filterX(f,g)=.5.*(1+cos((((2*pi)*frequencytotalX(g))-(2*pi*WHX))/(2*(SHX-WHX))));
%       else
%       filterX(f,g)=0;
%       end
%    end;
% end;
% 
% % finally apply the filter, and calculate the ifft
% FILTMAPX = -abs(ifft(filterX.*fftMAPX,numcolumns,2));
% FILTMAPX = FILTMAPX(1:oldrows,1:oldcolumns);
% figure(2)
% clf
% surf(FILTMAPX)
% axis equal
% lightangle(-90,1e-3)
% view(0,90);shading interp;colorbar
% title('X filter')
%%%%%%%%%%%% IN Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Filitere in Y s 
% %a vector with the frequencies of every harmonic is computed
% for f=1:((numrows/2)+1);
%    frequencyY(f)=sqrt(((f-1)/longy)^2);
% end
% 
% 
% %the matrix of the negative frequencies is also computed
% frequencyY2=fliplr(frequencyY);
% frequencytotalY=[frequencyY frequencyY2];
% frequencytotalY(1)=[];
% frequencytotalY(end)=[];
% frequencytotalplotY=frequencytotalY.*(2*pi);  %the frequency (1/wavelength) matrix is transformed to wavenumber (2*pi/wavelength) matrix
% % and now filtered
% filterY=zeros(size(fftMAPY));      %the filter matrix is set to zero
% for f=1:size(fftMAPY,1)
%    for g=1:size(fftMAPY,2);
%       if frequencytotalY(f)<WHY
%       filterY(f,g)=1;  
%       elseif frequencytotalY(f)<SHY
%       filterY(f,g)=0.5.*(1+cos((((2*pi)*frequencytotalY(f))-(2*pi*WHY))/(2*(SHY-WHY))));
%       else
%       filterY(f,g)=0;
%       end
%     end;
% end;
% % finally apply the filter, and calculate the ifft
% FILTMAPY = -abs(ifft(fftMAPY.*filterY,numrows,1));
% FILTMAPY = FILTMAPY(1:oldrows,1:oldcolumns);
% figure(1)
% clf
% surf(FILTMAPY)
% axis equal
% lightangle(0,1e-2)
% view(0,90);shading interp;colorbar
% title('Y filter')