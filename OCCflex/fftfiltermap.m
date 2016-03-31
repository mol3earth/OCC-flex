function [FILTMAP]=fftfiltermap(MAP,LAT,LONG,CUTOFF,GRADATION)
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

longx=ll2m([LAT(1) LAT(1)],[LONG(1) LONG(end)])*1e-3;
longy=ll2m([LAT(1) LAT(end)],[LONG(1) LONG(1)])*1e-3;

% create vectors for plotting spectrum

X = -[-longy:longy/numrows*2:0];
Y = -[-longx:longx/numcolumns*2:0];

WH=1/CUTOFF;           %smaller cut-off frequency (1/wavelength in Km) e. g. 0.01, i.e. 1/100 Km
SH=1/GRADATION;            %greater cut-off frequency (1/wavelength in Km) e. g. 0.012, i.e. 1/83.3 Km

fftMAP=fft2(MAP,numrows,numcolumns);  %the 2-D FFT of the gravity input matrix is computed after demeaning
bath_spectrum=abs(fftMAP);  %this computes the amplitude spectrum
%%
%the matrix with the frequencies of every harmonic is computed
for f=1:((numrows/2)+1);
   for g=1:((numcolumns/2)+1);
      frequency(f,g)=sqrt(((f-1)/longx)^2+((g-1)/longy)^2);
   end
end

%the matrix of the negative frequencies is also computed
frequency2=fliplr(frequency);
frequency3=flipud(frequency);
frequency4=fliplr(flipud(frequency));
entero=round(numcolumns/2);
if ((numcolumns/2) - entero)==0
   frequency2(:,1)=[];
   frequency3(1,:)=[];
   frequency4(:,1)=[];
   frequency4(1,:)=[];
   frequencytotal=[frequency frequency2;frequency3 frequency4];
else
   frequencytotal=[frequency frequency2;frequency3 frequency4];
end
frequencytotal(end,:)=[];
frequencytotal(:,end)=[];

frequencytotalplot=frequencytotal.*(2*pi);  %the frequency (1/wavelength) matrix is transformed to wavenumber (2*pi/wavelength) matrix
%%
%The high-cut filter is constructed
filter=frequencytotal.*0;      %the filter matrix is set to zero
for f=1:numrows;
   for g=1:numcolumns;
      if frequencytotal(f,g)<WH
      filter(f,g)=1;  
      elseif frequencytotal(f,g)<SH
      filter(f,g)=0.5.*(1+cos((((2*pi)*frequencytotal(f,g))-(2*pi*WH))/(2*(SH-WH))));
      else
      filter(f,g)=0;
      end
   end;
end;
% finally apply the filter, and calculate the ifft
FILTMAP = -abs(ifft2(fftMAP.*filter));
FILTMAP = FILTMAP(1:oldrows,1:oldcolumns);
%% Plot The Things
figure
clf
subplot(2,2,3)
surf(filter)   
view([0 90]);shading interp;colorbar
title('Amplitude spectrum of the Filter') %this is the title of the new graph

subplot(2,2,1)
surf(log10(bath_spectrum(1:numrows/2, 1:numcolumns/2)))   %and the spectrum is drawn only for visualization
view([0 90]);shading interp;colorbar
title('Amplitude spectrum of the Observed Bathymetry map') %this is the title of the new graph

subplot(2,2,2)
surf(frequencytotalplot(1:numrows/2, 1:numcolumns/2))   
view([0 90]);shading interp;colorbar
title('''Ideal'' Amplitude Spectrum of the Region')

subplot(2,2,4)
imshade(LONG,LAT,FILTMAP)   %and the spectrum is drawn only for visualization
view([0 90]);shading interp;colorbar
title(sprintf('High Pass (>%.0f km) Filtered Map',GRADATION))