function [H,k]=ezfft(h,dx)
%[H,k]=jfft(h,dx);
% provides positive-frequency fft of a function.
% vector k is such that dk = 2.*pi/(N*dx)
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h - vector of topography 
% dx - spacing of the topography points in the x direction
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H - 
%
%

% make sure power of 2. should have already been done though
N=power2(length(h));
huse=h;
%
% Possibly pad huse with last value, out to next power of 2:
if N~=length(h), disp('(ezfft): **Input not power of 2. Padding...'); end
%huse(length(h)+1:N)=interp1([length(h);N],[h(length(h));h(1)],length(h)+1:N);
%huse(length(h)+1:N)=table1([length(h) h(length(h));N h(1)],length(h)+1:N);

% this only does something if N~= length(h)
% otherwise, huse(length(h)...) is an empty matrix
huse(length(h)+1:N)=interp1([length(h) h(length(h))], [N h(1)],length(h)+1:N);

%disp([length(h),N,length(huse)])
%
%%% Define values for positive k
% N*dx should be 2e5
% making dk pi*1e-5
% this is our spacing in frequency space
dk = 2.*pi/(N.*dx);
% then we make a vector of half our length vector
% and scale it to our dk
k=dk .* [ 0:N./2]';
% we take half because we only take half below, of our fft

%%% now take the fast fourier transform of our topography. 
% this provides a power spectrum
hfft = fft(huse);
% we are only interested in the first half
H=hfft(1:(N./2 + 1));
%disp(length(H))
