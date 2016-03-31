function [h]=ezifft(H)
%[h]=ezifft(H);
%where H is complex for positive frequencies.
%length of H must be a n+1 where n is a power of 2.
hfft = H;
n=length(H);
if n-1 ~= power2(n-1),
  %disp('(ezifft): length not power2+1');
else
  N=2*(length(H)-1);
  hfft( N:-1:(N/2)+2 ) = conj(hfft(2:N/2));
  h=real(ifft(hfft));
end
