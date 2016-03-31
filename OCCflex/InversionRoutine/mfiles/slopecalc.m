function [fyt, ndx, nnx, yt] = slopecalc(heave,te,dx,crustthick,angle,profilelength,ift)

%%%%%%%%% 1. Input parameters for fault geometry
%%% heave/fh : horizontal extension of the fault 6e4 m from Schouten
%%% Te : Effective elastic thickness in m
% dx : minimum spacing, m default 1
% crustthick : crustal thickness, m 1000 default
% angle % fault angle, degrees 60 default
% profilelength=100000;   % profile length, m
%%% fh and te are now input in CCslopevsTe, as are all parameters

%%%%%%%%% 2. call the program that calculates the fault gemoetry, dofault.m
if ift == 0
    [nx,yt,ym,ndx]=dofault(angle,heave,profilelength,dx,crustthick);
else
    [nx,yt,ym,ndx]=dofaultWinfill(angle,heave,profilelength,dx,crustthick,ift);
end

% yt is seafloor topo
% ym is mantle topo
% nx is x array ndx is spacing 

%%% 3. To solve the flexure equation we need symmetric topography
% this is *probably* because we need to find the characteristic wavelengths
% using fft
nyt=[yt fliplr(yt)];
nym=[ym fliplr(ym)];
nnx=ndx*[1:length(nym)];

% plot up the initial topo if you want
    noplot=0;
    if(noplot)
    figure(1);clf

    subplot(2,1,1)
    plot(nnx,nyt,'b-',nnx,nym,'g-')
    hold on
    end

%%%%% 4. Flexure calcuation
w=flex(nyt,nym,ndx,te);    % call the flexure program

if(noplot)
%%%%% Plot the flexural response 
plot(nnx,w,'k-')
end%%if(noplot)

% calculate the 'flexed' topographies
fym=nym+w';
fyt=nyt+w'; % subtract the flexural response from the originial topography.

% the result is the curved topography
%
%  our w, is the amount of depression, that the topography experiences
% so we take our initial topography, and add the change (w)
% We get this,
%                      .                    .
% fym ----------------~'`.                ,' `-------------------
%                         `--------------'
%                         
%
% of course the exact shape changes with the parameters                      
%
mfyt=max(fyt);
%disp('**************************************************')
%disp(['breakaway top ',int2str(mfyt-fyt(end)),' meters'])

if(noplot)
    %%%%% Plot the flexed topography
    plot(nnx,fyt,'r-',nnx,fym,'m-');
    %%fill([nnx,fliplr(nnx)],[fyt,fliplr(fym)],'-g')
    plot([nnx(1),nnx(end)],[fyt(1),fyt(end)],'-k')
    plot([nnx(1),nnx(end)],[0,0],':k')

    %%%%%%%%% 5 . Calculate the tilt
    imax=find([fyt]==max([fyt]));imax=imax(1); % find summit of fault; choose first index if the two peaks are found
    hfyt=fyt(1:end/2);
    imin=find([hfyt]==min([hfyt]));imin=imin(1); % find base of fault; choose first index if the two peaks are found
    plot(nnx(imax),nyt(imax),'r.') % plot the point
    plot(nnx(imin),nyt(imin),'r.') % plot the point
    %%axis('equal')
    %%axis([0,lp,-6000,6000])
    grid
    %legend('Initial Seafloor Topography','Initial Mantle Topography','Flexural Response','Flexed Seafloor Topography','Flexed Mantle Topography')

    tilt=rad2deg(atan(diff(fyt(imax:imax+1))/ndx)); % calculate the tilt which will be maximum at this point

    %%%%%% alternatively
    %tilt2=abs(rad2deg(atan(diff(fyt)/ndx))); 
    if(0)
    tilt2=abs(rad2deg(atan(diff(real(w))/ndx)));
    else
    tilt2=(rad2deg(atan(diff(real(w))/ndx))); 
    end
    %ifau=find(tilt2>an*.3);tilt2(ifau)=tilt2(ifau)*0;
    tilt3=[tilt2(1) tilt2'];
end
if(noplot)   
subplot(2,1,2)
if(1)
plot(nnx,fyt,'r-',nnx,fym,'m-');
axis([0,lp,-6000,6000])
grid
elseif(0) % mlarson change, aug 9, don't plot any of this.
plot(nnx,tilt3,'r.');
axis([0,lp,-70,70])
end
end%%if(noplot)
