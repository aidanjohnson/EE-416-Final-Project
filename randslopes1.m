function [x,rcvd]  = randslopes1( T, sigman,pctdD,L,Aslp,Nslps2);
%function [x,r] = randslopes1( T, sigman,pctdD,L,Aslp,Nspls2);
%over [1:T] builds L Segments of differing Slopes 
%sigman AWGN noise power, L is about the number of slopes of duration D
%+-dD. normally I am setting Aslp=1 and Nslps2=1, so that we only have 1
%slope, but randomized over sign

Slopes = Aslp*[1:Nslps2];Slopes=sort([-Slopes Slopes]);
Nslopes=length(Slopes);
%Amp0 = 0;
Slopes
%pctdD=1;T=900;
D = floor( (1)*(T/L) );dD = floor( pctdD*D );
disp(['The duration are drawn from a uniform'])
disp([' Mean and plus/minus range are here:']) 
disp([D dD])
disp('  ');
%x=zeros(2*T,1);
disp('  ');
disp([ '    el start(el) ending(el) dur(el) slope(el)']);
disp(['    No Start Ending Duration Slope']);
disp('  ');
%amp =[1:L];[tmp,sindx]=sort(rand(L,1));
%amp=amp(sindx);%randomize
start=zeros(1,L);dur=start;slope=start;ending=start;

tlast = 0;xold=0;oldslope=+1;%updownslp =[];oldslope=+1;
for el = [1:L];% Loop over Segments
    start(el) = tlast+1;
    %amp(el) = sqrt( -2*P*log(rand(1,1)) );
    dur(el) = floor(D +dD*(rand(1,1)-1/2));%random duration of Segment el
    %oldstart = start(el)+dur(el);%This is the ending Time index
    ending(el) = start(el)+ (dur(el)-1);
    
    % Slopes are selcted from a Table of [1:Nslopes] using a Markov Process
    pSlope = 0.3;% Your pSlope Value may be different
    if(oldslope == +1);UrandSlope = (1-pSlope)*rand(1,1);end
    if(oldslope == -1);  UrandSlope = pSlope+ (1-pSlope)*rand(1,1);end 
    slope(el) = Slopes(1+floor(Nslopes*UrandSlope) );
    oldslope=sign(slope(el));
    
    %if(slope(el) >0);updownslp=[updownslp +1];else;updownslp=[updownslp -1];end
    %tindx = [start(el):1:oldstart];
    %disp([el start(el) oldstart dur(el) amp(el)]);
    %x(tindx) = amp(el)*ones(size(tindx));
    tindx = [start(el):1:ending(el)];
    tlast=ending(el);
    x(tindx) = xold + slope(el)*[1:dur(el)];xold=x(end);
    disp([el start(el) ending(el) dur(el) slope(el) ]);
      slp(tindx)=slope(el);
end
x = x(:);Lx = length(x);
%The array slp contains the 

disp(' ');format short e
Lending=length(ending);


t = [1:Lx]';

w =  randn(Lx,1);w = sigman*w;
rcvd=x+w;
maxX=max(abs(x));
if(maxX>=200);disp('REJECTED,since the signal has too few transistions ');end

subplot(211);
plot(t,x,'g','LineWidth',2);hold on;
plot(t,slp*maxX/Aslp,'m','LineWidth',2);hold off;grid%scaled
title('Green Signal and Magenta Scaled Slopes')
subplot(212);
plot(t,rcvd,'b','LineWidth',2);hold on;%off;grid
plot(t,x,'g','LineWidth',1);hold on;
%plot(t,rcvd,'b','LineWidth',2);hold off;grid
xATending = x(ending);xATending = xATending(:);

plot( ending(:),xATending ,'r+','MarkerSize',20);hold on;


hold off;grid
title('Noisy Received Signal')
subplot(111)

%whos;%enable whos to look at the variables and size

%disp('READ Comments Below')
% The red markers show where each transistion takes place
% You only need to find the transisions where the slope is reversed.



%save soln-1.txt x,r -ascii
%save quest-1.txt x -ascii