function mfsim06(  snroutdb, Ksig, Lt, s )
%function mfsim05(  snroutdb, Ksig, Lt, s);
% signal s can be arbitrary, you can enter any vector
% MF simulator 
%Typically we need about 10-15 dB for detection
%mfsim05( 16, 3, 1000, randn(20,1) );
%mfsim05( 16, 3, 1000, ones(20,1) );

%normalize the signal

s = s/sqrt(s'*s);Ls=length(s);% makes the given signal unit energy
snr = 10^(snroutdb/10); 
A=1;
sigma2N  = 1/(snr);%unit variance for noise
sigmaN = sqrt(sigma2N);% standard deviation

Tsig = 1+floor( (Lt-Ls)*rand(Ksig,1));
% This is random selection, and will likely give overlapping signals

hsig = zeros(Lt,1);
hsig(Tsig) = ones(Ksig,1);%These are the delay for each signal
Lhsig = length(hsig);
sall = filter(s,1,hsig);%the filter delays and superposes ieach copy
Lsall =length(sall);

%plot(A*sall,'r+');grid;pause

x = A*sall;%input
r = x + sigmaN*randn(Lt,1);%output
t = [1:Lt]';%time index


hmf = [s(end:-1:1) ];hmf = hmf(:);Lhmf = length(hmf);

y = filter(hmf,1,r);%y = conv(hmf, r);%This does the template matching

Ly = length(y);Lr = length(r);Lhsig=length(hsig);

%y = y(Lhmf:end);Ly = length(y);
absy = abs(y);%



y0 = filter(hmf,1,x);%Noise-free output%y0 = conv(hmf, A*sall);
%y0 = y0(Lhmf:end);
Ly0 = length(y0);
%absy0 = abs(y0);%Ly

figure(1)
subplot(311);
% time index t, rcvd sig r, output signal 
plot( t,r,'r');hold on;
plot( t,x,'b','LineWidth',3);hold off;grid
title('The Problem:  Input Blue    Rcvd Red')

subplot(312)
plot( t,x,'b',t,r,'r',t,y,'g');grid
title('The MF in Action:  Input Blue    Rcvd Red  MF output Green')
%indx = find(hsig>0);plot(t(indx)+Ls,hsig(indx),'r*');hold off;grid
%title('plot r/blue y/green')

subplot(313)
plot( t,x,'b',t,y,'g');hold on;
threshld=0.5;
indxMF = find( abs(y) > threshld );

  plot( t(indxMF),x(indxMF),'b+',t(indxMF),y(indxMF),'rX','MarkerSize',10);
  hold off; grid
title('Useful Timing Plot w/Markers:  Input Blue     MF output Green')
subplot(111)

%indx = find(hsig>0);plot(t(indx)+Ls,hsig(indx),'r*');hold off;grid
%title('plot r/blue y/green')

%subplot(212);
%plot( t,r,'b',t,y0,'g','LineWidth',3);grid
%title('plot r/blue y0/green')
%subplot(111)
