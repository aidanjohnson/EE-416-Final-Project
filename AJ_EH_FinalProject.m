% Aidan Johnson, 1431797
% Emi Harada,
% EE 416 Final Project - Group 4
% 2017-12-13

clear all;
close all;

% receives random transmitted signal
%[s,r] = randslopes1(1e+03, 20, 0.9 ,20, 1, 1);

% loads given test file (receives random transmitted signal)
Rx = load('SGroup4');
T = Rx.SimData.Nsigs; % number of signal transitions/segments
N = Rx.SimData.Nsim; % total signal length
r = Rx.SimData.rcvd; % received signal
f = filter(ones(1,T)/T, 1, r); % received signal filter with moving average

% template: symmetric triangular ramp up then down (each ramp length n)
n = N/(2*T); % duration of ramp with slope of +/- 1
mf = [linspace(0,n,n+1), linspace(n-1,0,n)];
% match filters received signal with template (convolves ~ correlate)
y_mf = conv(r,mf); 

% finds maxima/minima of convolution result (concave down/up)
[maxs,t_max] = findpeaks(y_mf);
[mins,t_min] = findpeaks(-1*y_mf);

% weeds out false positive detected transitions; endpoints (initial and
% final transitions) added manually
i_max = find(t_max > n/2 & t_max < length(r)-n/2);
i_min = find(t_min > n/2 & t_min < length(r)-n/2);

% shifts concave down (maxima) and up (minima) corner locations to 
% time scale of signal
dwncrnrs = t_max(i_max(:)) - n;
upcrnrs = t_min(i_min(:)) - n;

% combines and sorts transition time locations with the initial and final
% times manually assumed to be transition points
tcorners = sort([1; dwncrnrs; upcrnrs; length(r)]);

% reconstructs transmitted signal from the corners
z = zeros(length(r),1);
for i = 1:length(tcorners)-1
    % assumes linear line between each corner (maxima/minima) location
    t1 = tcorners(i);
    t2 = tcorners(i+1);

    t = t2-t1+1; % number of sample points for line connecting corners
    zi = []; % values of signal in between corners
    
    % final corner (penultimate with respect to the final data point)
    if i == length(tcorners)-1 
        if ismember(t1,upcrnrs)
            zi = z(t1) + linspace(0,t-1,t)*1;
        elseif ismember(tcorners(i),dwncrnrs)
            zi = z(t1) + linspace(0,t-1,t)*-1;
        end        
    elseif ismember(t2,upcrnrs) % a concave up transition
        zi = z(t1) + linspace(0,t-1,t)*-1;
    elseif ismember(t2,dwncrnrs) % a concave down transition
        zi = z(t1) + linspace(0,t-1,t)*1;
    end
    
    % reconstructs signal with corner connector
    z(t1+1:t2) = zi(2:end);
end

% outputs detection points (corners or concavity transistions)
% to ASCII text file
dlmwrite('detect-n.txt',tcorners);

% outputs matched filtered signal to ASCII text file
dlmwrite('filtsignal-n.txt',z);

% plots template matched filtered received signal with locations of corners
figure;
plot(y_mf);
hold on;
scatter(t_max,y_mf(t_max));
hold on;
scatter(t_min,y_mf(t_min));
hold off;
grid on;
title('correlation of noisy signal with corner template');
legend('matched filtered','maxima','minima','Location','best');

% plots detected corners/transitions locations overlaid with received signal
figure;
subplot(1,2,1);
plot(r);
hold on;
scatter(dwncrnrs,r(dwncrnrs));
hold on;
scatter(upcrnrs,r(upcrnrs));
hold off;
grid on;
title('received signal');

% plots transition locations overlaid with transmitted signal
subplot(1,2,2);
plot(z);
hold on;
scatter(dwncrnrs,z(dwncrnrs));
hold on;
scatter(upcrnrs,z(upcrnrs));
hold off;
grid on;
title('reconstructed signal');
legend('signal','transition + to -','transition - to +','Location','best');

% plots received and trasmitted signal for comparison
figure;
subplot(2,1,1);
plot(f);
hold on;
plot(r);
hold off;
grid on;
title('comparison of signals');
legend('filtered (moving average)','received','Location','best');

% plots slopes (the useful information of the signal) of plot above
dzdt = diff(z); % slope of reconstructed signal
dfdt = diff(f); % slope of low-pass filtered signal
subplot(2,1,2);
plot(dfdt);
grid on;
title('comparison of slopes');

% plots the reconstructed signal z (the received signal matched filtered to
% remove noise)
figure;
subplot(2,1,1);
plot(z);
grid on;
title('denoised received signal');
subplot(2,1,2);
plot(dzdt);
grid on;
title('signal slope');
ylim(1.25*[min(dzdt), max(dzdt)]);
