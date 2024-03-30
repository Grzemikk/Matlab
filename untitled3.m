%% Kolos gr 1

% ZAD 1

x = load("2023_tsd_xcorr_1d.txt");
Fs = 250;
tmax = 1/Fs * (length(x) - 1);
%t = (0: 1/Fs :(length(x)-1)/Fs)';
t = 0 : 1/Fs : tmax;
plot(t,x);

%trojkat
tc = -max(t) : 1/Fs : max(t);
tt = 0 : 1/Fs : 6;
xt = 1.9*(1-abs(tt-3)/3).*(abs(tt-3)<3);
% xct = xcorr(x, xt);
% nr = find( xct >= 0.9999*max(xct), 2, 'first');
% shift = tc(nr) + 3;
xct = xcorr(x.^4, xt.^4);
nr = find( xct >= 0.99999*max(xct), 2, 'first');
shift = tc(nr)';


% gauss
tg = -3 : 1/Fs : 3; 
xg = 1.5*exp(-(tg.^2)./2);
% xcg = xcorr( x, xg );
% nrg = find( xcg >= 0.9999*max(xcg), 1, 'first' );
% shift2 = tc(nrg) + 2;
xcg = xcorr(x.^0.25, xg.^0.25) + xcorr(1-x, 1-xg);
nrg = find( xcg >= 0.99999*max(xcg), 1, 'first');
shift2 = tc(nrg) + 3;


subplot(211), plot( t,x,'r',tt+shift, xt, 'g', tg+shift2, xg, 'b' );
subplot(212), plot( t, x, 'r', tg+shift2, xg, 'g');





% % trojkat
% 
% tt = 0: 1/Fs :6;
% xt = 1.9*(1-abs(tt-3)/3);
% tc = (-max(t): 1/Fs :max(t));
% xct = xcorr(x.^4, xt.^4);
% nr = find( xct >= 0.99999*max(xct), 2, 'first');
% pos_t = tc(nr)'
% 
% % gauss
% 
% tg = -3: 1/Fs : 3;
% xg = 1.5*exp(-(tg.^2)./2);
% xcg = xcorr(x.^0.25, xg.^0.25) + xcorr(1-x, 1-xg);
% nr = find( xcg >= 0.99999*max(xcg), 1, 'first');
% pos_g = tc(nr) + 3
% 
% subplot(211), plot(t,x, 'r', tt+pos_t, xt, 'b',tg+pos_g, xg, 'g');
% subplot(212), plot(tc, xct, 'r', tc, xcg, 'b');


%%
close all; clear clc;

% ZAD 2

Fs = 250;
t = -15 : 1/Fs : 15;
x1 = (1.5/5)*( mod(t,5));
x2 = 1.4*exp( -(t+2).^2/10);
liniowa = (((12-3)/30)/2)*t + 7.5;
x3 = 0.75*sin(2*pi*t.*liniowa );
liniowa2 = ((0.4-0.7)/30)*t + 0.55;
x4 = liniowa2.*sin(2*pi*t*18);
x5 = -0.2 + 1.3*randn(size(t));
x = x1 + x2 + x3 + x4 + x5;

f0 = 18;
f = linspace( -Fs/2, Fs/2, length(t) );
BP = 1./( 1 +  ( f./(f.^2 - f0.^2)).^2);
XT = fftshift( fft(x) );
WA = abs(XT);
xn = real(ifft(ifftshift( BP .* XT )));

% f0 = 18;
% f = linspace(-Fs/2, Fs/2, length(t));


subplot(211), plot( t,x,'g', t, xn, 'r' );
subplot(212), plot( f, WA, 'g', f, BP*400, 'r')
figure();
subplot(221), plot( t, x1 );
subplot(222), plot( t, x2 );
subplot(223), plot( t, x3 );
subplot(224), plot( t, x4 );


%%
close all; clear; clc;

%ZAD 3

a = load("2023_tsd_szum_1d.txt");
t = a(:,1)';
x = a(:,2)';
xs = a(:,3)';

dt = t(2) - t(1);
Fs = 1/dt;
disp(a);
t = 0 : 1/Fs : (length(t)-1)/Fs;

XT = fftshift( fft( xs ) );
WA = abs(XT);
f = linspace( -Fs/2, Fs/2, length(t) );

f0 = 22;
BP = 1./ ( 1+( 2*f ./(f.^2 - f0.^2)).^20);
plot(f,WA,'g', f, BP*1000,'r');

xn = real(ifft(ifftshift(BP.*XT)));

min_oc = Inf*ones(1,4);
n = length(t);
for N = 3:2:101
    LP = ones(1,N)/N;
    LP2 = fspecial("gaussian",[1 N], N/5 );
    xa = conv(xn,LP,'same');
    xg = conv(xn,LP2,'same');
    xm = medfilt1(xn, N);
    xm = medfilt1(xm, N+10);
    xw = wiener2(xn, [1 N]);
    xm = wiener2(xm, [1 N+5]);

    oc = (1000/n).*sqrt(sum((x-xa).^2));
    if oc < min_oc(1)
        min_oc(1) = oc;
        ka = N;
    end

    oc = (1000/n).*sqrt(sum((x-xg).^2));
    if oc < min_oc(2)
        min_oc(2) = oc;
        kg = N;
    end

    oc = (1000/n).*sqrt(sum((x-xm).^2));
    if oc < min_oc(3)
        min_oc(3) = oc;
        km = N;
    end

    oc = (1000/n).*sqrt(sum((x-xw).^2));
    if oc < min_oc(4)
        min_oc(4) = oc;
        kw = N;
    end
end

disp(min_oc);
disp( [ka kg km kw]);


%%


close all; clear; clc;
a = load('2023_tsd_szum_1d.txt');
t = a(:,1)';
x = a(:,2)';
xs = a(:,3)';
Fs = 1/(t(2) - t(1));

xt = fftshift(fft(xs));
wa = abs(xt);
f = linspace(-Fs/2, Fs/2, length(t));
W = 2;
n = 10;
f0 = 22;

bs = 1./(1+(f*W./(f.^2-f0.^2)).^(2*n));
plot(f,wa,'r',f, 1000*bs,'g');

xs2 = real(ifft(ifftshift(bs.*xt)));

min_oc = Inf * ones(1,4);
n = length(t);
for N = 3 : 2 : 101
    LP = ones(1,N)/N;
    LP2 = fspecial('gaussian', [1 N], N/5);

    xa = conv(xs2, LP, 'same');
    xg = conv(xs2, LP2, 'same');
    xm = medfilt1(xs2, N);
    xw = wiener2(xs2, [1 N]);
    xm = wiener2(xm, [1 N]);
    xm = medfilt1(xm, N + 10);
    xm = medfilt1(xm, N);
    xm = wiener2(xm, [1 N+10]);

    oc = (1000/n) .* sqrt(sum((x-xa).^2));
    if oc < min_oc(1)
        min_oc(1) = oc;
        ka = N;
    end

        oc = (1000/n) .* sqrt(sum((x-xg).^2));
    if oc < min_oc(2)
        min_oc(2) = oc;
        kg = N;
    end
    
        oc = (1000/n) .* sqrt(sum((x-xm).^2));
    if oc < min_oc(3)
        min_oc(3) = oc;
        km = N;
    end
    
        oc = (1000/n) .* sqrt(sum((x-xw).^2));
    if oc < min_oc(4)
        min_oc(4) = oc;
        kw = N;
    end
    
end

disp(min_oc)
disp( [ ka kg km kw]);


%%



















