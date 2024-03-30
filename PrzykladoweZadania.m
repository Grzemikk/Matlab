%% Filtracja 1
close all; clear; clc;

f1 = 1; f2 = 3; f3 = 5;
A1 = 1/f1; A2 = 1/f2; A3 = 1/f3;
Fs = 25;
t = 0 : 1/Fs : 50;

x1 = A1*sin(2*pi*t*f1);
x2 = A2*sin(2*pi*t*f2);
x3 = A3*sin(2*pi*t*f3);
x = x1 + x2 + x3;

XT = fftshift(fft(x));
WA = abs(XT);
f = linspace( -Fs/2, Fs/2, length(t) );

% pasmowozaporowy BUtterworth H(w) = 1/( 1 + ( (w*W/( w^2 -w0^2 ))^2n ) BP
f0 = 3;
BP = 1 ./ ( 1 + (0.5*f./(f.^2 - f0.^2)).^4);

xn = real(ifft(ifftshift(BP .* XT)));

subplot(211), plot(t,x,'r',t,xn,'g');
subplot(212), plot( f, WA, 'r', f, BP*400,'g');

%% 
% ZADANIE Z ZAGADNIEN OD DWORNIKA 3.2
close all; clear; clc;
%Stwórz sygnał (t=<0,20>s, Fs=50Hz) będący sumą składowych:
%• trójkątnej (tw=3, szerokość =4s, amp=1.5);
%• prostokąnej (szerokość = 2s, amp=1, środek=7s);
%• harmonicznej (f=17Hz, amp=0.1 - cały sygnał);
%• harmonicznej dla t>=10 s o narastającej liniowo częstotliwości (f_0=1Hz, f_k=11Hz)
%i amplitudzie 0.75

Fs = 50;
t = 0 : 1/Fs : 20;
x1 = 1.5*(1-abs(t-3)/2).*(abs(t-3)<2);
x2 = 1*(abs(t-7)<1);
x3 = 0.1*sin(2*pi*t*17);
x4 = 0.75*sin(2*pi*t.*(t-9)).*(t>=10);
x = x1 + x2 + x3 +x4;
plot(t,x);


%% Odszumianie 1
close all; clear; clc;
a = load( "dem_201.txt" );
x = a(:,1)';
x_szum = a(:,2)';
Fs = 1000;
disp(a);
t = 0 : 1/Fs : (length(x) - 1)/Fs;

% transformata fuoriera sprawdzam czy jest szum okresowy, jak jest to usuwa

XT = fftshift(fft(x_szum));
WA = abs(XT);
f = linspace( -Fs/2, Fs/2, length(x_szum) );

f0 = 420;
BP = 1 ./ ( 1 + (f./(f.^2 - f0.^2)).^2);

xn = ifft(ifftshift( BP .* XT ));

min_oc = Inf*ones(1,4);
n = length(t);
for N = 3:2:101
    LP = ones(1,N)/N;
    LP2 = fspecial( 'gaussian', [1 N], N/5 );
    xa = conv( xn, LP, 'same' );
    xg = conv( xn, LP2, 'same' );
    %xm = medfilt1( xn, N );
    xw = wiener2( xn, [1 N] );

    oc = sqrt( (1/n).*sum((x - xa).^2));
    if oc < min_oc(1)
        min_oc(1) = oc;
        ka = N;
    end

    oc = sqrt( (1/n).*sum((x - xg).^2));
    if oc < min_oc(2)
        min_oc(2) = oc;
        kg = N;
    end

    % oc = sqrt( (1/n).*sum((x - xm).^2));
    % if oc < min_oc(3)
    %     min_oc(3) = oc;
    %     km = N;
    % end

    oc = sqrt( (1/n).*sum((x - xw).^2));
    if oc < min_oc(4)
        min_oc(4) = oc;
        kw = N;
    end
end

disp( min_oc );
disp( [ ka kg kw ] );

% czemu nie działa xm ??????????



%% ODSZUMIANIE 2 

a = load( "szum_101.txt" );
x = a(:,1)';
x_szum = a(:,2)';
Fs = 1000;
t = 0 : 1/Fs : (length(x_szum) - 1)/Fs;

XT = fftshift(fft(x_szum));
WA = abs(XT);
f = linspace( -Fs/2, Fs/2, length(t) );

plot( f, WA );

f0 = 310;
BP = 1 ./ ( 1 + ( 20*f./(f.^2 - f0.^2)).^4);
xn = real(ifft( ifftshift( BP .* XT )));

min_oc = Inf*ones(1,4);
n = length(t);
for N = 3 : 2 : 101
    LP = ones(1, N)/N;
    LP2 = fspecial( 'gaussian', [1 N], N/5 );
    xa = conv( xn, LP, 'same' );
    xg = conv( xn, LP2, 'same' );
    xm = medfilt1( xn, N );
    xw = wiener2(xn, [1,N]);
    xw = wiener2(xw, [1,N+9]);
    xm = medfilt1(xw, N+22);

    
    oc = sqrt((1/n).*sum((x - xa).^2));
    if oc < min_oc(1)
        min_oc(1) = oc;
        ka = N;
    end

     oc = sqrt((1/n).*sum((x - xg).^2));
    if oc < min_oc(2)
        min_oc(2) = oc;
        kg = N;
    end

     oc = sqrt((1/n).*sum((x - xm).^2));
    if oc < min_oc(3)
        min_oc(3) = oc;
        km = N;
    end

     oc = sqrt((1/n).*sum((x - xw).^2));
    if oc < min_oc(4)
        min_oc(4) = oc;
        kw = N;
    end
end

disp( min_oc );
disp([ ka kg km kw]);

subplot(211), plot(t,x, 'g', t,xm,'r');
subplot(212), plot(f, WA, 'g', f, 50000*BP, 'r');

 
%% ODSZUMIANIE 3

a = load("szum_102.txt");
x = a(:,1)';
x_szum = a(:,2)';
Fs = 500;
t = 0 : 1/Fs : ( length(x_szum)-1 )/Fs;

XT = fftshift(fft(x_szum));
WA = abs(WA);
f = linspace( -Fs/2, Fs/2, length(t) );


f0 = 155;
BP = 1 ./ ( 1 + (6*f./(f.^2 - f0.^2)).^4);
W = real( ifft( ifftshift( BP .* XT)));

n = length(t);
min_oc = Inf*ones(1,4);

for N = 3 : 2 : 101
    LP = ones(1, N)/N;
    LP2 = fspecial( 'gaussian', [1 N], N/5 );
    xa = conv( W, LP, 'same' );
    xg = conv( W, LP2, 'same' );
    xm = medfilt1( W, N );
    xw = wiener2( W, [1 N] );
    xw = medfilt1( xw, N );
    xw = wiener2( xw, [1 N] );

    oc = sqrt( (1/n)*sum((x - xa).^2));
    if oc < min_oc(1)
        min_oc(1) = oc;
        ka = N;
    end

    oc = sqrt( (1/n)*sum((x - xg).^2));
    if oc < min_oc(2)
        min_oc(2) = oc;
        ka = N;
    end

    oc = sqrt( (1/n)*sum((x - xm).^2));
    if oc < min_oc(3)
        min_oc(3) = oc;
        ka = N;
    end

    oc = sqrt( (1/n)*sum((x - xw).^2));
    if oc < min_oc(4)
        min_oc(4) = oc;
        ka = N;
    end
    w = min( min_oc );
end

disp( min_oc );
disp( [ ka kg km kw ] );
disp(w);

subplot(211), plot( t, x_szum, 'r', t, xw, 'g' );
subplot(212), plot( f, WA, 'g', f, BP*5000, 'r');

%% KORELACJA 2
close all; clear; clc;

a = load( "kor_301.txt" );
t = a(:,1)';
x = a(:,2)';
dt = t(2) - t(1);
Fs = 1/dt;
disp(a);
t = 0 : 1/Fs : max(t);

% figura -> corr -> tc -> find -> shift
tc = -max(t) : 1/Fs : max(t);

% trojkat
tt = -10 : 1/Fs : 10;
x1 = 0.75*(1-abs(tt)/10).*(abs(tt)<10);
xc1 = xcorr( x.^4, x1.^4 );
nr = find( xc1 == max(xc1), 1,"first");
shift1 = tc(nr) + 10;

% trojkat 2
tt2 = - 14 : 1/Fs : 14;
x2 = 0.7*(1-abs(tt2)/14);
xc2 = xcorr( x.^2, x2.^2 );
nr = find( xc2 == max( xc2 ), 1, 'first' );
shift2 = tc(nr) + 14;

% prostokatny 
tt3 = - 10 : 1/Fs : 10;
x3 = 0.45.*(abs(tt3)<10);
xc3 = xcorr( x, x3 ) + xcorr( 1-x.^3, 1-x3.^3 );
nr = find( xc3 == max(xc3), 1, "first" );
shift3 = tc(nr) + 10;

% prostokatny2
tt4 = - 8 : 1/Fs : 8;
x4 = 0.55.*(abs(tt4)<8);
xc4 = xcorr( x, x4 );
nr = find( xc4 == max(xc4), 1, "first" );
shift4 = tc(nr) + 8;

% gauss 1
t5 = -20 : 1/Fs : 20;
x5 = 0.65*exp( -(t5).^2./32);
xc5 = xcorr( x, x5 );
nr = find( xc5 == max( xc5 ), 3, "first" );
shift5 = tc(nr);


% subplot(211),plot(t,x, t5+shift5, xx5);
% subplot(212),plot(tc,xc);


subplot(211), plot( t, x, 'g', t5+shift5, x5, 'r' );
subplot(212), plot( tc, xc5 );





disp( shift1 - 10 );





























