%% Daniel King 100921117 Part 3
% d) RMS value changes depending on R and C. In first case, with R=10 and
% C=10u, BW is 15.9kHz. In 2nd case, with R=50 and C=10u, BW is 3.2kHz.
% In 3rd case, with R=20 and C=20u, BW is 4kHz. It is clear that with more
% bandwidth, the noise is larger. In first the noise is basically doubled
% with roughly 4 times the BW.

clear
clf
R=20;
C=10*10^-6;
B = 1/(2*pi*R*C);

h = 10^-5;
t = 0:h:0.01;
E=1;
Erms=0;

steps = zeros(1,length(t));
x = 0;
steps1 = zeros(1,length(t));
x1 = 0;
steps2 = zeros(1,length(t));
x2 = 0;
steps3 = zeros(1,length(t));
x3 = 0;
stepsrms = zeros(1,length(t));
xrms = 0;
stepsrms1 = zeros(1,length(t));
xrms1 = 0;
stepsrms2 = zeros(1,length(t));
xrms2 = 0;
stepsrms3 = zeros(1,length(t));
xrms3 = 0;
analytic_const = E*(1-exp(-t/(R*C)));

Imax = E/(5*R);

for i=1:length(t)
    steps(i) = x;
    I = (Imax)*randn();   
    xjump = (1-h/(R*C))*x + (E/(R*C))*h + h*(I/C);
    x=xjump;
end

figure(1);
plot(t,steps);
hold on;
plot(t,analytic_const);
title('Analytic Solution no Noise and with Noise');
xlabel('Time (s)');
ylabel('Voltage (V)');

for i=1:length(t)
    stepsrms(i) = xrms;
    I = (Imax)*randn();   
    xjumprms = (1-h/(R*C))*xrms + (Erms/(R*C))*h + h*(I/C);
    xrms=xjumprms;
end

vrms = sqrt(mean(mean((stepsrms).^2)));
fprintf('RMS value with R=20, C=10u is %f\n' ,vrms);

fs=1/h;
f = linspace(0,fs,10e4);
fft_out = fft(steps,10e4)/length(t);

figure(2);
plot(f,abs(fft_out));
title('fft with R=20, C=10u');
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');
axis([0 10e3 0 100e-3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%R=10 C=10
R1=10;
C1=10*10^-6;
Imax1 = E/(5*R1);

for i=1:length(t)
    steps1(i) = x1;
    I = (Imax1)*randn();   
    xjump1 = (1-h/(R1*C1))*x1 + (E/(R1*C1))*h + h*(I/C1);
    x1=xjump1;
end
for i=1:length(t)
    stepsrms1(i) = xrms1;
    I = (Imax1)*randn();   
    xjumprms1 = (1-h/(R1*C1))*xrms1 + (Erms/(R1*C1))*h + h*(I/C1);
    xrms1=xjumprms1;
end


vrms1 = sqrt(mean(mean((stepsrms1).^2)));
fprintf('RMS value with R=10, C=10u is %f\n' ,vrms1);
figure(3);
fft_out1 = fft(steps1,10e4)/length(t);

plot(f,abs(fft_out1));
axis([0 10e3 0 100e-3]);
title('fft with R=10, C=10u');
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%R=50 C=10
R2=50;
C2=10*10^-6;
Imax2 = E/(5*R2);

for i=1:length(t)
    steps2(i) = x2;
    I = (Imax2)*randn();   
    xjump2 = (1-h/(R2*C2))*x2 + (E/(R2*C2))*h + h*(I/C2);
    x2=xjump2;
end
for i=1:length(t)
    stepsrms2(i) = xrms2;
    I = (Imax2)*randn();   
    xjumprms2 = (1-h/(R2*C2))*xrms2 + (Erms/(R2*C2))*h + h*(I/C2);
    xrms2=xjumprms2;
end


vrms2 = sqrt(mean(mean((stepsrms2).^2)));
fprintf('RMS value with R=50, C=10u is %f\n' ,vrms2);
figure(4);
fft_out2 = fft(steps2,10e4)/length(t);


plot(f,abs(fft_out2));
title('fft with R=50, C=10u');
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');
axis([0 10e3 0 100e-3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%R=20 C=20
R3=20;
C3=20*10^-6;
Imax3 = E/(5*R3);

for i=1:length(t)
    steps3(i) = x3;
    I = (Imax3)*randn();   
    xjump3 = (1-h/(R3*C3))*x3 + (E/(R3*C3))*h + h*(I/C3);
    x3=xjump3;
end
for i=1:length(t)
    stepsrms3(i) = xrms3;
    I = (Imax3)*randn();   
    xjumprms3 = (1-h/(R3*C3))*xrms3 + (Erms/(R3*C3))*h + h*(I/C3);
    xrms3=xjumprms3;
end


vrms3 = sqrt(mean(mean((stepsrms3).^2)));
fprintf('RMS value with R=20, C=20u is %f\n' ,vrms3);
figure(5);
fft_out3 = fft(steps3,10e4)/length(t);


plot(f,abs(fft_out3));
title('fft with R=20, C=20u');
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');
axis([0 10e3 0 100e-3]);