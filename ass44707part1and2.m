%% Daniel King 100921117 Parts 1 and 2
% a) Filter is a low pass filter. Transfer function is T(s)=1/(1+sRC).
%
% b) Bandwidth is given by: BW = 1/(2piRC) = 796Hz.
%
% c) Differential eqn is Vo' + (1/RC)Vo = E/RC.
%
% d) Analytic solution is Vo = E(1-e^(-t/RC)).
%
% e) Difference eqn is Vo(t+1) = (1-h/RC)Vo(t) + h(E/RC). f) Equation is
% explicit because value at future step is explicitly expressed and can be
% explicitly solved for. The future value does not appear on both sides of
% equation, which would make it implicitly dependent on itself.
%
% g) Depending on timestep chosen, dicrete solution can become unstable or
% undershoot analytic solution. If timestep is too large, value oscillates.
% Forward Euler solution is done. At a timestep of less than 2.5x10^-4,
% solution is stable.
%
% h) Rule of thumb timestep can be dt = T/100, where T is period of input.
%
% i) fft of LPF output with constant input is shown. At corner frequency
% value is an aymptote.
%
% Part 2:
%
% Johnson resistor noise eqn is given by vrms = (4kTRB)^0.5 where R=20Ohm,
% B=796Hz, k=1.38e-23J/K, and T=300K. vrms = 1.624*10^-8.

clear
clf
R=20;
C=10*10^-6;
B = 1/(2*pi*R*C);

h1 = 10^-5;
h2 = 5*10^-5;
h3 = 2.5*10^-4;
h4 = 10^-5;
t1 = 0:h1:0.01;
t2 = 0:h2:0.01;
t3 = 0:h3:0.01;
t4 = 0:h4:0.01;

E=1;
Esin = sin(2*pi*B*t4);

for i=1:10000
   TF(i) = 20*log10(1/sqrt( ((i*2*pi*R*C)^2) + 1));
end

x=0:1:9999;
figure(1);
semilogx(x,TF);
title('Frequency Response of Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
steps1 = zeros(1,length(t1));
x1 = 0;
steps2 = zeros(1,length(t2));
x2 = 0;
steps3 = zeros(1,length(t3));
x3 = 0;
steps4 = zeros(1,length(t4));
x4 = 0;

analytic_const = E*(1-exp(-t1/(R*C)));

for i=1:length(t1)
    steps1(i) = x1;
    xjump1 = (1-h1/(R*C))*x1 + (E/(R*C))*h1;
    x1=xjump1;
end

for i=1:length(t2)
    steps2(i) = x2;
    xjump2 = (1-h2/(R*C))*x2 + (E/(R*C))*h2;
    x2=xjump2;
end

for i=1:length(t3)
    steps3(i) = x3;
    xjump3 = (1-h3/(R*C))*x3 + (E/(R*C))*h3;
    x3=xjump3;
end


figure(2);
plot(t1,analytic_const);
hold on;
plot(t1,steps1);
hold on;
plot(t2,steps2);
hold on;
plot(t3,steps3);
hold on;
title('Transient Response of Filter, Const Input');
xlabel('Time (s)');
ylabel('Amplitude(V)');
figure(3);

for i=1:length(t4)
    steps4(i) = x4;
    xjump4 = (1-h4/(R*C))*x4 + (Esin(i)/(R*C))*h4;
    x4=xjump4;
end

plot(t4,steps4);
hold on;
plot(t4,Esin);
title('Transient Response of Filter, Sine Input f=796Hz');
xlabel('Time (s)');
ylabel('Amplitude(V)');

figure(4);
fft_out = fft(steps4,10e4)/length(t1);
fs=1/h1;
f = linspace(0,fs,10e4);
plot(f,abs(fft_out));
axis([0 10e3 0 100e-3]);
title('fft of Output');
xlabel('Frequency (Hz)');
ylabel('Magntidue (V)')
