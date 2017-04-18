clear
clf
R=20;
C=10*10^-6;
h=12.56*10^-6;
t = 0:h:0.01;
E=1;
x_FE = zeros(1,length(t));
xn = 0;


for i=1:10000
   TF(i) = 20*log10(1/sqrt( ((i*2*pi*R*C)^2) + 1));
end

x=0:1:9999;
figure(1);
semilogx(x,TF);
axis([0 10000 -20 0]);
f=796;


% Van = E*(1-exp(-t/(R*C)));
figure(2);
%plot(t,Van)
hold on;
% 
% for n=1:length(t)
%     x_FE(n) = xn;
%     xnp1 = (1-h/(R*C))*xn + (E/(R*C))*h;
%     xn=xnp1;
% end

% plot(t,x_FE);

Esin = sin(2*pi*796*t);

for n=1:length(t)
    x_FE(n) = xn;
    xnp1 = (1-h/(R*C))*xn + (Esin(n)/(R*C))*h;
    xn=xnp1;
end
plot(t,x_FE);
hold on;
plot(t,Esin);
hold on;
Van = Esin.*(1-exp(-t/(R*C)));
plot(t,Van);

