%% Daniel King 100921117, Part 2
% Temperature varies for each scattering of all particles. Average for all
% time steps is rough 300K. More with more particles averge speeds tends to
% vth. Most populous bin in histogram approaches vth. Mean free path and
% mean time between collisions is calculated based on steps/count of
% scattering. With more timesteps, MFP and MTBC converges to analytic
% result.







clear
clf

me=0.26*(9.11*10^-31); %eff mass
dt = 10*(10^-15); %sim for 1000 dt's
pscat = 1-exp(-(dt)/(0.2*10^(-12)));
kb = 1.3806*10^-23;
vth = sqrt((2*300*kb)/(me)); %thermal velocity, 2D
N = 1000; %no of particles
steps = 1000; %no of steps
p1 = 1; %partition lower bound
p2 = 3; %partition upper bound
countscat=0;

pcomp = zeros(1,steps);
vnet = zeros(1,N);
x = zeros(1,N);
y = zeros(1,N);
vx = zeros(1,N);
vy = zeros(1,N);
T = zeros(1,N);
for i = 1:N
   x(i) = rand*200*(10^-9);
   y(i) = rand*100*(10^-9);                 %initial conditions
   vx(i) = (vth).*randn(1,1)*(1/sqrt(2));
   vy(i) = (vth).*randn(1,1)*(1/sqrt(2));
   vnet(i) = sqrt(vx(i)^2 + vy(i)^2);
end

for s = 1:steps
    pcomp(s) = rand;
end

for j = 1:steps     %amount of timesteps
    x(1:N) = x(1:N) + (vx(1:N).*dt);
    y(1:N) = y(1:N) + (vy(1:N).*dt);    
    
    if(pscat > pcomp(j))
        countscat = countscat+1;
        for u = 1:N
            vx(u) = (vth).*randn(1,1)*(1/sqrt(2));
            vy(u) = (vth).*randn(1,1)*(1/sqrt(2));
            
        end
    end
    
    for k = 1:N   %boundary conditions
        if(y(k)<=0 || y(k)>=100*10^-9)
            vy(k) = -vy(k);
        end
        if(x(k)<=0)
            x(k) = x(k) + 200*10^-9;
        end
        if(x(k)>=200*10^-9)
            x(k) = x(k) - 200*10^-9;
        end
  
    end
    
    for q=p1:p2
       colorVec = hsv(5);
       plot(x(q),y(q),'o','color', colorVec(q,:));

    end
    
    T(j) = (me/(2*kb))*( mean(abs(vx).^2)+ mean(abs(vy).^2) );
    
    axis([0,200*10^-9,0,100*10^-9]);
    pause(0.00001);
    hold on;
    
end

figure(2);
h=linspace(0,1000*dt,1000);
plot(T);
% axis([0,1000,280,320]);
title('Temp vs. Time steps');
xlabel('Time Steps');
ylabel('Temp (K)');


figure(3);
hist(vnet,20);
title('Freq vs. Avg Velocity');
xlabel('Avg. Velocity (m/s)');
ylabel('Freq');

mfp = (steps/countscat)*dt*mean(vnet);
meantime = dt*(steps/countscat);

sprintf('avg. velocity is %0.5e m/s' ,mean(vnet))
sprintf('mean free path is %0.5e m and mean time between collisions is %0.5e s' ,mfp,meantime)