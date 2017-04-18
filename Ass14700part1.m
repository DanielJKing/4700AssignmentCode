%% Daniel King 100921117, Part 1
% Translation is only considered in the plane, therefore 2 degrees of
% freedom. The following equation holds:
%
% $$ \frac{1}{2}*m*vth^2 =\frac{2}{2}*k*T $$
% 
% Where m is 0.26*9.11*10^-31. Solving for vth gives:
% 
% $$ vth = (\frac{2*k*T}{m})^\frac{1}{2} $$
% 
% Mean free path is given by:
% 
% $$ \tau = vth*0.2ps $$


clear
clf

me=0.26*(9.11*10^-31); %eff mass
dt = 10*(10^-15); %sim for 1000 dt's
kb = 1.3806*10^-23;
vth = sqrt((2*300*kb)/(me)); %thermal velocity, 2D
N = 1000; %no of particles
steps = 1000; %no of steps
p1 = 1; %partition lower bound
p2 = 3; %partition upper bound

vnet = zeros(1,N);
x = zeros(1,N);
y = zeros(1,N);
vx = zeros(1,N);
vy = zeros(1,N);
T = zeros(1,N);
for i = 1:N
   x(i) = rand*200*(10^-9);
   y(i) = rand*100*(10^-9);                 %initial conditions
   vx(i) = (vth)*cos(2*pi*rand);
   vy(i) = (vth)*sin(2*pi*rand);
   vnet(i) = sqrt(vx(i)^2 + vy(i)^2);
end

for j = 1:steps     %amount of timesteps
    x(1:N) = x(1:N) + (vx(1:N).*dt);
    y(1:N) = y(1:N) + (vy(1:N).*dt);    
    
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

figure(2)
h=linspace(0,1000*dt,1000);
plot(h,T)
title('Temp vs. Time');
xlabel('Time (s)');
ylabel('Temp (K)');
axis([0,10^-14,0,400]);
