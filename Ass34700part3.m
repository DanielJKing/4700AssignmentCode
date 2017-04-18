%% Daniel King 100921117, Part 3
% Electron density is slightly more dense towards the right due to the
% acceleration of the E field to the right. Random rethermalizations
% however tend to make the particle density more uniform in the map.
% To make this simulation more accurate, Hall Effect could be included. A
% magnetic field could be added which would update the Lorentz Force on
% each electron to F=qE + qvB.
%
% The temperature dependence on conductivity
% could be accounted for, making the conductivity regions a function of
% temperature. Multiple scattering mechanisms could have been used which is
% the application of Matthiessen's Rule. Another effect to improve the
% accuracy is the thermoelectric effect or Seebeck effect. Voltage 
% differences can be modelled to effect temperature. In this model  that, 
% effect is not accounted for, but would make the simulation more realistic.

clear
clf

me=0.26*(9.11*10^-31); %eff mass
charge = 1.6*10^-19; %electron charge
dt = (10^-14); %sim for 1000 dt's
pscat = 1-exp(-(dt)/(0.2*10^(-12)));
kb = 1.3806*10^-23;
vth = sqrt((2*300*kb)/(me)); %thermal velocity, 2D
N = 1000; %no of particles
steps = 1000; %no of steps
p1 = 1; %partition lower bound
p2 = 3; %partition upper bound

partition=100;
count = zeros(100,100);
T = zeros(partition,partition);
velx = zeros(partition,partition);
vely = zeros(partition,partition);
pcomp = zeros(1,steps);
vnet = zeros(1,N);
x = zeros(1,N);
y = zeros(1,N);
vx = zeros(1,N);
vy = zeros(1,N);

Ex = zeros(100,100);
Ey = zeros(100,100);
E = [Ex,Ey];

ax = zeros(100,100);
ay = zeros(100,100);

E = fields(1,0.01);

for i = 1:100
    for j = 1:100
        ax(i,j) = E(i,j)*(charge/me);
    end
end

for i = 1:100
    for j = 1:100
        ay(i,j) = E(i,j+100)*(charge/me);
    end
end

for i = 1:N
   x(i) = rand*200*(10^-9);
   y(i) = rand*100*(10^-9);                 %initial conditions
   if(x(i)>=80*10^-9 && x(i)<=120*10^-9 && y(i)<=40*10^-9)
       x(i) = x(i) + 41*10^-9;
   end                     %if randomness puts particle in rectangle, then add 41nm to put it out
   if(x(i)>=80*10^-9 && x(i)<=120*10^-9 && y(i)>=40*10^-9) 
       x(i) = x(i) + 41*10^-9;
   end       
   
   vx(i) = (vth).*randn(1,1)*(1/sqrt(2));
   vy(i) = (vth).*randn(1,1)*(1/sqrt(2));
   vnet(i) = sqrt(vx(i)^2 + vy(i)^2);
end

for d = 1:steps
    pcomp(d) = rand;
end

for j = 1:steps     %amount of timesteps

    for a=1:100
        for b=1:100
            for c=1:partition
                if( x(c)<a*2*10^-9 && x(c)>(a-1)*2*10^-9 && y(c)<b*1*10^-9 && y(c)>(b-1)*1*10^-9)
                    vx(c) = vx(c) + ax(a,b).*dt;
                    vy(c) = vy(c) + ay(a,b).*dt;
                end
            end
        end
    end
    
    x(1:N) = x(1:N) + (vx(1:N).*dt);
    y(1:N) = y(1:N) + (vy(1:N).*dt);     

    if(pscat > pcomp(j))
        
        for u = 1:N
            vx(u) = (vth).*randn(1,1)*(1/sqrt(2));
            vy(u) = (vth).*randn(1,1)*(1/sqrt(2));
        end
    
    end
    
    
    
    for k = 1:N   %boundary conditions
        if( y(k)<=0 || y(k)>=100*10^-9)
            vy(k) = -vy(k);
        end
        if(x(k)<=0)
            x(k) = x(k) + 200*10^-9;
        end
        if(x(k)>=200*10^-9)
            x(k) = x(k) - 200*10^-9;
        end
        
        if( (y(k)>=60*10^-9 || y(k)<=40*10^-9) && x(k)+dt*vx(k)+dt*dt*ax(k)>=80*10^-9 && x(k)<=80*10^-9)%left side
            vx(k) = -vx(k); 
        end
        if( (y(k)>=60*10^-9 || y(k)<=40*10^-9) && x(k)+dt*vx(k)+dt*dt*ax(k)<=120*10^-9 && x(k)>=120*10^-9)%right side
            vx(k) = -vx(k); 
        end

        if( (y(k)+dt*vy(k)+dt*dt*ay(k)>=60*10^-9 || y(k)+dt*vy(k)+dt*dt*ay(k)<=40*10^-9) && x(k)<=120*10^-9 && x(k)>=80*10^-9)
            vy(k) = -vy(k);  %Horizontal components of rectangles
        end
        
    end
    for q=p1:p2
       colorVec = hsv(5);
       plot(x(q),y(q),'.','color', colorVec(q,:));
    end     

    axis([0,200*10^-9,0,100*10^-9]);
    pause(0.0001);
    hold on;

end


    rectangle('Position',[80*10^-9 60*10^-9 40*10^-9 40*10^-9]);
    rectangle('Position',[80*10^-9 0 40*10^-9 40*10^-9]);
    
    
for s=1:partition
    for r=1:partition
        for c=1:N
            if( x(c)<r*2*10^-9 && x(c)>(r-1)*2*10^-9 && y(c)<s*1*10^-9 && y(c)>(s-1)*1*10^-9) 
                count(s,r) = count(s,r)+1;
                velx(s,r) = velx(s,r) + vx(c);
                vely(s,r) = vely(s,r) + vy(c);
            end
        end
    end
end

for v=1:partition
   for z=1:partition
       T(v,z) = (me/(2*kb))*( mean(velx(v,z).^2)+ mean(vely(v,z).^2) );
   end
end

figure(2);
surf(count);
title('Electron Density');
view(2);
colorbar;
caxis([min(min(count)),  max(max(count))]);

figure(3);
surf(T);
title('Temperature Density');
view(2);
colorbar;
caxis([min(min(T)),  max(max(T))]);

sprintf('Avg temp is %0.5e K' ,mean(T(:)))