%% Daniel King 100921117, Part 1
% a) For a constant E field, using parallel plate capacitor model, V=Ed.
% If 0.1V is across the x direction which has a distance of 200nm,
% E = 0.1V/200nm = 500,000V/m.
%
% b) For the electric field, F=qE. Solving for F gives 
% F=(1.6x10^-19C)*(500,000V/m) = 8x10^-14N.
%
% c) From Newton's 2nd law, F=ma. Combining with F=qE gives qE=ma.
% Solving for acceleration gives a=qE/m. 
% m is effective mass of electron, and q is charge. Therefore
% a = (1.6x10^-19C)*(500,000V/m)/(0.26*9.11*10^-31kg) = 3.38x10^17m/s^2.
% Electron velocity increased by acceleration by adding acceleration*dt
% to velocity incrementally.
%
% d) Drift current density is given by J=qnuE, where u is electron
% mobility, and n is electron concentration of 10^-15cm^-2. Drift velocity
% is v=uE. Combining these gives J=qnv. This relates drift current density
% to avg. carrier velocity. Current increases linearly due to acceleration.
% As the acceleration speeds up the particles, v increases, increasing J
% and therefore I. When scaterring occurs, velocities are re thermalized,
% setting current back to steady state velocity value. Acceleration then
% increases velocity, increasing current until scattering occurs again. The
% process repeats resulting in a saw-tooth like I vs. dt graph.




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

Vx = 0.1;
Ex = Vx/(2*10^-7);
Fx = (1.6*10^-19)*Ex;
ax = Fx/me;

Vy = 0;
Ey = Vy/(2*10^-7);
Fy = (1.6*10^-19)*Ey;
ay = Fy/me;
Ix = zeros(1,steps);

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

for i = 1:N
   x(i) = rand*200*(10^-9);
   y(i) = rand*100*(10^-9);                 %initial conditions   
   vx(i) = (vth).*randn(1,1)*(1/sqrt(2));
   vy(i) = (vth).*randn(1,1)*(1/sqrt(2));
   vnet(i) = sqrt(vx(i)^2 + vy(i)^2);
end

for d = 1:steps
    pcomp(d) = rand;
end

for j = 1:steps     %amount of timesteps
    x(1:N) = x(1:N) + (vx(1:N).*dt);
    y(1:N) = y(1:N) + (vy(1:N).*dt);    
   
    if(pscat > pcomp(j))
        
        for u = 1:N
            vx(u) = (vth).*randn(1,1)*(1/sqrt(2));
            vy(u) = (vth).*randn(1,1)*(1/sqrt(2));
        end
    
    end
    
    
    
    vx(1:N) = vx(1:N) + ax*dt;
    vy(1:N) = vy(1:N) + ay*dt;
    Jx = (10^15)*(1.6*10^-19)*(mean(vx));
    Ix(j) = Jx*(10^-5);
    
   
    
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
        
       
        
    end
    for q=p1:p2
       colorVec = hsv(5);
       plot(x(q),y(q),'.','color', colorVec(q,:));

    end     

    axis([0,200*10^-9,0,100*10^-9]);
    pause(0.00001);
    hold on;

end

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

figure(4);
h=linspace(0,steps*dt,steps);
plot(Ix);
% axis([0,1000,280,320]);
title('Ix vs. Time steps');
xlabel('Time Steps');
ylabel('Ix (A)');

sprintf('Avg temp is %0.5e K' ,mean(T(:)))

