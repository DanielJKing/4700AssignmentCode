%% Daniel King 100921117, Part 2b,c,d
%First plot is current vs. meshsize. Meshsize refers to number of total
%points used. 10 data points are used. 10 Current values and 10 values for
%mesh size ranging from 100 total pts to 10000 total pts. This is somewhat
%slow, and should take about 30s. The absolute value of current increases
%with more pts. It is asymptotic towards an absolute value of 0.2. For
%bottle neck area, a fixed grid of 60x60 pts is used. The area of a single
%box in this plot ranges from 2x2 pts to 20x20 pts. The current stays
%roughly constant at an absolute value of 0.2 as the bottleneck increases.
%For the bottleneck conductivity, at fixed grid of 40x40 pts is used. In
%this plot the box conductivity ranges from 0.01 to 0.4. As bottleneck
%conductivity increases, current at the interface increases. This is
%because the region is less resistive, allowing more current to flow.


%graph of current vs. Mesh size

curr = zeros(10:1);
meshsize = zeros(10:1);

for i=1:10 
    curr(i) = Ass24700part2(10*i,10*i,1,0.01,4,4,0);
end

for i=1:10 
    meshsize(i) = 100*i*i;
end


figure(1)
plot(meshsize,curr);
xlabel('Meshsize (total pts)')
ylabel('Current (A)')
title('Current vs. Meshsize')

%graph of current vs. bottleneck areas

for i=1:10 
    curr(i) = Ass24700part2(60,60,1,0.01,2*i,2*i,0);
end

for i=1:10 
    bottleneckArea(i) = (2*i)*(2*i);
end

figure(2)
plot(bottleneckArea,curr);
xlabel('BottleNeck Area of one Square')
ylabel('current (A)')
title('Current vs. BottleNeck Area')

%graph of current vs. bottleneck conductivity

for i=1:40 
    curr(i) = Ass24700part2(40,40,1,0.01*i,10,10,0);
end

for i=1:40 
    bottleConductivity(i) = 0.01*i;
end

figure(3)
plot(bottleConductivity,curr);
xlabel('BottleNeck Conductivity')
ylabel('current (A)')
title('Current vs. BottleNeck Conductivity')