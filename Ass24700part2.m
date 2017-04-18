%% Daniel King 100921117, Part 2a
%Function called with commented values below. Graphs for parts b,c, and d
%are in the next document with comments. Current flow at interface
%calculated by summing perpendicular components of current density at each
%of the 6 sides of the bottleneck.

function [ Ex ] = Ass24700part2(nx, ny, condA, condB, bottleWidth, bottleHeight, plot)
%if plot=1, plots. if 0, no plot
% nx = 40; %no of x points. should be even
% ny = 40; %no of y points. should be even
% bottleWidth=10; %should be even
% bottleHeight=10; %should be even
% condA = 1; %conductivity of main region
% condB = 0.01; %conductivity of bottleneck
condAvg = (condA + condB)*0.5; %conductivity at interface between main region and bottleneck
cMap = zeros(nx,ny);

for i=0:bottleHeight-1  %conductivity of interface boundry
    cMap((nx/2)-bottleWidth/2,i+1) = condAvg;
    cMap((nx/2)+bottleWidth/2,i+1) = condAvg;
    cMap((nx/2)-bottleWidth/2,ny-i) = condAvg;
    cMap((nx/2)+bottleWidth/2,ny-i) = condAvg;
end
    
for i = 0:bottleWidth %conductivity of interface boundry
    cMap((nx/2)-bottleWidth/2+i,bottleHeight) = condAvg;
    cMap((nx/2)-bottleWidth/2+i,ny-bottleHeight) = condAvg;
end

for i=0:bottleWidth-2 %conductivity inside bottleneck
    for j=0:bottleHeight-2
        cMap((nx/2)-bottleWidth/2+i+1 ,j+1) = condB;
        cMap((nx/2)-bottleWidth/2+i+1 ,nx-j-1) = condB;
    end
    cMap((nx/2)-bottleWidth/2+i+1 ,nx) = condB;
end

for i=1:nx %conductivity outside bottleneck
    for j=1:ny
        if cMap(i,j) == 0
            cMap(i,j) = condA;
        end
    end
end
           
G = zeros(nx*ny); %initialize G matrix
B = zeros(1,nx*ny); %initialize B matrix

i=1;
for j=1:ny   %left bc
   n = i + (j-1)*nx;
   B(n) = 1;
end

i=nx;
for j=1:ny  %right bc
   n = i + (j-1)*nx;
   B(n) = 1;
end

j=1;
for i=1:nx  %bottom bc
   n = i + (j-1)*nx;
   B(n) = 0;
end

j=ny;
for i=1:nx  %top bc
   n = i + (j-1)*nx;
   B(n) = 0;
end

for i = 1:nx %build G matrix using conductivities. Taken from repo.
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
        elseif j == 1
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nyp = j + 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            ryp = (cMap(i, j) + cMap(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  ny
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = j - 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j))*0.5;
            rxp = (cMap(i, j) + cMap(i + 1, j))*0.5;
            rym = (cMap(i, j) + cMap(i, j - 1))*0.5;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))*0.5;
            rxp = (cMap(i,j) + cMap(i+1,j))*0.5;
            rym = (cMap(i,j) + cMap(i,j-1))*0.5;
            ryp = (cMap(i,j) + cMap(i,j+1))*0.5;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end

V1D = G\B';  %solution of potential in 1D form

V2D = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1)*ny;
        V2D(i,j) = V1D(n); %convert potential into 2D for plot
    end
end

for i = 1:nx %calculate E using central FD
    for j = 1:ny
        if i == 1
            Ex(i,j) = (V2D(i + 1,j) - V2D(i,j));
        elseif i == nx
            Ex(i,j) = (V2D(i,j) - V2D(i - 1,j));
        else
            Ex(i, j) = (V2D(i+1,j) - V2D(i-1,j)) * 0.5;
        end
        if j == 1
            Ey(i,j) = (V2D(i,j+1) - V2D(i,j));
        elseif j == ny
            Ey(i,j) = (V2D(i,j) - V2D(i,j-1));
        else
            Ey(i,j) = (V2D(i,j+1) - V2D(i,j-1)) * 0.5;
        end
    end
end

Ex = -Ex; %from E=-Grad(V)
Ey = -Ey;

Jx = cMap .* Ex; %from J=simgaE
Jy = cMap .* Ey;

x = linspace(0,1.5,nx);
y = linspace(0,1,ny);

if plot
    figure(1)
    surf(x,y,cMap')
    xlabel('x')
    ylabel('y')
    title('Conductivity')
    colorbar
    caxis([min(min(cMap)),  max(max(cMap))]);
    view(2)

    figure(2)
    surf(x,y,V2D')
    xlabel('x')
    ylabel('y')
    title('Potential')
    colorbar
    caxis([min(min(V2D)),  max(max(V2D))]);
    view(2)

    figure(3)
    quiver(x,y,Ex',Ey')
    xlabel('x')
    ylabel('y')
    title('E Field')
    view(2)

    figure(4)
    quiver(x,y,Jx',Jy')
    xlabel('x')
    ylabel('y')
    title('Current Density')
    view(2)
end
%sum current densitites at 6 boundries of bottleneck, taking only
%perpendicular components to find current at interface
C_left_upper = sum(Jx((nx/2)-bottleWidth/2,1:bottleHeight));
C_right_upper = sum(Jx((nx/2)+bottleWidth/2,1:bottleHeight));
C_left_lower = sum(Jx((nx/2)-bottleWidth/2,ny-bottleHeight:ny));
C_right_lower = sum(Jx((nx/2)+bottleWidth/2,ny-bottleHeight:ny));

C_upper_horizontal = sum(Jy((nx/2)-bottleWidth/2:(nx/2)+bottleWidth/2,bottleHeight));
C_lower_horizontal = sum(Jy((nx/2)-bottleWidth/2:(nx/2)+bottleWidth/2,ny-bottleHeight));
Curr = C_left_upper + C_right_upper + C_left_lower + C_right_lower + C_upper_horizontal + C_lower_horizontal;

end
