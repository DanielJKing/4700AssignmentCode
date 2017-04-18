%% Daniel King 100921117, Part 1
% This is the solution of potential for 2D region. Analytic
% solution is tended towards as number of points increases. With infinite
% pts, solution is smooth and analytic solution is obtained. Due to
% computational limits, only so many pts can be used. The analytic series
% should be stopped when a relative difference between additional series
% pts becomes sufficiently small. There is no sense continuing the series
% when additional series terms barely change the potential. Numerical solution is
% beneficial because it is very quick and easy to obtain. If one desires
% quick results without tedious analytic derivation, numerical approach is
% desireable. Analytic approach is useful for obtaining physical insight
% into a problem. Much physics can be learned from analytic derivation. It
% can be very difficult and time consuming to find an analytic solution
% however.

clear
clf

nx = 30; %no of x points
ny = 30; %no of y points
Lx = 1.5; %length of x
Ly = 1; %length of y

x = linspace(0,Lx,nx); %x points
y = linspace(0,Ly,ny); %y points

G = zeros(nx*ny,nx*ny); %initialize G matrix
B = zeros(1,nx*ny); %initialize B matrix

%interior pts G and B matrix coefficients (avg of neighbours)
for i=2:nx-1
   for j=2:ny-1
        n = i + (j-1)*nx;
        G(n,n) = -4;
        G(n,n-1) = 1;
        G(n,n+1) = 1;
        G(n,n-nx) = 1;
        G(n,n+nx) = 1;
        B(n,1) = 0;
   end
end

i=1;
for j=1:ny   %left bc
   n = i + (j-1)*nx;
   G(n,n) = 1;
   B(1,n) = 0.1;
end

i=nx;
for j=1:ny  %right bc
   n = i + (j-1)*nx;
   G(n,n) = 1;
   B(1,n) = 0;
end

j=1;
for i=1:nx  %bottom bc
   n = i + (j-1)*nx;
   G(n,n) = 1;
   B(1,n) = 0;
end

j=ny;
for i=1:nx  %top bc
   n = i + (j-1)*nx;
   G(n,n) = 1;
   B(1,n) = 0;
end

V1D = G\B'; %solution of potential in 1D form

for i=1:nx
   for j=1:ny
       n = i + (j-1)*nx;
       V2D(i,j) = V1D(n); %convert potential into 2D for plot
   end
end

mesh(x,y,V2D')
xlabel('x')
ylabel('y')
zlabel('V(x,y)')
colorbar;
caxis([min(min(V2D)),  max(max(V2D))]);