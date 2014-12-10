% SpiralWaves.m 
% v01
%
%
% 
% Martin Johansson, 2011
%

clear all, clc, close all

% - Parameters
epsilon = 1/14;
a = 0.75;
b = 0.06;
L = (80+2);
d = 0;
add_extra = true; % For enjoyment

% Initializations
u_ = -0.1:0.01:1.1;
v_ = -0.1:0.01:1.1;

u = zeros(L,L);
v = zeros(L,L);
u_new = zeros(L,L);

% Calculating nullclines
u1 = 0; u1 = linspace(u1,u1,length(u_));
u2 = 1; u2 = linspace(u2,u2,length(u_));
u3 = (v_+b)./a;
v1 = u_.^3;
u_nullclines = [u1; u2; u3];

% Storing the f(u,v) and g(u,v)
f = @(u,v) epsilon.^(-1).*u.*(1-u).*(u-(v+b)./a);
g = @(u,v) u.^3-v;

% Time step size
dt = 0.05;
% Time vector
t= 0:dt:70;

% Initial position and concentrations for u and v
u(2:60,30:34)= 0.9;
v(2:60,32:36) = 0.8;

f1 = figure(1)
       subplot(2,1,1)
       imagesc(u(2:L-1,2:L-1))
       title('Initial concentration of u','FontSize', 11, 'FontWeight', 'bold')
       subplot(2,1,2)        
       imagesc(v(2:L-1,2:L-1))
       title('Initial concentration of v','FontSize', 11, 'FontWeight', 'bold')
movegui(f1,'northwest')

% Starting simulation
for i = 1:length(t)

   % Setting the zero-flux boundaries
   u(1,1:L) = u(2,1:L);     
   u(L,1:L) = u(L-1,1:L);
   u(1:L,L) = u(1:L,L-1);   
   u(1:L,1) = u(1:L,2);
     
   % Updating u and v with Euler forward and laplacian discretization
   u_new = u + dt.*(f(u,v) + 4.*del2(u));
   v = v + dt.*g(u,v);  
   u = u_new;
   
   % Plotting every tenth time step
   if (0 == mod(i,10) )       
       hold on
       
       % Plotting the u and v concentration
       f2 = figure(2);
       hold on
       colormap();
       subplot(2,1,1);
            imagesc(u(2:L-1,2:L-1))
            annotation(figure(2),'ellipse',...
                [0.282 0.82 0.02 0.02],'FaceColor',[1 0 0],...
                'Color',[1 0 0]);
       subplot(2,1,2);
            imagesc(v(2:L-1,2:L-1))
            annotation(figure(2),'ellipse',...
                [0.282 0.348 0.02 0.02],'FaceColor',[1 0 0],...
                'Color',[1 0 0]);
            
       set(f2, 'Position', [400 400 560 420])
       pause(0.06)
       
       % Concentration behaviour over time for a fixed point in space
       f3 = figure(3);
       plot(u(20,20),v(20,20),'o')
       axis([-0.1 1.1 -0.1 1.1])
       plot(u_nullclines,v_,u_,v1)
       set(f3, 'Position', [1000 400 560 420])
       hold off        
   end
   
   % Adding some extra v at time 25 and 50 for enjoyment
   if ( 0 == mod(i,500) && add_extra )       
       v5(45:55,45:55) = 0.9;       
   end
   
   % Displaying elapsed time units
   clc; disp('Time: '); disp(t(i))      
end

% Adding labels and dots
figure(2)
subplot(2,1,1)
title('Concentration of u','FontSize', 11, 'FontWeight', 'bold')
xlabel('x','FontWeight','bold','FontSize',14)
ylabel('y','FontWeight','bold','FontSize',14)
colorbar()

subplot(2,1,2)
title('Concentration of v','FontSize',11,'FontWeight','bold')
xlabel('x','FontWeight','bold','FontSize',14)
ylabel('y','FontWeight','bold','FontSize',14)
colorbar

figure(3)
xlabel('u','FontWeight','bold','FontSize',16)
ylabel('v','FontWeight','bold','FontSize',16)

figure(4)
plot(2:L-1,u(2:L-1,20),'k',2:L-1,v(2:L-1,20),'r')
xlabel('y','FontWeight','bold','FontSize',16)
ylabel('Concentration','FontWeight','bold','FontSize',16)
legend('u', 'v')

% Clean up a bit
clear ans i u1 u2 u3 u_new f1 f2 f3
