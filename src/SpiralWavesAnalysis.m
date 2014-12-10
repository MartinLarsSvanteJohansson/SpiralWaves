% SpiralWavesAnalysis.m
% v01
%
% Note: This script requires the Symbolic Math Toolbox
%
% Identifies the homogenous steady states and check their linear 
% stability to homogenous perturbations
%
% Martin Johansson, 2011
%

clear all, clc, close all

syms a b

% Parameters
epsilon = 1/12;
a = 0.75;
b = 0.06;
L = 80+2;
d = 0;

u = -0.2:0.01:1.1;
v = -0.2:0.01:1.1;

 f = @(u,v) epsilon.^(-1).*u.*(1-u).*(u-(v+b)./a);
 g = @(u,v) u.^3-v;

% Nullclines
u1 = 0; u1 = linspace(u1,u1,length(u));
u2 = 1; u2 = linspace(u2,u2,length(u));
u3 = (v+b)./a;
v1 = u.^3;
u_nullclines = [u1; u2; u3];

plot(u_nullclines,v,u,v1)
xlabel('u','FontWeight','bold','FontSize',16)
ylabel('v','FontWeight','bold','FontSize',16)
legend('u = 0', 'u = 1', 'u = (v+b)/a', 'v = u^3', 'location', 'SouthEast')
axis image

% Finding the intersection between nullclines
v_ = meshgrid(0:0.0001:0.5,1);
func_ = (v_+b)./a - (v_).^(1/3);
[Y,I_1] = min(abs(func_));

v__ = meshgrid(0.5:0.0001:1,1);
func_ = (v__+b)./a - (v__).^(1/3);
[Y,I_2] = min(abs(func_));

sprintf('v1 = %g, u1 = %g; v2 = %g , u2 = %g ',v_(I_1), (v_(I_1)+b)/a,...
    v__(I_2),(v__(I_2)^(1/3)))

hold on
u = -0.1:0.1:1.1;
v = -.1:0.1:1.1;

f_vec = zeros(length(u),length(v));
g_vec = zeros(length(u),length(v));
for i = 1:length(v)
    f_vec(i,:) = f(u,v(i));
    g_vec(i,:) = g(u,v(i));
end

quiver(u,v,f_vec,g_vec)

% Calculating Jacobian

syms x y a b

f_ = epsilon.^(-1).*x.*(1-x).*(x-(y+b)./a);
g_ = x.^3-y;

df_u = diff(f_,x);
df_v = diff(f_,y);
dg_u = diff(g_,x);
dg_v = diff(g_,y);

% Our jacobian
J = [df_u df_v ; dg_u dg_v];

% Check the four found nullclines intersection points
x = 0; y = 0;
lambda_1 = real(eig(eval(J)))
x=1; y=1;
lambda_2 = real(eig(eval(J)))
x = 0.0806667; y = 0.0005;
lambda_3 = real(eig(eval(J)))
x = 0.822832 ; y = 0.5571;
lambda_4 = real(eig(eval(J)))
