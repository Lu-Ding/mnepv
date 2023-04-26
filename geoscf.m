% Illustration of geometry of SCF for mNEPv from numerical radius computation 
%
% (Ding.Lu @ uky.edu, dated 04-20-2023)
%

close all; clear;
warning('off', 'map:removing:combntns');
set(0,'defaultTextInterpreter','latex'); 
rng(0);

% Set data matrices
n = 10;
B1 = orth(randn(n,n)); B1 = B1 *diag(randn(n,1))* B1';
B2 = orth(randn(n,n)); B2 = B2 *diag(randn(n,1))* B2';
a = 1; b = 1.2;
B1 = B1*a - 1*eye(n); 
B2 = B2*b + .8*eye(n);

% Plot the numerical range
figure
numrange(B1+1i*B2, [0,0,0] + 0.9); hold on 

% Coefficient functions
phi = @(x,y) (x).^2 + (y).^2;
nx = @(x,y) [2*x; 2*y];

% Starting point, manually picked on plots
x0= [ -.7; -.15];
nx0 = nx(x0(1),x0(2)); nx0 = nx0/norm(nx0)*0.5;
quiver(x0(1), x0(2), nx0(1), nx0(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1, 'Marker', 'o', 'MarkerFaceColor', 'k' );
text(x0(1)+0.2,x0(2)-.1, '$g(x_0)$', 'Fontsize', 16)

% Draw first tangent line
H = B1*nx0(1) + B2*nx0(2);
[v,e] = eigs(H,1,'largestreal');
x1 = [real(v'*B1*v); real(v'*B2*v)]; % touching point
t = linspace(-1.8, 1.8, 2);
tg = x1 + [0,1; -1,0]*nx0*t;
plot(tg(1,:), tg(2,:), 'k', 'Linewidth', 1.5)

% Draw second point
nx1 = nx(x1(1),x1(2)); nx1 = nx1/norm(nx1)*0.5;
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1, 'Marker', 'o', 'MarkerFaceColor', 'k');
text(x1(1)+0.2,x1(2)-0.1, '$g(x_1)$', 'Fontsize', 16)


% Draw second tangent line
H = B1*nx1(1) + B2*nx1(2);
[v,e] = eigs(H,1,'largestreal');
x2 = [real(v'*B1*v); real(v'*B2*v)]; % touching point
tg = x2 + [0,1; -1,0]*nx1*t;
plot(tg(1,:), tg(2,:), 'k', 'Linewidth', 1.5)

% Draw 3rd point
nx2 = nx(x2(1),x2(2)); nx2 = nx2/norm(nx2)*0.5;
quiver(x2(1), x2(2), nx1(1), nx1(2), '-k', 'Linewidth',1, 'MaxHeadSize', 1, 'Marker', 'o', 'MarkerFaceColor', 'k');
quiver(x2(1), x2(2), nx2(1), nx2(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1, 'Marker', 'o', 'MarkerFaceColor', 'k');
text(x2(1)+0.1, x2(2)-0.1, '$g(x_2)$', 'Fontsize', 16)


% Draw global max point with tangent line
xmx = [-1.1, 2.75];
plot(xmx(1), xmx(2), 'Marker', 'p', 'Markersize', 15, 'MarkerFaceColor', 'k');
text(xmx(1)+0.1, xmx(2)+0.4, '$g(x_*)$', 'Fontsize', 16)
% 
nxx = nx(xmx(1),xmx(2)); nxx = nxx/norm(nxx)*0.5; nxx2 = nxx/norm(nxx)*0.7;
quiver(xmx(1), xmx(2), nxx(1), nxx(2), '-k', 'Linewidth',1, 'MaxHeadSize', 1, 'Marker', 'o', 'MarkerFaceColor', 'k');
quiver(xmx(1), xmx(2), nxx2(1), nxx2(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1, 'Marker', 'o', 'MarkerFaceColor', 'k');
tg = xmx' + [0,1; -1,0]*nxx*t;
plot(tg(1,:), tg(2,:), 'k', 'Linewidth', 1.5)


xlabel('$y(1)$', 'Fontsize', 15)
ylabel('$y(2)$', 'Fontsize', 15, 'rotation', 0, 'HorizontalAlignment','right')


% Draw contours
x = linspace(-3.1,1.4, 100);
y = linspace(-.8,3.6, 100);
[xx, yy] = meshgrid(x,y);
%lv = [23, 12, 6, 3];
lv = [phi(x0(1),x0(2)), phi(x1(1), x1(2)), phi(x2(1),x2(2)), phi(xmx(1),xmx(2))];
contour(xx, yy, phi(xx,yy), lv, '--k')
axis equal;
axis([-3.1,1.4,-.8,3.6]);

% Origin
plot(0,0,'+k','linewidth',2,'MarkerSize',10);
