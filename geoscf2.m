% Illustration of stable and non-stable eigenvectors 
% 
% (Ding.Lu @ uky.edu, dated 04-20-2023)
%

close all; clear;
warning('off', 'map:removing:combntns');
set(0,'defaulttextInterpreter','latex') 
rng(1);
n = 30;


% Set data matrices
B = gallery('grcar', 30); 
B = (B+B')/2*.8 + (B-B')/2;
B = exp(1i*pi/4)*B; % rotate a bit for illustration
B1 = (B+B')/2;
B2 = (B-B')/(2i);


% Coefficient functions
nx = @(x) 2*x;
H2 = @(y) y(1)*B1 + y(2)*B2;


% Run SCF to find the global maximizer 
A{1} = B1; A{2} = B2; 
fun{1} = @(y) y(1)^2+y(2)^2;
fun{2} = @(y) 2*[y(1); y(2)];
fun{3} = @(y) [2; 2];
x00 = randn(n,1)+1i*randn(n,1); 
x0 = x00/norm(x00);
x0 = scf(A, fun, x0, 0, 300);


% Plot numerical range and the largest contour
figure(1); 
numrange(B1+1i*B2, [0,0,0] + 0.9); hold on              % - numrange
plot( 0, 0, '+k', 'Markersize', 12, 'linewidth', 2)
rd = norm([x0'*B1*x0,  x0'*B2*x0]);                     % - largest contour
plot(rd*exp(1i*linspace(0,2*pi,100)), '--k');

plot(rd*.6*exp(1i*linspace(0,2*pi,100)), '--k');
plot(rd*.3*exp(1i*linspace(0,2*pi,100)), '--k');

xlabel('$y(1)$', 'Fontsize', 15)
ylabel('$y(2)$', 'Fontsize', 15, 'rotation', 0, 'HorizontalAlignment','right')


% Mark solution point I 
x1 = [x0'*B1*x0; x0'*B2*x0];
nx0 = nx(x1); nx0 = nx0/norm(nx0);
nx1 = nx0*1.5;
plot(x1(1),x1(2), 'Marker', 'p', 'Markersize', 15, 'MarkerFaceColor', 'k');
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);

d = [1; -1]/sqrt(2);
zz = [x0'*B1*x0; x0'*B2*x0];
zz = zz - 2*(zz'*d)*d;
plot(zz(1), zz(2), '-k', 'Marker', 'p', 'Markersize', 15,  'MarkerFaceColor', 'k');

znx0 = nx0; znx0 = znx0 - 2*(znx0'*d)*d;
znx1 = nx1; znx1 = znx1 - 2*(znx1'*d)*d;
quiver(zz(1), zz(2), znx0(1), znx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(zz(1), zz(2), znx1(1), znx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);


% Mark arrows nearby point 1
x1 = [1; -6];
nx0 = x1/norm(x1);
HH = H2(x1);
[xx,~] = eigs(HH, 1, 'lr');
x1 = [xx'*B1*xx; xx'*B2*xx];
nx1 = nx(x1); nx1 = nx1/norm(nx1);
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);

zx1 = x1; zx1 = zx1 - 2*(zx1'*d)*d;
znx0 = nx0; znx0 = znx0 - 2*(znx0'*d)*d;
znx1 = nx1; znx1 = znx1 - 2*(znx1'*d)*d;
quiver(zx1(1), zx1(2), znx0(1), znx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(zx1(1), zx1(2), znx1(1), znx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);


% Mark arrows nearby point 2
x1 = [1; 0];
nx0 = x1/norm(x1);
HH = H2(x1);
[xx,~] = eigs(HH, 1, 'lr');
x1 = [xx'*B1*xx; xx'*B2*xx];
nx1 = nx(x1); nx1 = nx1/norm(nx1);
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);

zx1 = x1; zx1 = zx1 - 2*(zx1'*d)*d;
znx0 = nx0; znx0 = znx0 - 2*(znx0'*d)*d;
znx1 = nx1; znx1 = znx1 - 2*(znx1'*d)*d;
quiver(zx1(1), zx1(2), znx0(1), znx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(zx1(1), zx1(2), znx1(1), znx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);


% Mark arrows nearby point 3
HH = H2([1,1]/2);
ee = eigs(HH, 2, 'bothendsreal');
x1 = [ee(1); ee(1)];
plot(x1(1), x1(1), '-k', 'Marker', 'p', 'Markersize', 15, 'Linewidth', 1);
nx0 = -1*[cos(pi/4); sin(pi/4)];
nx1 = nx(x1); nx1 = nx1/norm(nx1)*1.5;
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);

x1 = [ee(2); ee(2)];
plot(x1(1), x1(1), '-k', 'Marker', 'p', 'Markersize', 15, 'Linewidth', 1);
nx0 = [cos(pi/4); sin(pi/4)];
nx1 = nx(x1); nx1 = nx1/norm(nx1)*1.5;
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);


% Mark arrows nearby point 4
alpha = pi/4 + pi/10;
x1 = [cos(alpha); sin(alpha)];
nx0 = x1/norm(x1);
HH = H2(x1);
[xx,~] = eigs(HH, 1, 'lr');
x1 = [xx'*B1*xx; xx'*B2*xx];
nx1 = nx(x1); nx1 = nx1/norm(nx1);
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);

zx1 = x1; zx1 = zx1 - 2*(zx1'*d)*d;
znx0 = nx0; znx0 = znx0 - 2*(znx0'*d)*d;
znx1 = nx1; znx1 = znx1 - 2*(znx1'*d)*d;
quiver(zx1(1), zx1(2), znx0(1), znx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(zx1(1), zx1(2), znx1(1), znx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);


% Mark arrows nearby point 5
x1 = [-.905; .244]; % picked from figure
nx0 = -1*[cos(pi/4); sin(pi/4)];
nx1 = nx(x1); nx1 = nx1/norm(nx1);
quiver(x1(1), x1(2), nx0(1), nx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(x1(1), x1(2), nx1(1), nx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);

zx1 = x1; zx1 = zx1 - 2*(zx1'*d)*d;
znx0 = nx0; znx0 = znx0 - 2*(znx0'*d)*d;
znx1 = nx1; znx1 = znx1 - 2*(znx1'*d)*d;
quiver(zx1(1), zx1(2), znx0(1), znx0(2), '-k', 'Linewidth', 1, 'MaxHeadSize', 1);
quiver(zx1(1), zx1(2), znx1(1), znx1(2), ':k', 'Linewidth', 1, 'MaxHeadSize', 1);


axis([-4,4,-4,4]);
