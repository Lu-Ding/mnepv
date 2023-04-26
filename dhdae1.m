% Example 2: mNEPv from distance problems of linear dHDAE systems
% 
% (Ding.Lu @ uky.edu, dated 04-20-2023)
%

close all; clear;
warning('off', 'map:removing:combntns');
set(0,'defaultTextInterpreter','latex'); 
rng(0);

maxit = 10000;
tol = 1.0E-13;
tolnt = 0.1;


% Set testing problem 
m = 3; n = 30; 
B1 = orth(rand(n,n)); B1 = B1 *diag(rand(n,1)+1.0E-6)* B1';
B2 = orth(rand(n,n)); B2 = B2 *diag(rand(n,1)+1.0E-6)* B2';    
J = rand(n,n); J = (J-J'); J =J/norm(J); 
AA = (J*J - B1*B1 - B2*B2);

A{1} = AA; A{2} = B1; A{3} = B2; 
fun{1} = @(y) y(1) + (y(2)^2 + y(3)^2)/2;
fun{2} = @(y) [1; y(2); y(3)];
fun{3} = @(y) [0; 1; 1];
H2 = @(d) d(1)*AA + d(2)*B1 + d(3)*B2;


% Generate initial vectors using supporting points of the numerical
% range, along sampled directions over equally spaced angles 
THETA = linspace(0,2*pi, 40);
ETA = linspace(0, pi, 20);
XX=[]; YY=[]; ZZ=[];

for i = 1:length(THETA)

	theta = THETA(i);

	for j = 1:length(ETA)

		phi = ETA(j);
		v0 = [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)]; % normal dir
    	HH = H2(v0);
    	[x0,e0] = eigs(HH, 1, 'largestreal');
		V00{i,j} = x0;

		% Save boundary points 
		a = x0'*AA*x0; b = x0'*B1*x0; c = x0'*B2*x0;
		XX = [XX;a]; YY = [YY;b]; ZZ = [ZZ;c];

	end

end


% Run SCF w/wo refinement with different initial vector in X
its_scf = zeros(length(THETA), length(ETA));
its_acc =its_scf;
fun_scf =its_scf;
fun_acc =its_scf;

for i = 1:length(THETA)

	for j = 1:length(ETA)

    	x00 = V00{i,j};   

    	% 1. Run SCF without refinement
		x0 = x00/norm(x00);
		[x0, fx0, its, hist] = scf(A, fun, x0, 0, maxit);

    	% 2. Run SCF with refinement
		x0 = x00/norm(x00);
		[x0, fx02, its2, hist2] = scf(A, fun, x0, 1, maxit);

		% 
		its_scf(i,j) = its;
		its_acc(i,j) = its2;
		fun_scf(i,j) = fx0;
		fun_acc(i,j) = fx02;

	end

end


% ----------------------------------------------------------
% - Figure 1: Show numerical range and computed solution   -
% ----------------------------------------------------------
figure(1);
KK = boundary(XX,YY,ZZ);
trisurf(KK,XX,YY,ZZ, 'Facecolor',[.9,.9,.9],'FaceAlpha',0.9, 'EdgeColor', 'k','linewidth', 1); hold on;
plot3(XX, YY, ZZ, '.', 'Markersize', 6, 'color', [.1,.1,.1], 'linewidth', 1); hold on;
xlabel('$y(1)$', 'Fontsize', 16);
ylabel('$y(2)$',  'Fontsize', 16); 
zlabel('$y(3)$', 'rotation', 0, 'Fontsize', 16);

% Plot level surface  
fs = fun_scf(1);
fxfun = @(y,z) fs-(y.^2+z.^2)/2;
fyfun = @(y,z) y;
fzfun = @(y,z) z;
fsurf(fxfun,fyfun,fzfun, [0,.8,0,.8],'FaceAlpha',0.3,'EdgeColor',[.7,.7,.7],'FaceColor',[.8,.8,.8]);
ylim([0,1]); zlim([0,1]);
grid on; axis equal;

% Plot convergence history at solutions found 
% - Get all the different solutions computed (may be more than one)
ffx = sort(fun_scf(:)); 
dffx = diff(ffx);
idx = find(abs(dffx) > 1.0E-3); 
ffx = [ffx(idx), ffx(end)]; 

% - For each solution found, show convergence history of SCF from a
% 	particular x00(THETA/ETA)
corder0 = lines(length(ffx));
corder0 = flipud(corder0);
for jj = 1:length(ffx)

	% Get initial vector that convergences ffx
	[i,j] = ind2sub(size(V00), find(abs(fun_scf - ffx(jj))<1.0E-3, 1, 'first'));
	x00 = V00{i,j};   
	x0 = x00/norm(x00);

	% Run SCF 
	[x0, fx0, its, hist] = scf(A, fun, x0, 0, maxit);
	x0hist = hist.Xhist;

	% Record convergence history over the numerical range
	X0 = []; Y0 = []; Z0 = []; 
	for ii = 1:size(x0hist,2);
		x0 = x0hist(:,ii);
		a = x0'*AA*x0; b = x0'*B1*x0; c = x0'*B2*x0;
		X0 = [X0;a]; Y0 = [Y0;b]; Z0 = [Z0;c];
	end

	% Draw convergence history for the first k steps
	k = 5;
	plot3(X0(1:k), Y0(1:k), Z0(1:k),'--o', 'color',corder0(jj,:), 'Markersize',10, 'linewidth', 2, 'MarkerFaceColor', 'w' );
	% - starting point
	plot3(X0(1), Y0(1), Z0(1),'--o', 'color',corder0(jj,:), 'Markersize', 12, 'linewidth', 2, 'MarkerFacecolor',corder0(jj,:));
	% - end point
	plot3(a, b, c,'p', 'color', corder0(jj,:), 'Markersize', 20, 'linewidth', 2, 'MarkerFacecolor','w');

end

% ----------------------------------------------------------
% - Figure 2: Show number of SCF iterations at grid points -
% ----------------------------------------------------------
figure(2);
[TM, PM] = meshgrid(THETA,ETA);
plot3(TM, PM, its_scf, 'or'); hold on;
plot3(TM, PM, its_acc, 'xb', 'MarkerSize',8);

xlabel('$\theta$', 'Fontsize', 16);
ylabel('$\eta$',  'Fontsize', 16); 
zlabel('iteration number', 'Fontsize', 16);
grid on; 


% ----------------------------------------------------------
% - Figure 3: Show objective values at grid points
% ----------------------------------------------------------
figure(3)
plot3(TM, PM, fun_scf, 'or'); hold on;
plot3(TM, PM, fun_acc, 'xb','MarkerSize',8);
xlabel('$\theta$', 'Fontsize', 16);
ylabel('$\eta$',  'Fontsize', 16); 
zlabel('$f(x_*)$', 'Fontsize', 16);
grid on; 


% Show convergence history from 6 starting vectors sampled along the axis
V0 = [eye(3), -eye(3)]; % 6 directions (column)  along axis 
corder = lines(size(V0,2)+3);
corder = flipud(corder);

for i = 1: size(V0, 2)

	% 0. Sample initial vector along direction V0(:,i)
	v0 = V0(:,i);
 	HH = H2(v0);
    [x00, ~ ] = eigs(HH, 1, 'largestreal');

    % 1. Run SCF without refinement
    x0 = x00/norm(x00);
    [x0, fx0, its, hist] = scf(A, fun, x0, 0, 30, tol);

    % 2. Run SCF with refinement
    x0 = x00/norm(x00);
    [x0, fx0, its2, hist2] = scf(A, fun, x0, 1, 30, tol, tolnt);

    % 3. Plot history
    RESD = hist.RESD; OBJFX = hist.OBJFX; x0hist = hist.Xhist; 
    RESD2 = hist2.RESD; OBJFX2 = hist2.OBJFX; x0hist2 = hist2.Xhist; 

    % 3.1. Figure 5: convergence of objective function value
    figure(5);
    plot([0:its-1], OBJFX, '--o','color', corder(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
    plot([0:its2-1], OBJFX2, '--x','color', corder(i,:), 'Markersize', 10, 'linewidth', 1); hold on

    % 3.2. Figure 6: convergence of residual norm 
    figure(6);
    semilogy([0:its-1], RESD,'--o','color', corder(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
    semilogy([0:its2-1],RESD2,'--x','color', corder(i,:), 'Markersize', 10, 'linewidth', 1); hold on

end


figure(5)
%title('function value', 'interpreter', 'latex')  
xlabel('$k$','Fontsize', 16);
ylabel('$F(x_k)$', 'Fontsize', 16); %, 'rotation', 0)
xlim([-0.5,6.5])
ylim([-1.5,0])

figure(6)
%title('relative residual norm', 'interpreter', 'latex')  
xlabel('$k$', 'Fontsize', 16);
ylabel('res$(x_k)$ ','Fontsize', 16);
xlim([-0.5,18.5])
ylim([1.0E-16,1.0E0])
