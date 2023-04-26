% Example 1: mNEPv from numerical radius computation  
% 
% (Ding.Lu @ uky.edu, dated 04-20-2023)
%

close all; clear all;
warning('off', 'map:removing:combntns');
set(0,'defaultTextInterpreter','latex'); 

rng(0);
maxit = 100000;
ntheta = 100; % sampled theta

% Set testing matrices: a randomly generated 4-by-4 example
n = 4;
B = [0.6 + 0.6i  -0.2 + 2.5i  -1.9 - 0.2i  -0.3 + 2.5i;
  	-0.1 + 2.3i  -0.3 - 2.6i  -1.3 + 0.4i  -1.2 + 1.3i;
  	-2.0 + 0.0i  -1.6 + 0.6i  -2.1 - 0.4i   1.3 + 1.2i;
  	-0.1 + 2.0i  -1.6 + 1.4i   1.5 + 1.0i  -0.1 - 2.3i ]; 
B1 = (B+B')/2; B2 = (B-B')/(2i);
A{1} = B1; A{2} = B2;
fun{1} = @(y) y(1)^2 + y(2)^2;	%fx
fun{2} = @(y) [2*y(1); 2*y(2)]; %dfx
fun{3} = @(y) [2;  2];			%ddfx

% First run of SCF at 4 sampled boundary points along x and y axis 
THETA = [0, pi/2, pi, 3/2*pi];
ncolor0 = length(THETA);
corder0 = lines(ncolor0+3); % specified color order
corder0 = flipud(corder0);

for i = 1:length(THETA)

	% 0. Sample initial vector along angle Theta(i)
    t = THETA(i);
    HH = cos(t)*B1 + sin(t)*B2;
    [VV,EE] = eigs(HH, 1, 'largestreal');
    x00 = VV;
    
    % 1. Run SCF without refinement
    x0 = x00/norm(x00);
    [x0, fx0, its, hist] = scf(A, fun, x0, 0, maxit);
    x00hist{i} = x0; 
	fxhist(i) = fx0; 

    % 2. Run SCF with refinement
    x0 = x00/norm(x00);
    [x0, fx02, its2, hist2] = scf(A, fun, x0, 1, maxit);

    % 3. Plot history
    RESD = hist.RESD; 
	OBJFX = hist.OBJFX; 
    RESD2 = hist2.RESD; 
	OBJFX2 = hist2.OBJFX; 
    
    % 3.1. Figure 2: convergence of objective function value
    figure(2);
    plot([0:its-1], OBJFX, '--o','color', corder0(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
    plot([0:its2-1], OBJFX2, '--x','color', corder0(i,:), 'Markersize', 10, 'linewidth', 1); hold on

    % 3.2. Figure 3: convergence of residual norm 
    figure(3);
    semilogy([0:its-1], RESD,'--o','color', corder0(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
    semilogy([0:its2-1],RESD2,'--x','color', corder0(i,:), 'Markersize', 10, 'linewidth', 1); hold on
end

figure(2)
xlabel('$k$')
ylabel('$F(x_k)$'); %, 'rotation', 0)
title('function value', 'interpreter', 'latex')  
xlim([-0.5,10.5])

figure(3)
xlabel('$k$')
ylabel('res$(x_k)$ ')
xlim([-0.5,15.5])
title('relative residual norm', 'interpreter', 'latex')  

% Figure 1: the numerical range
figure(1); 
numrange(B1+1i*B2, [0,0,0] + 0.9); hold on
plot( 0, 0, '+k', 'Markersize', 12, 'linewidth', 2)
   
% Second run of SCF at sampled points along the boundary of numerical range		
THETA = linspace(0,2*pi,ntheta);
corder = lines(ntheta+3); % specified color order
corder = flipud(corder);

for i = 1:length(THETA)

	% 0. Sample initial vector along angle Theta(i)
    t = THETA(i);
    HH = cos(t)*B1 + sin(t)*B2;
    [VV,EE] = eigs(HH, 1, 'largestreal');
    x00 = VV;
    
    % 1. Run SCF without refinement
    x0 = x00/norm(x00);
    [x0, fx0, its, hist] = scf(A, fun, x0, 0, maxit);

    % 2. Run SCF with refinement
    x0 = x00/norm(x00);
    [x0, fx02, its2, hist2] = scf(A, fun, x0, 1, maxit);

    % x. Record history
    its_scf(i) = its;
    its_acc(i) = its2;
    fun_scf(i) = fx0;
    fun_acc(i) = fx02;

    % 3. Figure 1: solution in the numerical range
    figure(1);
    ic = find(abs(fxhist-fx0)<1.0E-6, 1, 'first'); % find corresponding color 
    x0 = hist2.Xhist(:,1); % initial point
    plot( x0'*B1*x0,  x0'*B2*x0, 'o', 'Markersize', 6, 'color', corder0(ic,:), 'linewidth', 2, 'MarkerFacecolor', corder0(ic,:))

% %
% % UNCOMMENT IF YOU WANT TO TRACK CONVERGENCE HISTORY
% %
%    RESD = hist.RESD; OBJFX = hist.OBJFX; x0hist = hist.Xhist; 
%    RESD2 = hist2.RESD; OBJFX2 = hist2.OBJFX; x0hist2 = hist2.Xhist; 
%    % 3.2. Figure 5: convergence of objective function value
%    figure(6);
%    plot([0:its-1], OBJFX, '--o','color', corder(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
%    plot([0:its2-1], OBJFX2, '--x','color', corder(i,:), 'Markersize', 10, 'linewidth', 1); hold on
%
%    % 3.3. Figure 6: convergence of residual norm 
%    figure(7);
%    semilogy([0:its-1], RESD,'--o','color', corder(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
%    semilogy([0:its2-1],RESD2,'--x','color', corder(i,:), 'Markersize', 10, 'linewidth', 1); hold on

end


% Figure 1: solution in the numerical range
figure(1)
for i = 1:ncolor0
    x0 = x00hist{i};    
    rd = norm([x0'*B1*x0,  x0'*B2*x0]); 
    plot(rd*exp(1i*linspace(0,2*pi,100)), '--k');
    plot( x0'*B1*x0,  x0'*B2*x0, 'pk', 'Markersize', 12, 'linewidth', 1, 'MarkerFacecolor','w')
end
xlabel('$y(1)$', 'Fontsize', 16)
ylabel('$y(2)$', 'Fontsize', 16, 'rotation', 0)
axis([-5,5,-5,5])
text(-4, 3.5, 'III', 'Fontsize', 16)
text(-4, -4, 'I', 'Fontsize', 16)
text(2, -4.5, 'II', 'Fontsize', 16)

% Figure 4: summary of iteration numbers 
figure(4);
dff = diff(fun_scf);
idxx = find(abs(dff) > 1.0E-2, 3, 'first'); % find breaking point
plot(THETA, its_acc, 'xb','MarkerSize',8,'linewidth',1,'markerfacecolor','b'); hold on;
plot(THETA, its_scf, 'or','MarkerSize',6,'linewidth',1,'markerfacecolor','w'); hold on;
for i = 1:length(idxx)
    theta = THETA(idxx(i));
    plot([theta,theta], [0,max(its_scf)*1.05],':k');
end
xlim([0,6.4])
ylim([0,max(its_scf)*1.05])
xlabel('$\theta$', 'Fontsize', 16);
ylabel('iteration number', 'Fontsize', 16);
legend('plain SCF','accelerated','Location','northeast')
text(.1, 120, 'II', 'Fontsize', 16)
text(1.5, 120, 'III', 'Fontsize', 16)
text(3.8, 120, 'I', 'Fontsize', 16)
text(5.5, 120, 'II', 'Fontsize', 16)

% Figure 5: summary of objective values at solution 
figure(5)
plot(THETA, fun_acc, 'xb','MarkerSize',8,'linewidth',1,'markerfacecolor','b'); hold on;
plot(THETA, fun_scf, 'or','MarkerSize',6,'linewidth',1,'markerfacecolor','w'); hold on;
for i = 1:length(idxx)
    theta = THETA(idxx(i));
    plot([theta,theta], [0,max(fun_scf)*1.1],':k');
end
xlim([0,6.4])
ylim([min(fun_scf)*0.95,max(fun_scf)*1.05])
xlabel('$\theta$', 'Fontsize', 16);
ylabel('$F(x_*)$', 'Fontsize', 16);
legend('plain SCF','accelerated','Location','northeast')
