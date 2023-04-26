% Example 4: mNEPv from Tensor rank-1 approximation.
%
% (Ding.Lu @ uky.edu, dated 04-20-2023)
%

close all; clear
warning('off', 'map:removing:combntns');
set(0,'defaultTextInterpreter','latex'); 

rng(0);
maxit = 10000;
tol = 1.0E-13;

% Select testing tensors:
% 	'n' = New Orleans; 'p' = Princeton; 'r' = Reuters; otherwise = random
% 	RNDSTART = number of random starting vectors; 

tensor_name= 'p'; 
RNDSTART = 100; 

if strfind(tensor_name, 'n')
    load ./data/fb20.mat;
    FNAME = 'New Orleans';

elseif strfind(tensor_name, 'p')
    load ./data/princeton.mat;
    FNAME = 'Princeton';

elseif strfind(tensor_name, 'r')
    load ./data/reuters.mat;
    FNAME = 'Reuters';

else % synthetic
    n = 100;
    m = 10;
    for i = 1:m
        B1 = orth(randn(n,n)); B1 = B1 *diag(randn(n,1))* B1';
        A{i} = B1;
    end
    FNAME = 'Synthetic';
end

% Size and number of matrices
n = size(A{1},1); 	
m = size(A,2);  	

% Objective functions, 1st and 2nd derivatives
fun{1} = @(y) norm(y,2)^2;
fun{2} = @(y) 2*(y);
fun{3} = @(y) 2*ones(m,1);

% Random starting vectors: collected in matrix X
X = abs(randn(n,RNDSTART));

% Run SCF with or without refinement from different initial vectors in X
for i = 1:size(X,2)

	% 0. Initial vectors
    x00 = X(:,i); 

    % 1. Run SCF without refinement
    tic
    x0 = x00/norm(x00);
    [x0, fx0, its0, hist0] = scf(A, fun, x0, 0, maxit, tol, 1.0E-1, 1, listA);
    t0 = toc
	resd0 = hist0.RESD(end);


    % 2. Run SCF with refinement
    tic;
    x1 = x00/norm(x00);
    [x1, fx1, its1, hist1] = scf(A, fun, x1, 1, maxit, tol, 1.0E-1, 1, listA);
    t1 = toc
	resd1 = hist1.RESD(end);
   
    disp('time: SCF / SCF acc')
    disp([t0, t1])
    disp('fval: SCF / SCF acc')
    disp([fx0, fx1])
    disp('Resd: SCF /SCF acc ')
    disp([resd0, resd1])

    % Record timing and history
    T0(i) = t0; T1(i) = t1; 
    ITS0(i) = its0; ITS1(i) = its1;
    FVAL0(i) = fx0; FVAL1(i) = fx1;
    HIST0{i} = hist0; HIST1{i} = hist1;

end

% Plot convergence history (at most 8 curves);
nlines = min(size(X,2), 8); 
corder = lines(nlines+3);
corder = flipud(corder);

for i = 1:nlines  

    hist0 = HIST0{i}; hist1 = HIST1{i}; 
    its0 = ITS0(i); its1 = ITS1(i);
    RESD0 = hist0.RESD; OBJFX0 = hist0.OBJFX; x0hist0 = hist0.Xhist; 
    RESD1 = hist1.RESD; OBJFX1 = hist1.OBJFX; x0hist1 = hist1.Xhist; 

    % Figure 1: convergence of objective function value
    figure(1);
    plot([0:its0-1], OBJFX0, '--o','color', corder(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
    plot([0:its1-1], OBJFX1, '--x','color', corder(i,:), 'Markersize', 10, 'linewidth', 1); hold on

    % Figure 2: convergence of relative residual norm 
    figure(2);
    semilogy([0:its0-1], RESD0,'--o','color', corder(i,:), 'Markersize', 10, 'linewidth', 1, 'DisplayName', 'SCF'); hold on
    semilogy([0:its1-1], RESD1,'--x','color', corder(i,:), 'Markersize', 10, 'linewidth', 1); hold on
end

% Set labels
figure(1); % Objective function values
xlabel('$k$')
ylabel('$F(x_k)$'); %, 'rotation', 0)
title(FNAME, 'interpreter', 'latex')  
%xlim([-0.5,10.5])

figure(2); % Relative residual norms
xlabel('$k$')
%ylabel('res$(x_k)$ ')
title(FNAME, 'interpreter', 'latex')  
ylim([1.0E-16,1.0E0])
xlim([-0.5, max(ITS0) + 0.5])

% Display timing statistics in figure(1)
dim = [.42, 0.6, .3, .3];
tt0 = mean(T0); dtt0 = max(abs(T0-mean(T0)));
tt1 = mean(T1); dtt1 = max(abs(T1-mean(T1)));
str = {'Timing in seconds:',...
	['\phantom{accl. }SCF ~~', num2str(tt0,'%.2f'),' ($\pm$', num2str(dtt0,'%.2f'),')'],...
	['accl. SCF ~~', num2str(tt1,'%.2f'),' ($\pm$', num2str(dtt1,'%.2f'), ')']};
annotation('textbox',dim,'interpreter', 'latex','String',str,'FitBoxToText','on', 'FontSize', 13, 'EdgeColor', 'k');


disp('average time (var): SCF / SCF acc / RTR')
[mean(T0), max(abs(T0-mean(T0)));
 mean(T1), max(abs(T1-mean(T1)))]'

disp('average iters (var): SCF / SCF acc / RTR')
[mean(ITS0), max(abs(ITS0-mean(ITS0)));
 mean(ITS1), max(abs(ITS1-mean(ITS1)))]'

disp('average fval (var): SCF / SCF acc / RTR')
[mean(FVAL0), max(abs(FVAL0-mean(FVAL0)));
 mean(FVAL1), max(abs(FVAL1-mean(FVAL1)))]'

return

% END OF EXAMPLE 2
