% Example 3: mNEPv from distance problems of quadratic dHDAE systems
% 
% (Ding.Lu @ uky.edu, dated 04-20-2023)
%

close all; clear;
rng(0);

maxit = 100000;
tol = 1.0E-13;
tolnt = 0.1;

% Set true to run manopt comparison, requiring manopt installed.
doManopt = false;

% Set testing problem sizes
NN = [100,500,1000,2000,3000];
nx0 = 50; % number of initial points = 2 * nx0 
Ns = length(NN);
m = 4; 

%   
T1 = zeros(2*nx0, Ns); 	% SCF timing
T2 = T1; 				% SCF Acc timing
T3 = T1; 				% Manopt timing
ITS1 = T1; ITS2 = T1; ITS3 = T1; % iterations
FUN1 = T1; FUN2 = T1; FUN3 = T1; % objfun

% Start the main loop
for ii = 1:Ns

	% Define testing matrices of size n
    n = NN(ii)
    mm = [1:n]; dd = mm; kk =mm; 
    M=zeros(n); K=zeros(n); D=zeros(n);
    for j=1:n-1
        M(j,j)=mm(j);
        K(j,j)=kk(j)+kk(j+1);
        K(j,j+1)=-kk(j+1);
        K(j+1,j)=K(j,j+1);
        D(j,j)=dd(j)+dd(j+1);
        D(j,j+1)=-dd(j+1);
        D(j+1,j)=D(j,j+1);
    end
    M(n,n)=mm(n); K(n,n)=kk(n); D(n,n)=dd(n);
    B1 = M/n; B2=D/n; B3=K/n; 
    G = randn(n);  G = (G-G'); G=G/norm(G);
    
	% Define NEPv 
    AA = (G*G - B1*B1 - B2*B2 -B3*B3);
    A{1} = AA; A{2} = B1; A{3} = B2; A{4} = B3;
    fun{1} = @(y) y(1) + (y(2)^2 + y(3)^2 + y(4)^2)/2;
    fun{2} = @(y) [1; y(2); y(3); y(4)];
    fun{3} = @(y) [0; 1; 1; 1];
    
    H2 = @(d) d(1)*AA + d(2)*B1 + d(3)*B2 + d(4)*B3;
	Hx = @(x) H2( fun{2}([x'*AA*x, x'*B1*x, x'*B2*x, x'*B3*x]) );
	objfun = @(x) (x'*AA*x+ (x'*B1*x)^2/2+ (x'*B2*x)^2/2+ (x'*B3*x)^2/2);
	gradfun= @(x) 2*AA*x + 2*(x'*B1*x)*B1*x + 2*(x'*B2*x)*B2*x + 2*(x'*B3*x)*B3*x;
	hessfun= @(x,y) 2*AA*y + 2*(x'*B1*x)*B1*y + 2*(x'*B2*x)*B2*y + 2*(x'*B3*x)*B3*y + 4*(x'*B1*y)*B1*x + 4*(x'*B2*y)*B2*x + 4*(x'*B3*y)*B3*x;
    
    % Generate initial vectors using supporting points for the numerical
	% range, along nx0 random directions
    V0 = randn(m,nx0);
    X = []; 
    for i = 1:size(V0,2)
        v0 = V0(:,i);
        HH = H2(v0);
        [VV,EE] = eig(HH);
        % get the max end
        [ee,idx] = max(diag(EE));
        x0 = VV(:,idx);
        x0 = x0/norm(x0);
        X = [X, x0];
        % get the min end
        [ee,idx] = min(diag(EE));
        x0 = VV(:,idx);
        X = [X, x0];
    end

    % Run SCF w/wo refinement with different initial vector in X
    kk = size(X,2);
    for i = 1:kk
        x00 = X(:,i); 

    	% 1. Run SCF without refinement
        x0 = x00/norm(x00);
        tic;
        [x0, fx0, its, hist] = scf(A, fun, x0, 0, maxit, tol);
        t1 = toc; 
		T1(i,ii) = t1; ITS1(i,ii) = its; FUN1(i,ii) = fx0;

    	% 2. Run SCF with refinement
        x0 = x00/norm(x00);
        tic;
        [x0, fx02, its2, hist2] = scf(A, fun, x0, 1, maxit, tol, tolnt);
        t2 = toc;
	    T2(i,ii) = t2; ITS2(i,ii) = its2; FUN2(i,ii) = fx02; 

		% 3. Run Manopt comparison, if required
		if doManopt
			manifold =  spherefactory(n); 
			problem.M = manifold;
			problem.cost  = @(x) -1*objfun(x);
			problem.egrad = @(x) -1*gradfun(x); 
			problem.ehess = @(x, y) -1*hessfun(x,y);
			%checkgradient(problem); 	%checked
			%checkhessian(problem); 	%checked
			options.verbosity = 0;
			options.tolgradnorm = tol * (1+norm(Hx(x0), 1)); % Note NEPv rel_resd = grad * norm(H(x_*),1)
			x0 = x00/norm(x00);

			tic;
			[x0, xcost, info, options] = trustregions(problem,x0,options);
			t3 = toc;

        	T3(i,ii) = t3; ITS3(i,ii) = length(info); FUN3(i,ii) = -xcost; 

			disp(['Timing SCF = ',num2str(t1), ' Accl.= ',num2str(t2), ' Manopt= ',num2str(t3)]);
        	disp(['fval/fscf: ', num2str([fx0,fx02,-xcost]/fx0)]);

		else
			disp(['Timing SCF = ',num2str(t1), ' Accl.= ',num2str(t2)]);
        	disp(['fval/fscf: ', num2str([fx0,fx02]/fx0)]);

		end

    end

end


disp(['-------'])
disp(['--Timing for n = ', num2str(NN)]);
disp(['- SCF (mean): ', num2str(mean(T1))]);
disp(['- SCF (div) : ', num2str(max(abs(T1-mean(T1))))]);
disp(['- SCFacc (mean): ', num2str(mean(T2))]);
disp(['- SCFacc (div) : ', num2str(max(abs(T2-mean(T2))))]);
if doManopt
	disp(['- manopt (mean): ', num2str(mean(T3))]);
	disp(['- manopt (div) : ', num2str(max(abs(T3-mean(T3))))]);
end
