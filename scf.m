function [xx, fxx, its, hist] = scf(A, fun, x0, doacc, maxit, tol, tolnt, iterSol, listA)
% SCF iteration with or without acceleration.
%
% INPUT:
% 	A:			cell containing m Hermitian matrices
% 	fx:       	fx{1}=function, fx{2} = derivative, fx{3} = 2nd derivative
% 	x0:       	starting vector
% 	doacc:    	apply acceleration if true
%
% (optional)
% 	maxit:		maximum iteration 
% 	tol:		residual tolerance  
% 	tolnt: 		local acceleration threshold
% 	iterSol:	apply iterative linear/eigen solver if true
%	listA: 		struct with fields vA, ri, ci, n. This is used if all
%				matrices Aj share the same sparsity pattern S, with
%				row/column indices ri/ci. Entries at S of the j-th
%				matrix Aj are saved as a column vector vA{j} = vec(Aj(S)).
%
% OUTPUT:
% 	xx:       	solution eigenvector
% 	fxx:      	objective function value F(xx)
% 	its:      	iteration number
% 	hist:     	history file
%
% (Ding.Lu @ uky.edu, dated 04-20-2023)
%

% Set default parameters
if nargin<5, maxit = 30; end
if nargin<6, tol = 1.0E-13; end
if nargin<7, tolnt = 1.0E-1; end
if nargin<8, iterSol = false; end
if nargin<9, listA=[]; end
warning("off", 'MATLAB:nearlySingularMatrix');

% Set local parameters and function handles 
m = length(A); n = length(A{1});
fx = fun{1}; dfx = fun{2}; ddfx = fun{3};
rhofun = @(x) rhofun0(x, A, m);
if ~iterSol || isempty(listA)
	H = @(x) Hfun0(dfx(rhofun(x)), A, m);
else
	H = @(x) Hfun0s(dfx(rhofun(x)), listA, m);
end
XK = @(x) XKfun(x, A, m); 

% -------------------------
% -- Begin the main loop --
% -------------------------
OBJFX = []; RESD = []; Xhist = []; Rhohist = [];

for j=1:maxit
    % Set matrix for the new iteration
    HH = H(x0);
    rho0 = rhofun(x0);
    fx0 = fx(rho0);

    % Check residual
    ee = real(dot(x0,HH*x0));
    rk = HH*x0 - ee*x0;
    resd = norm(rk)/norm(HH,1);

	% uncomment to change shift every 6 iterations
    %if mod(j,6)==1, eee = ee; end 
    
    % Run local acceleration if needed 
    if doacc && resd < tolnt && j > 1
        Xk = XK(x0); 
        Xk = Xk - x0*(x0'*Xk);
		% Apply direct or iterative solvers
        if ~iterSol 
            Wk = Xk*diag(2*ddfx(rho0))*Xk';
            HHk = (HH - ee*eye(size(HH)) + Wk);
            xx = HHk\x0;
        else
            CKk = ddfx(rho0); 
            AFUN = @(x) HH*x - ee*x + Xk*(2*CKk.*(x'*Xk)');
            xx = minres(AFUN, x0, min(resd^2,1.0E-3), 2000);

			% uncomment to use other linear solvers
            %xx = gmres(AFUN, x0, 10, min(resd^2,1.0E-3), 2000);
            %xx = symmlq(AFUN, x0, min(resd^2,1.0E-3), 2000);
        end

        xx = xx/norm(xx);
        fk = fx(rhofun(xx));

		% Accept or reject updates by checking monotonicity of objective values
        if fk >= fx0 - (abs(fx0)+1)*10*eps
            x0 = xx;
            fx0 = fk;
            HH = H(x0);
            rho0 = rhofun(x0);
            ee = real(dot(x0,HH*x0));
            rk = HH*x0 - ee*x0;
            resd = norm(rk)/norm(HH,1);
        end
    end

	% Record computation histories 
    RESD = [RESD, resd];
    OBJFX = [OBJFX, fx(rho0)];
    Xhist = [Xhist, x0];
    Rhohist = [Rhohist, rho0];

	% Check stopping criteria 
    if resd < tol, break; end

    % Run SCF iteration, using direct or iterative solver
    if ~iterSol
        [VV,EE] = eig(HH);
    else
        [VV,EE] = eigs(HH, 1, 'largestreal', 'StartVector', x0, 'Tolerance', min([1.0E-3,resd^2]), 'SubspaceDimension', min([10,n]));
    end
    [~,idx] = max(real(diag(EE)));
    x0 = VV(:,idx); 
    x0 = x0/norm(x0);
end

% Prepare to return
xx = x0; fxx=fx0;
its = j;
hist.RESD = RESD;
hist.OBJFX = OBJFX;
hist.Xhist = Xhist;
hist.Rhohist = Rhohist;

% END OF SCF
end


% ----------------------------
% -- Define local functions --
% ----------------------------

% Routine to evaluate Xk = [A1*x, A2*x, ..., Am*x]
function Xk = XKfun(x, A, m)
Xk = [];
for i = 1:m
    Xk = [Xk, A{i}*x];
end
end

% Routine to evaluate Rayleigh quotients y = [x'*A1*x, x'*A2*x, ..., x'*Am*x]
function y = rhofun0(x, A, m)
y = zeros(m,1);
for i = 1:m
    y(i) = real(dot(x,A{i}*x));
end
% a one-line equivalent: y = cellfun(@(A) real(x'*(A*x)), A);
end

% Routine to evaluate H(y)
function H = Hfun0(y, A, m)
H = y(1)*A{1};
for i = 2:m
    H = H + y(i)*A{i};
end
end

% Routine to evaluate H(y), for A{i} saved in vector forms
function H = Hfun0s(y, listA, m)
vA = listA.vA;
H = y(1)*vA{1};
for i = 2:m
    H = H + y(i)*vA{i};
end
H = sparse(listA.ri, listA.ci, H, listA.n, listA.n);
end
