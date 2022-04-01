
function [xopt,fopt,niter,gnorm,dx] = grad_descent(varargin)
% grad_descent.m mostra il metodo del gradiente per il problema di
% ottimizzazione non vincolata per la funzione f(x_1,x_2)=x_1^2+x_1x_2+3x_2^2.
% Ho usato come referenza il sito: nl.mathworks.com/matlabcentral/fileexchange.

%punto iniziale e tolleranza per la terminazione
x0 = [3 3]';
tol = 1e-6;

% massimo numero di iterazioni, minimo errore, passo del metodo
maxiter = 1000;
dxmin = 1e-6;
alpha = 0.2;

% inizializza il metodo del gradiente
gnorm = inf; x = x0; niter = 0; dx = inf;

% definizione della funzione obiettivo per il grafico
f = @(x1,x2) x1.^2 + x1.*x2 + 3*x2.^2;

% linee di livello della funzione obiettivo
figure(1); clf; ezcontour(f,[-5 5 -5 5]); axis equal; hold on

% ridefinizione della funzione obiettivo per l'ottimizzazione
f2 = @(x) f(x(1),x(2));

% implementazione del metodo del gradiente
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calcola il gradiente
    g = grad(x);
    gnorm = norm(g);
    % passo successivo
    xnew = x - alpha*g;
    % controllo sul passo successivo
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    % plot del passo corrente
    plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh
    % nuove condizioni di terminazione
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    
end
xopt = x;
fopt = f2(xopt);
niter = niter - 1;

% funzione gradiente
function g = grad(x)
g = [2*x(1) + x(2)
    x(1) + 6*x(2)];








