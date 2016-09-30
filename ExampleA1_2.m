close all
clear options

% Add directory to current path
addpath('..')
addpath('/Users/nataliesheils/Dropbox (uwamath)/InterfaceProblems/Heat_Imperfect/MultDiff-master')
addpath('/Users/nataliesheils/Dropbox (uwamath)/InterfaceProblems/Heat_Imperfect/MultDiff-master/analytical')
addpath('/Users/nataliesheils/Documents/MATLAB')


% Parameters
m     = 3;                          % Number of layers
kappa = ones(m,1);                   % Diffusivities 
l0    = 0.0;                         % Left end of slab
lm    = 1.0;                         % Right end of slab
dx    = (lm-l0)/m; l = dx:dx:lm-dx;  % Location of interfaces
u0    = @(x) x.^3;                   % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,0.0};   % Boundary condition (x = l0)
Rbnd  = {'Dirichlet',1.0,0.0,0.0};   % Boundary condition (x = lm)
tspan = [.001,.01,.1];               % Times at which to compute solution
options.NX    = 25;                  % Number of places to evaluate solution
options.N = 100;                     % Number of eigenvalues
tic
[umd,xmd] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect',options);
toc
[ug,xg] = multdiff_global(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect',options);

% Plot
figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
plot(xmd,umd,'b','LineWidth',2.0)
plot(xf,u,'r--','LineWidth',2.0)
plot(xg,ug,'go')
axis([0,1,0,max(max(u))])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExA_2_compare.pdf')