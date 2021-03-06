close all
clear all

% Add directory to current path
addpath('..')

% Parameters
n     = 100-1;                     % Number of interfaces
sigma = ones(n+1,1);               % Diffusivities 
for j=1:n+1;
    sigma(j)=sqrt(1.1+sin(j));
end
xj = linspace(0,1,n+2);
xj = xj(2:n+2);                   % Location of interfaces
u0    = @(x) x;                   % Initial condition
beta  = [0, 1, 1, 0];             % Boundary conditions
f1  = @(t) 0.;                    % RHS Boundary condition
f2  = @(t) 0.;                    % LHS Boundary condition
tspan = [0.1,0.3,5.];             % Times at which to compute solution
options.NX    = 15;               % Number of places to evaluate solution
options.NN    = 10;               % Integration bounds
options.Ny    = 200;              % Number of points to use in integration
H     = .5*ones(1,n);             % Contact coefficients
tic
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1,f2,tspan,'Imperfect',H, options);
toc

% Plot
figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9])
    hold on
end
for j=1:length(tspan)
    plot(xf,u(:,j),'r-','LineWidth',2.0)
end
axis([0,1,0,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExF.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xf,u(:,j),'r-','LineWidth',2.0)
    axis([0,1,0,1.1])
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExF_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end