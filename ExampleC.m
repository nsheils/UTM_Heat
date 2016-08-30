close all
clear all
clc

% Add directory to current path
addpath('..')

% Parameters
n     = 10-1;                      % Number of interfaces
sigma = ones(n+1,1);               % Diffusivities 
j=1;
while 2*j<=n+1
    sigma(2*j)=sqrt(.1);
    j=j+1;
end
xj    = linspace(0,1,n+2);
xj    = xj(2:n+2);                % Location of interfaces
u0    = @(x) zeros(size(x));      % Initial condition
beta  = [1, 0, 1, 1,0, 0];        % Boundary conditions
tspan = [0.0005, .01, .2, 1.];    % Times at which to compute solution
options.NX    = 15;               % Number of places to evaluate solution
options.NN    = 20;               % Integration bounds
options.Ny    = 200;              % Number of points to use in integration
tic
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Perfect', options);
toc

% Plot
figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9])
    hold on
end
for j=1:length(tspan)
    plot(xf,u(:,j),'m-','LineWidth',2.0)
end
axis([0,1,0,2])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExC.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xf,u(:,j),'m-','LineWidth',2.0)
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    axis([0,1,0,2])
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExC_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end