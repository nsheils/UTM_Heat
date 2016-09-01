close all
clc

% Add directory to current path
addpath('..')

% Parameters
n     = 200-1;                     % Number of interfaces
sigma = ones(n+1,1);               % Diffusivities 
for j=1:n+1;
    sigma(j)=sqrt(1.1+sin(j));
end
xj = linspace(0,1,n+2);
xj = xj(2:n+2);                   % Location of interfaces
u0    = @(x) ones(size(x));       % Initial condition
beta  = [1, 0, 1, 0];             % Boundary conditions
f1  = @(t) .5;                    % RHS Boundary condition 1
f2  = @(t) 0;                     % RHS Boundary condition 2
tspan = [.01,0.1,0.3,5.];         % Times at which to compute solution
options.NX    = 15;               % Number of places to evaluate solution
options.NN    = 20;               % Integration bounds
options.Ny    = 200;              % Number of points to use in integration
tic
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1, f2, tspan,'Perfect',options);
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
%axis([0,1,-0.1,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
axis([0,1,0,1.01]);
saveas(gcf,'ExE.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xf,u(:,j),'m-','LineWidth',2.0)
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    axis([0,1,0,1.01]);
    filename=['ExE_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end