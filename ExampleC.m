close all
clear all

% Add directory to current path
addpath('..')

% Parameters
n     = 4-1;                      % Number of interfaces
sigma = sqrt([.2, .01, .1, 1.]);  % Diffusivities 
xj = linspace(0,1,n+2);
xj = xj(2:n+2);                   % Location of interfaces
u0    = @(x) ones(size(x));       % Initial condition
beta  = [1, 0, 1, 1];             % Boundary conditions
f1  = @(t) cos(t);                % LHS Boundary condition
f2  = @(t) 0.;                    % RHS Boundary condition
tspan = [0.5,1,2,10];             % Times at which to compute solution
options.NX    = 15;               % Number of places to evaluate solution
options.NN    = 10;               % Integration bounds
options.Ny    = 200;              % Number of points to use in integration
tic
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1, f2, tspan,'Perfect',options);
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
axis([0,1,-1,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExC.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xf,u(:,j),'r-','LineWidth',2.0)
    axis([0,1,-1,1.1])
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExC_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end

%%For Paper
figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9])
    hold on
end
for j=1:length(tspan)
    plot(xf,u(:,j),'r-','LineWidth',2.0)
end
axis([0,1,-1,1.1])
ax=gca;
ax.XTickLabel={};
ax.YTickLabel={};
saveas(gcf,'ExC_p.pdf')