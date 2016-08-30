close all %, clear all, clc

% Add directory to current path
addpath('..')

% Parameters
m     = 4;                           % Number of layers
kappa = [.2, .01, .1, 1];                % Diffusivities 
l0    = 0.0;                           % Left end of slab
lm    = 1.0;                           % Right end of slab
l    = [.25, .5, .75];    % Location of interfaces
u0    = @(x) ones(size(x));            % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,0.2};       % Boundary condition (x = l0)
Rbnd  = {'Robin',.4,1.0,0.1};     % Boundary condition (x = lm)
tspan = [0.02,0.1,0.5,1.0,2.0,10.0];  % Times at which to compute solution
tic
[umd,xmd] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Perfect');
toc

% Plot
figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
plot(xmd,umd,'b','LineWidth',2.0)
plot(xf,u,'m--','LineWidth',2.0)
axis([0,1,-0,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExB_compare.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xmd,umd(:,j),'b','LineWidth',2.0)
    plot(xf,u(:,j),'m--','LineWidth',2.0)
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    axis([0,1,-0,1.1])
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExB_compare_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end