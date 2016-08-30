close all, clc

% Add directory to current path
addpath('..')

% Parameters
m     = 200;                          % Number of layers
kappa = ones(m,1);               % Diffusivities
j=1;
for j=1:m;
    kappa(j)=(1.1+sin(j));
end
l0    = 0.0;                         % Left end of slab
lm    = 1.0;                         % Right end of slab
dx    = (lm-l0)/m; l = dx:dx:lm-dx;  % Location of interfaces
u0    = @(x) ones(size(x));         % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,.5};   % Boundary condition (x = l0)
Rbnd  = {'Dirichlet',1.,0.0,0};   % Boundary condition (x = lm)
tspan = [.01,0.1,0.3,5.];      % Times at which to compute solution
H     = .5*ones(1,n);             % Contact coefficients
tic
[umd,xmd] = multdiff(m,kappa,l0,lm,l,u0,Lbnd,Rbnd,tspan,'Imperfect',H);
toc

% Plot
figure;
for i = 1:m-1,
    plot([l(i),l(i)],[-0.1,1.1],'Color',[0.9,0.9,0.9])
    hold on
end
plot(xmd,umd,'b','LineWidth',2.0)
plot(xf,u,'m--','LineWidth',2.0)
axis([0,1,0,1.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExF_compare.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xmd,umd(:,j),'b','LineWidth',2.0)
    plot(xf,u(:,j),'m--','LineWidth',2.0)
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    axis([0,1,-.1,1.1])
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExF_compare_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end