close all

% Add directory to current path
addpath('..')

% Parameters
m     = 10;                         % Number of layers
kappa = ones(m,1);                  % Diffusivities 
j=1;
while 2*j<=m
    kappa(2*j)=.1;
    j=j+1;
end
l0    = 0.0;                        % Left end of slab
lm    = 1.0;                        % Right end of slab
l     = [0.1:.1:.9];                  % Location of interfaces
u0    = @(x) zeros(size(x));        % Initial condition
Lbnd  = {'Dirichlet',1.0,0.0,1.0};  % Boundary condition (x = l0)
Rbnd  = {'Neumann',0,1,0};      % Boundary condition (x = lm)
tspan = [0.007,2,10];    % Times at which to compute solution
H     = .5*ones(m-1,1);
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
plot(xf,u,'r--','LineWidth',2.0)
axis([0,1,-.1,2])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExD_compare.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xmd,umd(:,j),'b','LineWidth',2.0)
    plot(xf,u(:,j),'r--','LineWidth',2.0)
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    axis([0,1,-.1,2])
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExD_compare_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end

%%For Paper
figure;
for i = 1:m-1, 
    plot([l(i),l(i)],[-0.1,2],'Color',[0.9,0.9,0.9])
    hold on
end
plot(xmd,umd,'b','LineWidth',2.0)
plot(xf,u,'r--','LineWidth',2.0)
axis([0,1,-.1,2])
ax=gca;
ax.XTickLabel={};
ax.YTickLabel={};
saveas(gcf,'ExD_p.pdf')