close all
clear options
%clc

% Add directory to current path
addpath('..')

% Parameters
n     = 2;                        % Number of interfaces
sigma = [1,1,1];                  % Diffusivities
xj    = [1/3,2/3,1];              % Location of interfaces
u0    = @(x) x.^3;                % Initial condition
beta  = [1,0,1,0];                % Boundary conditions
f1  = @(t) 0.;                    % LHS Boundary condition
f2  = @(t) 0.;                    % RHS Boundary condition
tspan = [.001,.01,.1];            % Times at which to compute solution
options.NX    = 25;               % Number of places to evaluate solution
options.NN    = 10;               % Integration bounds
options.Ny    = 200;              % Number of points to use in integration
tic
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1,f2,tspan,'Perfect',options);
toc

%EXACT solution.  Only works for Dirichlet data
clear f1; clear f2;
f1=0; f2= 0.;
nspace=1:1000;
ue=cell(length(tspan),1);
c=zeros(1,length(nspace));
intx=linspace(0,1,200);
for j=1:length(nspace)
    c(1,j)=2*trapz(intx,(u0(intx)-f1.*(1-intx)-f2.*intx).*sin(j*pi*intx));
end
for tau=1:length(tspan)
    ue{tau}= @(x) f1*(1-x)+f2*x+sum(exp(-nspace.^2.*pi^2.*tspan(tau)).*sin(nspace.*pi.*x).*c);
end
uexact=zeros(length(xf),length(tspan));
for j=1:length(tspan)
    for y=1:length(xf)
        uexact(y,j)=ue{j}(xf(y));
    end
end

% Plot
figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
    hold on
end
plot(xf,uexact,'black-','LineWidth',2.0)
plot(xf,u,'r--','LineWidth',2.0)
axis([0,1,0,max(max(u))])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExA_2.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xf,uexact(:,j),'black-','LineWidth',2.0)
    plot(xf,u(:,j),'r--','LineWidth',2.0)
    axis([0,1,0,max(max(u))])
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExA_2_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end

%%For Paper
figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
    hold on
end
plot(xf,uexact,'black-','LineWidth',2.0)
plot(xf,u,'r--','LineWidth',2.0)
axis([0,1,0,max(max(u))])
ax=gca;
ax.XTickLabel={};
ax.YTickLabel={};
saveas(gcf,'ExA_2_p.pdf')
