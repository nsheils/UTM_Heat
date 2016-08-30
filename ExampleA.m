close all
clear all
clc

% Add directory to current path
addpath('..')

% Parameters
n     = 2;                        % Number of layers
sigma = [1,1,1];                  % Diffusivities 
xj    = [1/3,2/3,1];              % Location of interfaces
u0    = @(x) zeros(size(x));      % Initial condition
beta  = [1,0,1,0];                % Boundary conditions
f1  = @(t) .5;                    % RHS Boundary condition 1
f2  = @(t) 1;                     % RHS Boundary condition 2
%tspan = [0.001,0.01,0.1,0.2,1.0];% Times at which to compute solution
tspan = [1,2, 4, 8];
options.NX    = 15;               % Number of places to evaluate solution
options.NN    = 20;               % Integration bounds
options.Ny    = 200;              % Number of points to use in integration
tic
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1,f2,tspan,'Perfect',options);
toc


%EXACT solution.  Only works for Dirichlet data with 0 initial condition
clear f1; clear f2;
f1=.5; f2= 1.;
nspace=1:1000;
ue=cell(length(tspan),1);
for tau=1:length(tspan)
    ue{tau}= @(x) f1*(1-x)+f2*x+sum(2./(nspace.*pi).*(f2*(-1).^(nspace+1)-f1).*exp(-nspace.^2.*pi^2.*tspan(tau)).*sin(nspace.*pi.*x));
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
for j=1:length(tspan)
    plot(xf,u(:,j),'m--','LineWidth',2.0)
    plot(xf,uexact(:,j),'b-','LineWidth',2.0)
end
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
saveas(gcf,'ExA.pdf')

for j=1:length(tspan)
    figure
    for i = 1:n+1,
        plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
        hold on
    end
    plot(xf,u(:,j),'m--','LineWidth',2.0)
    plot(xf,uexact(:,j),'b-','LineWidth',2.0)
    xlabel('$x$','Interpreter','LaTeX','FontSize',20)
    ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
    title(['$t=$ ',num2str(tspan(j))],'Interpreter','LateX','FontSize',20)
    set(gca,'FontSize',14,'Layer','top')
    filename=['ExA_t', num2str(tspan(j))];
    filename(filename==['.'])=[];
    saveas(gcf,[filename '.pdf'])
end