## UTM_Heat: A MatLab code for solving the one-dimensional multilayer heat equation

``UTM_Heat`` Solves the one-dimensional heat equation using the Unified Transform Method (UTM) for a medium with `n` interfaces.  The code is applicable to both perfect and imperfect contact at the interaces between adjacent layers with any boundary (including time-dependent) conditions as well as general initial conditions.


## References
UTM_Heat is an implementation of the Unified Transform Method (UTM) as
laid out by N. E. Sheils in the paper "Imperfect Heat."

This work was inspired by the work of E.J. Carr and I.W. Turner and compares solutions to their ``MultDiff`` code as laid out in:

E. J. Carr and I. W. Turner, A semi-analytical solution for multilayer diffusion in a 
composite medium consisting of a large number of layers, Applied Mathematical Modelling (2016), 
http://dx.doi.org/10.1016/j.apm.2016.02.041


## Examples

### Example A
```
n     = 2;                         % Number of layers
sigma = [1,1,1];                   % Diffusivities 
xj    = [1/3,2/3,1];               % Location of interfaces
u0    = @(x) zeros(size(x));       % Initial condition
beta  = [1,0,.5,1,0,1.];           % Boundary conditions
%tspan = [0.001,0.01,0.1,0.2,1.0]; % Times at which to compute solution
tspan = [1,2, 4, 8];
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Perfect');
```

In this case we can find the exact solution as well

```
f1=beta(3)/beta(1); f2= beta(6)/beta(4);
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
```

We plot the solution using
```
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
```

<figure><img src="UTM_Heat/ExA.pdf"></figure>

### Example B
```
n     = 4-1;                      % Number of interfaces
sigma = sqrt([.2, .01, .1, 1.]);  % Diffusivities 
xj = linspace(0,1,n+2);
xj = xj(2:n+2);                   % Location of interfaces
u0    = @(x) ones(size(x));       % Initial condition
beta  = [1, 0, .2, .4, 1, .1];    % Boundary conditions
tspan = [.02,0.1,0.5,1,2,10];     % Times at which to compute solution
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Perfect');
```

We plot the solution using
```
figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9])
    hold on
end
for j=1:length(tspan)
    plot(xf,u(:,j),'m-','LineWidth',2.0)
end
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
```

<figure><img src="UTM_Heat/ExB.pdf"></figure>


### Example C
```
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
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Perfect');

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
```

<figure><img src="UTM_Heat/ExC.pdf"></figure>

### Example D
```
n     = 10-1;                      % Number of interfaces
sigma = ones(n+1,1);               % Diffusivities 
j=1;
while 2*j<=n+1
    sigma(2*j)=sqrt(.1);
    j=j+1;
end
xj    = linspace(0,1,n+2);
xj    = xj(2:n+2);                % Location of interfaces
u0    = @(x) zeros(size(x));       % Initial condition
beta  = [1, 0, 1, 0,1, 0];        % Boundary conditions
tspan = [0.007, 2., 10.];         % Times at which to compute solution
H     = .5*ones(1,n);             % Contact coefficients
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Imperfect', H);

figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9])
    hold on
end
for j=1:length(tspan)
    plot(xf,u(:,j),'m-','LineWidth',2.0)
end
axis([0,1,-.1,4])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
```

<figure><img src="UTM_Heat/ExD.pdf"></figure>


## License

See `LICENSE.md` for licensing information.