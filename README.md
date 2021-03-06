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
u0    = @(x) x.^3;                 % Initial condition
beta  = [1,0,1,0];                 % Boundary conditions
f1  = @(t) 0;                      % LHS Boundary condition
f2  = @(t) 1.;                     % RHS Boundary condition
tspan = [.01,.1,1];                % Times at which to compute solution
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1,f2,tspan,'Perfect';
```

In this case we can find the exact solution as well

```
clear f1; clear f2;
f1=0; f2= 1.;
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
```

We plot the solution using
```
figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u)),max(max(u))],'Color',[0.9,0.9,0.9], 'LineWidth',1.)
    hold on
end
plot(xf,uexact,'black-','LineWidth',2.0)
plot(xf,u,'r--','LineWidth',2.0)
axis([0,1,0,1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
```

<figure><img src="https://github.com/nsheils/UTM_Heat/blob/master/ExA.png"></figure>

### Example B
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
beta  = [1, 0, 1,0];              % Boundary conditions
f1  = @(t) 1;                     % LHS Boundary condition
f2  = @(t) 0;                     % RHS Boundary condition
tspan = [0.0005, .01, .2, 1.];    % Times at which to compute solution
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1,f2,tspan,'Perfect');

figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9])
    hold on
end
for j=1:length(tspan)
    plot(xf,u(:,j),'r-','LineWidth',2.0)
end
axis([0,1,0,2])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')
```

<figure><img src="https://github.com/nsheils/UTM_Heat/blob/master/ExB.png"></figure>

### Example C
```
n     = 4-1;                      % Number of interfaces
sigma = sqrt([.2, .01, .1, 1.]);  % Diffusivities 
xj = linspace(0,1,n+2);
xj = xj(2:n+2);                   % Location of interfaces
u0    = @(x) ones(size(x));       % Initial condition
beta  = [1, 0, 1, 1];             % Boundary conditions
f1  = @(t) cos(t);                % LHS Boundary condition
f2  = @(t) 0;                     % RHS Boundary condition
tspan = [0.5,1,2,10];             % Times at which to compute solution
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1, f2, tspan,'Perfect',options);

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

```

<figure><img src="https://github.com/nsheils/UTM_Heat/blob/master/ExC.png"></figure>

### Example D
```
n     = 10-1;                      % Number of interfaces
sigma = ones(n+1,1);               % Diffusivities 
j=1;
while 2*j<=n+1
    sigma(2*j)=sqrt(.1);
    j=j+1;
end
xj = linspace(0,1,n+2);
xj = xj(2:n+2);                   % Location of interfaces
u0    = @(x) zeros(size(x));      % Initial condition
beta  = [1, 0, 0, 1];             % Boundary conditions
f1  = @(t) 1;                	  % LHS Boundary condition
f2  = @(t) 0;                     % RHS Boundary condition
tspan = [0.007,2,10];             % Times at which to compute solution
H     = .5*ones(1,n);             % Contact coefficients
options.NX    = 15;               % Number of places to evaluate solution
options.NN    = 10;               % Integration bounds
options.Ny    = 200;              % Number of points to use in integration
[u,xf] = UTM_Heat(n,sigma,xj,u0,beta,f1, f2, tspan,’Imperfect’,H,options);

figure;
for i = 1:n+1,
    plot([xj(i),xj(i)],[min(min(u))-.1,max(max(u))+.1],'Color',[0.9,0.9,0.9])
hold on
end
for j=1:length(tspan)
    plot(xf,u(:,j),'r-','LineWidth',2.0)
end
axis([0,1,-1,max(max(u))+.1])
xlabel('$x$','Interpreter','LaTeX','FontSize',20)
ylabel('$u(x,t)$','Interpreter','LaTeX','FontSize',20)
set(gca,'FontSize',14,'Layer','top')

```

<figure><img src="https://github.com/nsheils/UTM_Heat/blob/master/ExD.png"></figure>


## License

See `LICENSE.md` for licensing information.