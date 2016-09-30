function [u,xf]=UTM_Heat(n,sigma,xj,u0,beta,f1,f2,tspan,interface,varargin)
% UTM_Heat Solves the one-dimensional multilayer diffusion problem using
%                 the Unified Transform Method (UTM).
%
%   UTM_Heat solves the heat equation in a one-dimensional composite slab
%   of finite length with multiple layers. The code is applicable to both
%   perfect and imperfect contact at the interfaces between adjacent layers
%   and any boundary conditions at the ends of the slab.
%
%   UTM_Heat is an implementation of the Unified Transform Method (UTM) as
%   laid out by N. E. Sheils in the paper "Imperfect Heat"
%
%
%   Description:
%   -----------------------------------------------------------------------
%   UTM_Heat solves the heat equation in each layer (x_{i-1} < x < x_{i}):
%
%      du_(i)/dt = d/dx * (kappa(i) * du_(i)/dx),   i = 1,...,n-1,
%
%   subject to the following initial and boundary conditions:
%
%      u_(i)(x,t) = u0(x)                                    at t = 0
%      beta1 * u_(1)(x,t) + beta2 * du_(1)/dx(x,t) = f1(t)   at x = 0
%      beta3 * u_(n)(x,t) + beta4 * du_(n+1)/dx(x,t) = f2(t)   at x = x_{n+1}
%
%   where u_(i) is the solution in layer i, sigma(i)=sqrt(kappa(i)) is the
%   square root ofdiffusivity in layer i (constant) and beta1, beta2,
%   beta3, and beta4 are constants.
%
%   Either perfect or imperfect contact is imposed at the interfaces.
%
%    - Perfect contact
%       u_(i)(x_i,t) = u_(i+1)(x_i,t)
%       kappa(i) * u_(i)(x_i,t) = kappa(i+1) * u_(i+1)(x_i,t)
%
%    - Imperfect contact
%       kappa(i)*du_(i)/dx(x_i,t) = H(i)*(u_(i+1)(x_i,t)-u_(i)(x_i,t))
%       kappa(i+1)*du_(i+1)/dx(x_i,t) = H(i)*(u_(i+1)(x_i,t)-u_(i)(x_i,t))
%
%   Usage:
%   -----------------------------------------------------------------------
%   [U,XF] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Perfect')
%   [U,XF] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Perfect',options)
%   [U,XF] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Imperfect',H)
%   [U,XF] = UTM_Heat(n,sigma,xj,u0,beta,tspan,'Imperfect',H,options)
%
%   Input Arguments:
%   -----------------------------------------------------------------------
%   n           Number of interfaces (n+1 layers). Must be an integer
%               greater than or equal to 2.
%   sigma       A vector of length n+1 containing the square root of the
%               diffusivity values in each layer such that the diffusivity
%               in Layer i is given by sigma(i).^2 (i = 1,...,n+1).
%   xj          A vector of length n+1 of the x coordinates of the
%               locations of the interfaces as well as the right boundary
%               of the slab.  The left boundary is assumed to be at 0.
%   u0          A function handle specifying the initial condition. The
%               function uint = u0(X) should accept a vector argument x and
%               return a vector result uint. Use array operators .*, ./ and
%               .^ in the definition of u0 so that it can be evaluated with
%               a vector argument.
%   beta        A vector of four values specifying the boudary conditons
%               beta=(beta1,beta2,beta3,beta4)
%   f1          A function handle specifying the RHS of the boundary
%               condition.
%   f2          A function handle specifying the RHS of the boundary
%               condition.
%   tspan       A vector specifying the times at which a solution is
%               requested. To obtain solutions at specific times
%               t0,t1,...,tf, use TSPAN = [t0,t1,...,tf].
%   interface   Internal boundary conditions at interfaces between adjacent
%               layers. inteface can be either 'Perfect' or 'Imperfect'.
%   H           A vector of length n containing the contact
%               transfer coeffecients at the interfaces between adjacent
%               layers such that the coefficient between layer i and layer
%               i+1 is given by H(i) (i = 1,..,n).
%               * Applicable to imperfect contant only.
%   options     An (optional) set of solver options. Fields in the
%               structure options are
%                - NX   number of divisions within each slab. U(:,j) gives
%                       the solution at xf =
%                       xj(i-1):(xj(i)-xj(i-1))/NX:xj(i) and t = tspan(j).
%                       NX does not change the accuracy of the code, just
%                       the places where the solution is evaluated.
%                       [NX = 15 by default]
%                - NN   Integration bounds (-NN, NN)
%                       [NN = 10 by default]
%                - Ny   number of points to use in integration.
%                       [Ny = 200 by default]
%
%   Output Arugments:
%   -----------------------------------------------------------------------
%   u   Matrix of solution values. u(:,j) gives the solution on the entire
%       slab (0 <= x < x_{n+1}) at time t = tspan(j) and at the grid points
%       returned in the output vector xf.
%   xf  Vector of grid points at which solution is given. xf is a vector
%       of length (n+1)*NX.

%   Example:
%   -----------------------------------------------------------------------
%   u0 = @(x) zeros(size(x));
%   [u,x] = UTM_Heat(2,[1,0.1,1],[0.3,0.7, 1.0],u0,[1,0,0,1,1,.5],
%           [0.02,0.05,0.2,0.5], @(t) .1, @(t) 1.,'Perfect');
%
% -------------------------------------------------------------------------
% Check inputs
% -------------------------------------------------------------------------
if nargin < 9
    error('Not enough input arguments.');
elseif nargin == 9
    if strcmp(interface,'Imperfect')
        error('H must be specified for imperfect contact at interfaces.');
    end
    options = struct;
elseif nargin == 10
    if strcmp(interface,'Perfect')
        options = varargin{1};
    elseif strcmp(interface,'Imperfect')
        H = varargin{1};
        options = struct;
    end
elseif nargin == 11
    if strcmp(interface,'Perfect')
        error('Too many input arguments for interface = ''Perfect''.');
    elseif strcmp(interface,'Imperfect')
        H = varargin{1};
        if length(H)~=n
            error('H must have n values');
        end
        options = varargin{2};
    end
else
    error('Too many input arguments.');
end

% Number of layers
if round(n) ~= n || n < 1
    error('n must be an integer greater than or equal to 1.')
end

% Diffusivities
if length(sigma) ~= n+1 || sum(sigma > 0) ~= n+1
    error('sigma must be a vector of length n+1 with sigma(i)>0.')
end

% Interfaces
if length(xj) ~= n+1 || all(diff(xj)>0)==0
    error('xj must be a vector of length n+1 with with increasing values.')
end

% Initial condition
if ~isa(u0,'function_handle') || nargin(u0) ~= 1
    error('u0 must be a function handle of the form uint = u0(x).');
end

% Boundary conditions
if length(beta) ~= 4
    error('beta must be a vector of length 4')
end
if ~isa(f1,'function_handle') || nargin(f1) ~= 1
    error('f1 must be a function handle of the form f1 = f1(t).');
end
if ~isa(f2,'function_handle') || nargin(f2) ~= 1
    error('f2 must be a function handle of the form f2 = f2(t).');
end

% Time vector
tlength = length(tspan);
if sum(tspan > 0) ~= tlength
    error('tspan must have entries that are greater than or equal to 0.')
end

% Internal boundary conditions at interfaces
if strcmp(interface,'Perfect') || strcmp(interface,'Imperfect')
else
    error('interface must be either ''Perfect'' or ''Imperfect''.')
end
% Check options structure
if ~isa(options,'struct')
    error('options must be a structure.')
end
Names = {'NX', 'NN', 'Ny'};
fn = fieldnames(options);
for i = 1:length(fn)
    j = strcmp(fn(i),Names);
    if sum(j) == 0
        error('Invalid option ''%s''.',fn{i});
    end
end
% Number of divisions within each slab
if isfield(options,'NX')
    NX = options.NX;
    if round(NX) ~= NX && NX < 1
        error('options.NX must be an integer greater than or equal to 1.')
    end
else
    NX = 15; % Default
end
% Range to use in integration
if isfield(options,'NN')
    NN = options.NN;
    if round(NN) ~= NN && NN < 5
        error('options.NN must be an integer greater than or equal to 5.')
    end
else
    NN = 10; % Default
end
% Number of points to use in integration
if isfield(options,'Ny')
    Ny = options.Ny;
    if round(Ny) ~= Ny && Ny < 5
        error('options.Ny must be an integer greater than or equal to 5.')
    end
else
    Ny = 200; % Default
end

% Check boundary conditions are implemented correctly
if beta(1) == 0 && beta(2) == 0
    error('Boundary condition is incorrect at left boundary.')
end
if beta(3) == 0 && beta(4) == 0
    error('Boundary condition is incorrect at right boundary.')
end

% -------------------------------------------------------------------------
% Grid spacing within each slab
% -------------------------------------------------------------------------
xgrid = zeros(NX+1,n+1);
% Slab 1 (First slab)
xgrid(:,1) = 0:xj(1)/NX:xj(1);
% Slabs 2,...,n+1
for i = 2:n+1
    xgrid(:,i) = xj(i-1):(xj(i)-xj(i-1))/NX:xj(i);
end
xf = reshape(xgrid,(NX+1)*(n+1),1);

% -------------------------------------------------------------------------
% Preliminaries
% -------------------------------------------------------------------------
A11=cell(n+2,n+2);
A11(:)={@(nu) 0};
A12=A11;
A21=A11;
A22=A11;

%Build the u0hat(k) functions
%CALLING THESE FUNCTIONS TAKES ALL THE TIME
u0hat=cell(n+1,1);
u0hat{1}= @(k) integral(@(x) exp(-1i.*k.*x).*u0(x),0,xj(1));
for j=2:n+1
   u0hat{j}= @(k) integral(@(x) exp(-1i.*k.*x).*u0(x),xj(j-1),xj(j));
end

Y=cell(2*n+4,length(tspan));
Y(:)={0};
thetaspace=linspace(-NN,NN,Ny+1);
%parameterize D+ and D-
kp= @(theta) 1i.*sin(pi/8-1i.*theta);
km= @(theta) -1i.*sin(pi/8-1i.*theta);

myYp=cell(2*n+4,length(tspan));
myYm=myYp;
myAp=cell(2*n+4,2*n+4);
myAm=myAp;

Apnu=cell(length(thetaspace));
Apnu(:)={zeros(2*n+4,2*n+4)};
Amnu=Apnu;
Ypnu=cell(length(thetaspace),1);
Ypnu(:)={zeros(2*n+4,length(tspan))};
Ymnu=Ypnu;

Xp=cell(length(thetaspace),length(tspan));
Xm=Xp;
% -------------------------------------------------------------------------
% IMPERFECT INTERFACE conditions
% -------------------------------------------------------------------------
if strcmp(interface,'Imperfect')
    
    % ---------------------------------------------------------------------
    % Build matrix A to solve for unknown functions
    % ---------------------------------------------------------------------
    %boundary conditions
    A11{1,1}=@(nu) beta(2);
    A11{1,2}=@(nu) beta(1);
    A22{n+2,n+1}=@(nu) beta(3);
    A22{n+2,n+2}=@(nu) beta(4);
    
    A11{2,1}= @(nu) -sigma(1).^2;
    A11{2,2}= @(nu) -1i.*sigma(1).*nu;
    A21{1,1}= @(nu) -sigma(1).^2;
    A21{1,2}= @(nu) 1i.*sigma(1).*nu;
    for j=1:n
        A11{j+1,j+2}=@(nu) arrayfun(@(j) H(j).*exp(-1i.*nu.*xj(j)./sigma(j)),j);
        A11{j+2,j+2}=@(nu) arrayfun(@(j) -(H(j)+1i.*sigma(j+1).*nu).*exp(-1i.*nu.*xj(j)./sigma(j+1)),j);
        A12{j+1,j}=@(nu) arrayfun(@(j) (1i.*sigma(j).*nu-H(j)).*exp(-1i.*nu.*xj(j)./sigma(j)),j);
        A12{j+2,j}=@(nu) arrayfun(@(j) H(j)*exp(-1i.*nu.*xj(j)./sigma(j+1)),j);
        A21{j,j+2}=@(nu) arrayfun(@(j) H(j)*exp(1i.*nu.*xj(j)/sigma(j)),j);
        A21{j+1,j+2}=@(nu) arrayfun(@(j) (1i.*sigma(j+1).*nu-H(j)).*exp(1i.*nu.*xj(j)./sigma(j+1)),j);
        A22{j,j}=@(nu) arrayfun(@(j) -(1i.*sigma(j).*nu+H(j)).*exp(1i.*nu.*xj(j)./sigma(j)),j);
        A22{j+1,j}=@(nu) arrayfun(@(j) H(j)*exp(1i.*nu.*xj(j)./sigma(j+1)),j);
    end
    
    A12{n+2,n+1}=@(nu) 1i.*sigma(n+1).*nu.*exp(-1i.*nu.*xj(n+1)./sigma(n+1));
    A12{n+2,n+2}=@(nu) sigma(n+1).^2.*exp(-1i.*nu.*xj(n+1)./sigma(n+1));
    A22{n+1,n+1}=@(nu) -1i.*sigma(n+1).*nu.*exp(1i.*nu.*xj(n+1)./sigma(n+1));
    A22{n+1,n+2}=@(nu) sigma(n+1).^2*exp(1i.*nu.*xj(n+1)./sigma(n+1));
    
    A=cell(2*n+4,2*n+4);
    for j=1:n+2
        for k=1:n+2
            A{j,k}=A11{j,k};
            A{j+n+2,k}=A21{j,k};
            A{j,k+n+2}=A12{j,k};
            A{j+n+2,k+n+2}=A22{j,k};
        end
    end
    
    % Note this Y is scaled
    for tau=1:length(tspan)
        % FOURIER TRANSFORM OF condition on left
        if f1(inf)==f1(0) && f1(10)==f1(.2) %If f1(s) is constant we can do the integration ourselves
            Y{1,tau}= @(nu) f1(tspan(tau)).*(1.-exp(-tspan(tau).*(nu).^2))./(nu).^2;
        else
            Y{1,tau}= @(nu) integral(@(s) exp(nu.^2.*(s-tspan(tau))).*f1(s),0,tspan(tau));
        end
        % FOURIER TRANSFORM OF condition on right
        if f2(inf)==f2(0) && f2(10)==f2(.2) %If f2(s) is constant we can do the integration ourselves
            Y{2*n+4,tau}= @(nu) f2(tspan(tau)).*(1.-exp(-tspan(tau).*(nu).^2))./(nu).^2;
        else
            Y{2*n+4,tau}= @(nu) integral(@(s) exp(nu.^2.*(s-tspan(tau))).*f2(s),0,tspan(tau));
        end
        for j=1:n+1
            Y{j+1,tau}=@(nu) arrayfun(@(j) -exp(-nu.^2*tspan(tau)).*u0hat{j}(nu./sigma(j)),j);
            Y{n+2+j,tau}=@(nu) arrayfun(@(j) -exp(-nu.^2*tspan(tau)).*u0hat{j}(-nu./sigma(j)),j);
        end
    end
    
     %Evaluate everything in the transformed theta interval.
    for j=1:2*n+4
        %Each column of myY corresponds to a given time
        for tau=1:length(tspan)
            myYp{j,tau}=arrayfun(Y{j,tau},kp(thetaspace)); 
            myYm{j,tau}=arrayfun(Y{j,tau},km(thetaspace));
            for th=1:length(thetaspace)
                Ypnu{th}(j,tau)=myYp{j,tau}(th);
                Ymnu{th}(j,tau)=myYm{j,tau}(th);
            end
        end
        for k=1:2*n+4
            myAp{j,k}=arrayfun(A{j,k},kp(thetaspace));
            myAm{j,k}=arrayfun(A{j,k},km(thetaspace));
            for th=1:length(thetaspace)
                Apnu{th}(j,k)=myAp{j,k}(th);
                Amnu{th}(j,k)=myAm{j,k}(th);
            end
        end
    end
    
    badp=zeros(length(thetaspace),1);
    badm=zeros(length(thetaspace),1);
    for tau=1:length(tspan)
        for j=1:length(thetaspace)
            if sum(sum(isinf(Apnu{j})))==0
                Xp{j,tau}=sparse(Apnu{j})\Ypnu{j}(:,tau);
            else
                Xp{j,tau}=zeros(2*n+4,1);
                badp(j)=1;
            end
            if sum(sum(isinf(Amnu{j})))==0
                Xm{j,tau}=sparse(Amnu{j})\Ymnu{j}(:,tau);
            else
                Xm{j,tau}=zeros(2*n+4,1);
                badm(j)=1;
            end
        end
    end
    
    g0p=cell(n+1,length(tspan));
    g0p(:)={zeros(length(thetaspace),1)};
    g0m=g0p;
    h0p=g0p;
    h0m=g0p;
    g1p=cell(length(tspan),1);
    g1p(:)={zeros(length(thetaspace),1)};
    h1m=g1p;
    for k=1:length(thetaspace)
        for tau=1:length(tspan)
            g1p{tau}(k)=Xp{k,tau}(1);
            for j=2:n+2
                g0p{j-1,tau}(k)=Xp{k,tau}(j);
                g0m{j-1,tau}(k)=Xm{k,tau}(j);
                h0p{j-1,tau}(k)=Xp{k,tau}(j+n+1);
                h0m{j-1,tau}(k)=Xm{k,tau}(j+n+1);
            end
            h1m{tau}(k)=Xm{k,tau}(2*n+4);
        end
    end
    
    %%Add a fit to fix when the X values are automatically set to 0.
    h0m_fitR=cell(length(tspan),1);
    h0m_fitI=h0m_fitR;
    h1m_fitR=h0m_fitR;
    h1m_fitI=h0m_fitR;
    g0p_fitR=h0m_fitR;
    g0p_fitI=h0m_fitR;
    g1p_fitR=h0m_fitR;
    g1p_fitI=h0m_fitR;
    g0m_fitR=h0m_fitR;
    g0m_fitI=h0m_fitR;
    min_m=find(badm==0,1);
    max_m=find(badm==0,1,'last');
    min_p=find(badp==0,1);
    max_p=find(badp==0,1,'last');
    for tau=1:length(tspan)
        h0m_fitR{tau}=fit(transpose(min_m:max_m),real(h0m{n+1,tau}(min_m:max_m)),'pchipinterp');
        h0m_fitI{tau}=fit(transpose(min_m:max_m),imag(h0m{n+1,tau}(min_m:max_m)),'pchipinterp');
        h1m_fitR{tau}=fit(transpose(min_m:max_m),real(h1m{tau}(min_m:max_m)),'pchipinterp');
        h1m_fitI{tau}=fit(transpose(min_m:max_m),imag(h1m{tau}(min_m:max_m)),'pchipinterp');
        g0p_fitR{tau}=fit(transpose(min_p:max_p),real(g0p{1,tau}(min_p:max_p)),'pchipinterp');
        g0p_fitI{tau}=fit(transpose(min_p:max_p),imag(g0p{1,tau}(min_p:max_p)),'pchipinterp');
        g1p_fitR{tau}=fit(transpose(min_p:max_p),real(g1p{tau}(min_p:max_p)),'pchipinterp');
        g1p_fitI{tau}=fit(transpose(min_p:max_p),imag(g1p{tau}(min_p:max_p)),'pchipinterp');
        g0m_fitR{tau}=fit(transpose(min_m:max_m),real(g0m{1,tau}(min_m:max_m)),'pchipinterp');
        g0m_fitI{tau}=fit(transpose(min_m:max_m),imag(g0m{1,tau}(min_m:max_m)),'pchipinterp');
    end
    for tau=1:length(tspan)
        for k=[1:min_p,max_p:length(thetaspace)]
            g0p{1,tau}(k)=g0p_fitR{tau}(k)+1i*g0p_fitI{tau}(k);
            g1p{tau}(k)=g1p_fitR{tau}(k)+1i*g1p_fitI{tau}(k);
        end
        for k=[1:min_m,max_m:length(thetaspace)]
            g0m{1,tau}(k)=g0m_fitR{tau}(k)+1i*g0m_fitI{tau}(k);
            h0m{n+1,tau}(k)=h0m_fitR{tau}(k)+1i*h0m_fitI{tau}(k);
            h1m{tau}(k)=h1m_fitR{tau}(k)+1i*h1m_fitI{tau}(k);
        end
    end
    
    % ---------------------------------------------------------------------
    % Solve for u
    % ---------------------------------------------------------------------
    usoln=cell(n+1,length(tspan));
    for s=1:length(tspan)
        usoln{1,s}=@(x) 1/(2*pi)*integral(@(k) exp(1i.*k.*x-(sigma(1).*k).^2.*tspan(s)).*arrayfun(u0hat{1},k),-Inf,Inf)...
            -1/(2*pi*sigma(1)).*trapz(thetaspace,(exp(1i.*km(thetaspace).*(x-xj(1))/sigma(1)).*(H(1).*transpose(g0m{2,s}(:))+(1i.*sigma(1).*km(thetaspace)-H(1)).*transpose(h0m{1,s}(:)))).*(-cos(pi/8-1i.*thetaspace)))...
            -1/(2*pi).*trapz(thetaspace,(exp(1i.*kp(thetaspace).*x/sigma(1)).*(sigma(1)*transpose(g1p{s}(:))+1i.*kp(thetaspace).*transpose(g0p{1,s}(:)))).*(cos(pi/8-1i.*thetaspace)));
        for j=2:n
            usoln{j,s}= @(x) 1/(2*pi)*integral(@(k) exp(1i.*k.*x-(sigma(j).*k).^2.*tspan(s)).*arrayfun(u0hat{j},k),-Inf,Inf)...
                -1/(2*pi*sigma(j)).*trapz(thetaspace,(exp(1i.*km(thetaspace).*(x-xj(j))/sigma(j)).*(H(j).*transpose(g0m{j+1,s}(:))+(1i.*sigma(j).*km(thetaspace)-H(j)).*transpose(h0m{j,s}(:)))).*(-cos(pi/8-1i.*thetaspace)))...
                -1/(2*pi*sigma(j)).*trapz(thetaspace,(exp(1i.*kp(thetaspace).*(x-xj(j-1))/sigma(j)).*((H(j-1)+1i.*sigma(j).*kp(thetaspace)).*transpose(g0p{j,s}(:))-H(j-1).*transpose(h0p{j-1,s}(:)))).*(cos(pi/8-1i.*thetaspace)));
        end
        usoln{n+1,s}= @(x) 1/(2*pi)*integral(@(k) exp(1i.*k.*x-(sigma(n+1).*k).^2.*tspan(s)).*arrayfun(u0hat{n+1},k),-Inf,Inf)...
            -1/(2*pi).*trapz(thetaspace,(exp(1i.*km(thetaspace).*(x-xj(n+1))/sigma(n+1)).*(sigma(n+1).*transpose(h1m{s}(:))+1i.*km(thetaspace).*transpose(h0m{n+1,s}(:)))).*(-cos(pi/8-1i.*thetaspace)))...
            -1/(2*pi*sigma(n+1)).*trapz(thetaspace,(exp(1i.*kp(thetaspace).*(x-xj(n))/sigma(n+1)).*((H(n)+1i.*sigma(n+1).*kp(thetaspace)).*transpose(g0p{n+1,s}(:))-H(n).*transpose(h0p{n,s}(:)))).*(cos(pi/8-1i.*thetaspace)));
    end
    
    u=zeros((n+1)*(NX+1),length(tspan));
    for s=1:length(tspan)
        for m=1:NX+1
            for j=1:n+1
                u((j-1)*(NX+1)+m,s)=real(usoln{j,s}(xgrid(m,j)));
            end
        end
    end
% -------------------------------------------------------------------------
% PERFECT INTERFACE conditions
% -------------------------------------------------------------------------
elseif strcmp(interface,'Perfect')
    
    % ---------------------------------------------------------------------
    % Build matrix A to solve for unknown functions
    % ---------------------------------------------------------------------
    %boundary conditions
    A11{1,1}=@(nu) beta(1);
    A12{1,1}=@(nu) beta(2);
    A21{n+2,n+2}=@(nu) beta(3);
    A22{n+2,n+2}=@(nu) beta(4);
    A11{2,1}=@(nu) -1i.* sigma(1).*nu;
    A12{2,1}=@(nu) -sigma(1)^2;
    A21{1,1}= @(nu) 1i.*sigma(1).*nu;
    A22{1,1}= @(nu) -sigma(1)^2;
    for j=1:n+1
        A11{j+1,j+1}= @(nu) arrayfun(@(j) 1i.*sigma(j).*nu.*exp(-1i.*nu.*xj(j)./sigma(j)),j);
        A21{j,j+1}= @(nu) arrayfun(@(j) -1i.*sigma(j).*nu.*exp(1i.*nu.*xj(j)./sigma(j)),j);
    end
    for j=1:n
        A11{j+2,j+1}= @(nu) arrayfun(@(j) -1i.*sigma(j+1).*nu.*exp(-1i.*nu.*xj(j)./sigma(j+1)),j);
        A12{j+1,j+1}= @(nu) arrayfun(@(j) sigma(j+1)^2.*exp(-1i.*nu.*xj(j)./sigma(j)),j);
        A12{j+2,j+1}= @(nu) arrayfun(@(j) -sigma(j+1)^2.*exp(-1i.*nu.*xj(j)./sigma(j+1)),j);
        A21{j+1,j+1}= @(nu) arrayfun(@(j) 1i.*sigma(j+1).*nu.*exp(1i.*nu.*xj(j)./sigma(j+1)),j);
        A22{j,j+1}= @(nu) arrayfun(@(j) sigma(j+1)^2.*exp(1i.*nu.*xj(j)./sigma(j)),j);
        A22{j+1,j+1}= @(nu) arrayfun(@(j) -sigma(j+1)^2.*exp(1i.*nu.*xj(j)./sigma(j+1)),j);
    end
    
    A12{n+2,n+2}= @(nu) sigma(n+1)^2.*exp(-1i.*nu.*xj(n+1)./sigma(n+1));
    A22{n+1,n+2}= @(nu) sigma(n+1)^2.*exp(1i.*nu.*xj(n+1)./sigma(n+1));
    
    A=cell(2*n+4,2*n+4);
    for j=1:n+2
        for k=1:n+2
            A{j,k}=A11{j,k};
            A{j+n+2,k}=A21{j,k};
            A{j,k+n+2}=A12{j,k};
            A{j+n+2,k+n+2}=A22{j,k};
        end
    end
    
    % Note this Y is scaled
    for tau=1:length(tspan)
        % FOURIER TRANSFORM OF condition on left
        if f1(inf)==f1(0) && f1(10)==f1(.2) %If f1(s) is constant we can do the integration ourselves
            Y{1,tau}= @(nu) f1(tspan(tau)).*(1.-exp(-tspan(tau).*(nu).^2))./(nu).^2;
        else
            Y{1,tau}= @(nu) integral(@(s) exp(nu.^2.*(s-tspan(tau))).*f1(s),0,tspan(tau));
        end
        % FOURIER TRANSFORM OF condition on right
        if f2(inf)==f2(0) && f2(10)==f2(.2) %If f2(s) is constant we can do the integration ourselves
            Y{2*n+4,tau}= @(nu) f2(tspan(tau)).*(1.-exp(-tspan(tau).*(nu).^2))./(nu).^2;
        else
            Y{2*n+4,tau}= @(nu) integral(@(s) exp(nu.^2.*(s-tspan(tau))).*f2(s),0,tspan(tau));
        end
        %Y doesn't play nicely with previously defined u0hats so we plug in
        % full formula here.
        Y{1+1,tau}=@(nu) -integral(@(x) exp(-nu.^2*tspan(tau)-1i.*nu./sigma(1).*x).*u0(x),0,xj(1));
        Y{n+2+1,tau}=@(nu) -integral(@(x) exp(-nu.^2*tspan(tau)+1i.*nu./sigma(1).*x).*u0(x),0,xj(1));
        for j=2:n+1
            Y{j+1,tau}=@(nu) arrayfun(@(j) -integral(@(x) exp(-nu.^2*tspan(tau)-1i.*nu./sigma(j).*x).*u0(x),xj(j-1),xj(j)),j);
            Y{n+2+j,tau}=@(nu) arrayfun(@(j) -integral(@(x) exp(-nu.^2*tspan(tau)+1i.*nu./sigma(j).*x).*u0(x),xj(j-1),xj(j)),j);
        end
    end
    
    %Evaluate everything in the transformed theta interval.
    for j=1:2*n+4
        %Each column of myY corresponds to a given time
        for tau=1:length(tspan)
            myYp{j,tau}=arrayfun(Y{j,tau},kp(thetaspace)); 
            myYm{j,tau}=arrayfun(Y{j,tau},km(thetaspace));
            for th=1:length(thetaspace)
                Ypnu{th}(j,tau)=myYp{j,tau}(th);
                Ymnu{th}(j,tau)=myYm{j,tau}(th);
            end
        end
        for k=1:2*n+4
            myAp{j,k}=arrayfun(A{j,k},kp(thetaspace));
            myAm{j,k}=arrayfun(A{j,k},km(thetaspace));
            for th=1:length(thetaspace)
                Apnu{th}(j,k)=myAp{j,k}(th);
                Amnu{th}(j,k)=myAm{j,k}(th);
            end
        end
    end
    
    badp=zeros(length(thetaspace),1);
    badm=zeros(length(thetaspace),1);
    for tau=1:length(tspan)
        for j=1:length(thetaspace)
            if sum(sum(isinf(Apnu{j})))==0
                Xp{j,tau}=sparse(Apnu{j})\Ypnu{j}(:,tau);
            else
                Xp{j,tau}=zeros(2*n+4,1);
                badp(j)=1;
            end
            if sum(sum(isinf(Amnu{j})))==0
                Xm{j,tau}=sparse(Amnu{j})\Ymnu{j}(:,tau);
            else
                Xm{j,tau}=zeros(2*n+4,1);
                badm(j)=1;
            end
        end
    end
    
    g0p=cell(n+1,length(tspan));
    g0p(:)={zeros(length(thetaspace),1)};
    g0m=g0p;
    g1p=g0p;
    g1m=g0p;
    h0n1m=cell(length(tspan),1);
    h0n1m(:)={zeros(length(thetaspace),1)};
    h1n1m=h0n1m;
    for tau=1:length(tspan)
        for k=1:length(thetaspace)
            for j=1:n+1
                g0p{j,tau}(k)=Xp{k,tau}(j);
                g0m{j,tau}(k)=Xm{k,tau}(j);
                g1p{j,tau}(k)=Xp{k,tau}(j+n+2);
                g1m{j,tau}(k)=Xm{k,tau}(j+n+2);
            end
            h0n1m{tau}(k)=Xm{k,tau}(n+2);
            h1n1m{tau}(k)=Xm{k,tau}(2*n+4);
        end
    end
    
    %%Add a fit to fix when the X values are automatically set to 0.
    h0n1m_fitR=cell(length(tspan),1);
    h0n1m_fitI=h0n1m_fitR;
    h1n1m_fitR=h0n1m_fitR;
    h1n1m_fitI=h0n1m_fitR;
    g0p_fitR=h0n1m_fitR;
    g0p_fitI=h0n1m_fitR;
    g1p_fitR=h0n1m_fitR;
    g1p_fitI=h0n1m_fitR;
    g0m_fitR=h0n1m_fitR;
    g0m_fitI=h0n1m_fitR;
    min_m=find(badm==0,1);
    max_m=find(badm==0,1,'last');
    min_p=find(badp==0,1);
    max_p=find(badp==0,1,'last');
    for tau=1:length(tspan)
        h0n1m_fitR{tau}=fit(transpose(min_m:max_m),real(h0n1m{tau}(min_m:max_m)),'pchipinterp');
        h0n1m_fitI{tau}=fit(transpose(min_m:max_m),imag(h0n1m{tau}(min_m:max_m)),'pchipinterp');
        h1n1m_fitR{tau}=fit(transpose(min_m:max_m),real(h1n1m{tau}(min_m:max_m)),'pchipinterp');
        h1n1m_fitI{tau}=fit(transpose(min_m:max_m),imag(h1n1m{tau}(min_m:max_m)),'pchipinterp');
        g0p_fitR{tau}=fit(transpose(min_p:max_p),real(g0p{1,tau}(min_p:max_p)),'pchipinterp');
        g0p_fitI{tau}=fit(transpose(min_p:max_p),imag(g0p{1,tau}(min_p:max_p)),'pchipinterp');
        g1p_fitR{tau}=fit(transpose(min_p:max_p),real(g1p{1,tau}(min_p:max_p)),'pchipinterp');
        g1p_fitI{tau}=fit(transpose(min_p:max_p),imag(g1p{1,tau}(min_p:max_p)),'pchipinterp');
        g0m_fitR{tau}=fit(transpose(min_m:max_m),real(g0m{1,tau}(min_m:max_m)),'pchipinterp');
        g0m_fitI{tau}=fit(transpose(min_m:max_m),imag(g0m{1,tau}(min_m:max_m)),'pchipinterp');
    end
    for tau=1:length(tspan)
        for k=[1:min_p,max_p:length(thetaspace)]
            g0p{1,tau}(k)=g0p_fitR{tau}(k)+1i*g0p_fitI{tau}(k);
            g1p{1,tau}(k)=g1p_fitR{tau}(k)+1i*g1p_fitI{tau}(k);
        end
        for k=[1:min_m,max_m:length(thetaspace)]
            g0m{1,tau}(k)=g0m_fitR{tau}(k)+1i*g0m_fitI{tau}(k);
            h0n1m{tau}(k)=h0n1m_fitR{tau}(k)+1i*h0n1m_fitI{tau}(k);
            h1n1m{tau}(k)=h1n1m_fitR{tau}(k)+1i*h1n1m_fitI{tau}(k);
        end
    end 
    
    % ---------------------------------------------------------------------
    % Solve for u
    % ---------------------------------------------------------------------
    usoln=cell(n+1,length(tspan));
    for s=1:length(tspan)
        usoln{1,s}= @(x) 1/(2*pi)*integral(@(k) exp(1i.*k.*x-(sigma(1).*k).^2*tspan(s)).*arrayfun(u0hat{1},k), -Inf, Inf)...
            -1/(2*pi).*trapz(thetaspace,exp(1i.*km(thetaspace).*(x-xj(1))/sigma(1)).*(sigma(2)^2/sigma(1).*transpose(g1m{2,s})+1i.*km(thetaspace).*transpose(g0m{2,s})).*(-cos(pi/8-1i.*thetaspace)))...
            -1/(2*pi).*trapz(thetaspace,exp(1i.*kp(thetaspace).*(x)/sigma(1)).*(sigma(1).*transpose(g1p{1,s})+1i.*kp(thetaspace).*transpose(g0p{1,s})).*(cos(pi/8-1i.*thetaspace)));
        if n>1
            for j=2:n
                usoln{j,s}= @(x) 1/(2*pi)*integral(@(k) exp(1i.*k.*x-(sigma(j).*k).^2.*tspan(s)).*arrayfun(u0hat{j},k),-Inf, Inf)...
                    -1/(2*pi).*trapz(thetaspace,(exp(1i.*km(thetaspace).*(x-xj(j))/sigma(j)).*(sigma(j+1)^2/sigma(j).*transpose(g1m{j+1,s}(:))+1i.*km(thetaspace).*transpose(g0m{j+1,s}(:)))).*(-cos(pi/8-1i.*thetaspace)))...
                    -1/(2*pi).*trapz(thetaspace,(exp(1i.*kp(thetaspace).*(x-xj(j-1))/sigma(j)).*(sigma(j).*transpose(g1p{j,s}(:))+1i.*kp(thetaspace).*transpose(g0p{j,s}(:)))).*(cos(pi/8-1i.*thetaspace)));
            end
        end
        usoln{n+1,s}= @(x) 1/(2*pi)*integral(@(k) exp(1i.*k.*x-(sigma(n+1).*k).^2.*tspan(s)).*arrayfun(u0hat{n+1},k),-Inf,Inf)...
            -1/(2*pi).*trapz(thetaspace,exp(1i.*km(thetaspace).*(x-xj(n+1))/sigma(n+1)).*(sigma(n+1).*transpose(h1n1m{s})+1i.*km(thetaspace).*transpose(h0n1m{s})).*(-cos(pi/8-1i.*thetaspace)))...
            -1/(2*pi).*trapz(thetaspace,exp(1i.*kp(thetaspace).*(x-xj(n))/sigma(n+1)).*(sigma(n+1).*transpose(g1p{n+1,s})+1i.*kp(thetaspace).*transpose(g0p{n+1,s})).*(cos(pi/8-1i.*thetaspace)));
    end
    
    u=zeros((n+1)*(NX+1),length(tspan));
    for s=1:length(tspan)
        for m=1:NX+1
            for j=1:n+1
                u((j-1)*(NX+1)+m,s)=real(usoln{j,s}(xgrid(m,j)));
            end
        end
    end
    
end