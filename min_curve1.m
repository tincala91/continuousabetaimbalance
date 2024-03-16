function [b,status] = min_curve1(b0,CSF,PET)

% Estimates parameters of a hyperbolic regression model 
% (a modified version of the Michaelis-Menten equation)
% 
% FORMAT [b,status] = min_curve1(b0,CSF,PET)
% 
%         b0     - 3x1 array containing the initial parameters of the hyperbolic curve
%                  according to the format [a b c]
%                  
%         CSF    - 1xN array containing the observed inidividual CSF values
%         
%         PET   - 1xN arrat containing the observed individual PET uptake values 
%         
%         b      - 3x1 array containing the optimized parameters 
%         
%         status - 1x1 array containing values 0 or 1
%                  0 - routine failed
%                  1 - routine successfull
%
%--------------------------------------------------------------
% This function fits a hyperbolic regression model (a modified version of the Michaelis-Menten equation)
% to an observed data distribution. The model is proposed in the context of describing
% the relationship between CSF Ab42 concentration and Amyloid PET burden. The model was fitted by minimizing
% the sum of the Euclidean distance of the experimental points to the fitted line. Convergence is assumed when 
% the value of the objective function does not change above 10^-5 from one iterative step to the next one or 
% when the size of the step is smaller than 10^-5. This fitting strategy based on Euclidean distances was 
% preferred over the traditional fitting strategy based on minimizing the sum of the squared vertical distance, 
% to better describe the distribution of observations with lower PET uptake, whose shape is particularly 
% elongated over the y-axis. 
%
% Inputs: starting parameters to initialize the hyperbolic model; observed
% values for both CSF and PET measures.
%
% Outpus: optimized parameters of the hyperbolic model; status
% (failed/success)
% --------------------------------------------------------------
% Author: Juan Doming Gispert, BBRC-Fundació Pasqual Maragall
% Date: 15.03.2024
%
% The code is provided under GPLv3 license 
%
% Any publication based on this code should cite:
% Mastenbroek, Sala et al., A continuous amyloid-β CSF/PET imbalance model
% to capture Alzheimer’s disease heterogeneity (2024). Neurology
% 
% --------------------------------------------------------------
% 
% Ref:
% Mastenbroek, Sala et al., A continuous amyloid-β CSF/PET imbalance model
% to capture Alzheimer’s disease heterogeneity (2024). Neurology
%
%---------------------------------------------------------------

%A scaling factor is computed so as the balance the weight of each variable on the
%model estimation 
ScalingFactor = std(CSF)/std(PET); 

close all

%Set number of points of the fitted curve to be estimated
Npoints = 10000;
%PET values (typically PET) are scaled 
PET=PET*ScalingFactor;

%generate vector of 10000 points (equally spaced) from min to max PET
%value
x=linspace(min(PET),max(PET),Npoints); 

%declare function xyfit with parameters (inputs) x and b, that will call
%eval_func with inouts x b scaingfactor
xyfit = @(x,b) eval_func(x,b,ScalingFactor); 

%errors are estimated as Euclidean distance of each point to the fitted curve
err = @(b)dist2_points_curve([PET CSF],xyfit(x',b));

%displays relevant parameters in fitted model across interactions
opts = optimset('Display','iter');

%finds the solution that minimizes eucliean distances of curve to line

LB = [0 -10000 -20];
UB =[10000 0 0];

[b,~,status] = fminsearchbnd(err,b0,LB,UB,opts)

end

%parametrization of the curve

function xyfit = eval_func(x,b,ScalingFactor)
tol = 0.01;

xfit = x(find(x>b(3)+tol));

yfit = b(1)+((b(2)*xfit/ScalingFactor)./((xfit/ScalingFactor)-b(3)));

xyfit=[xfit yfit];
end


%given data points and curve, it computes distance between the two and
%stores (overall) them into an output d
function d=dist2_points_curve(points,curve);

for i = 1:length(points)
    for j = 1:length(curve)
        %computes the squared sum of all (absolute) distances between points and curve
        d(i,j)=sqrt(sum((points(i,:)-curve(j,:)).^2));
    end
end
d=sum(min(d,[],2).^2);
end

function [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBND: FMINSEARCH, but with bound constraints by transformation
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
% 
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example usage:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06
% size checks
xsize = size(x0);
x0 = x0(:);
n=length(x0);
if (nargin<3) || isempty(LB)
  LB = repmat(-inf,n,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = repmat(inf,n,1);
else
  UB = UB(:);
end
if (n~=length(LB)) || (n~=length(UB))
  error 'x0 is incompatible in size with either LB or UB.'
end
% set default options if necessary
if (nargin<5) || isempty(options)
  options = optimset('fminsearch');
end
% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.n = n;
% note that the number of parameters may actually vary if 
% a user has chosen to fix one or more parameters
params.xsize = xsize;
params.OutputFcn = [];
% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(n,1);
for i=1:n
  k = isfinite(LB(i)) + 2*isfinite(UB(i));
  params.BoundClass(i) = k;
  if (k==3) && (LB(i)==UB(i))
    params.BoundClass(i) = 4;
  end
end
% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      if x0(i)<=LB(i)
        % infeasible starting value. Use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(x0(i) - LB(i));
      end
      
      % increment k
      k=k+1;
    case 2
      % upper bound only
      if x0(i)>=UB(i)
        % infeasible starting value. use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(UB(i) - x0(i));
      end
      
      % increment k
      k=k+1;
    case 3
      % lower and upper bounds
      if x0(i)<=LB(i)
        % infeasible starting value
        x0u(k) = -pi/2;
      elseif x0(i)>=UB(i)
        % infeasible starting value
        x0u(k) = pi/2;
      else
        x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
      end
      
      % increment k
      k=k+1;
    case 0
      % unconstrained variable. x0u(i) is set.
      x0u(k) = x0(i);
      
      % increment k
      k=k+1;
    case 4
      % fixed variable. drop it before fminsearch sees it.
      % k is not incremented for this variable.
  end
  
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=n
  x0u(k:n) = [];
end
% were all the variables fixed?
if isempty(x0u)
  % All variables were fixed. quit immediately, setting the
  % appropriate parameters, then return.
  
  % undo the variable transformations into the original space
  x = xtransform(x0u,params);
  
  % final reshape
  x = reshape(x,xsize);
  
  % stuff fval with the final value
  fval = feval(params.fun,x,params.args{:});
  
  % fminsearchbnd was not called
  exitflag = 0;
  
  output.iterations = 0;
  output.funcCount = 1;
  output.algorithm = 'fminsearch';
  output.message = 'All variables were held fixed by the applied bounds';
  
  % return with no call at all to fminsearch
  return
end
% Check for an outputfcn. If there is any, then substitute my
% own wrapper function.
if ~isempty(options.OutputFcn)
  params.OutputFcn = options.OutputFcn;
  options.OutputFcn = @outfun_wrapper;
end
% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(@intrafun,x0u,options,params);
% undo the variable transformations into the original space
x = xtransform(xu,params);
% final reshape to make sure the result has the proper shape
x = reshape(x,xsize);
% Use a nested function as the OutputFcn wrapper
  function stop = outfun_wrapper(x,varargin);
    % we need to transform x first
    xtrans = xtransform(x,params);
    
    % then call the user supplied OutputFcn
    stop = params.OutputFcn(xtrans,varargin{1:(end-1)});
    
  end
end % mainline end
% ======================================
% ========= begin subfunctions =========
% ======================================
function fval = intrafun(x,params)
% transform variables, then call original function
% transform
xtrans = xtransform(x,params);
% and call fun
fval = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});
end % sub function intrafun end
% ======================================
function xtrans = xtransform(x,params)
% converts unconstrained variables into their original domains
xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.n
  switch params.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      % just in case of any floating point problems
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = params.LB(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end
end % sub function xtransform end
