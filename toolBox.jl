using Plots
using LinearAlgebra
using LaTeXStrings
using Optim
using Printf
using Dierckx
using DSP


function gravity(n,example=1,a=0,b=1,d=0.25)
# Set up abscissas and matrix.
   dt = 1. /n;
   ds = (b-a)/n;
   t = dt *((1:n) .- 0.5);
   s = a .+ ds*((1:n) .- 0.5);
   A = dt*d ./(d^2 .+ [(S-T).^2 for S in s,T in t]).^(3/2);

 # Set up solution vector and right-hand side.
   nt = Int(round(n/3));
   nn = Int(round(n*7/8));
   x = ones(n,1);

   if example ==1
      x = sin.(pi*t) + 0.5*sin.(2*pi*t);
   elseif example == 2
      x[1:nt] = 2*ones(nt,1);
   else
      println("Example has to be 1, 2")
      return
   end

   return A, A*x,x;

end

function picard(U,s,b)
   n, = size(s);
   beta = abs.(U[:,1:n]'*b);
   eta=beta ./s;
  
   plot(1:n,s,yaxis=:log,marker=:hexagon,labels=L"\sigma_i",xlabel='i',title="Picard plot")
   plot!((1:n)[beta.!=0],beta[beta.!=0],yaxis=:log,marker=:cross,labels=L"|u_i^Tb|")
   display(plot!((1:n)[eta.!=0],eta[eta.!=0],yaxis=:log,marker=:circle,labels=L"|u_i^Tb|/\sigma_i"))
end


function tikhonov(U,s,V,b,lambda,x_0)
#TIKHONOV Tikhonov regularization.

# Computes the Tikhonov regularized solution x_lambda, given the SVD or
# GSVD as computed via csvd or cgsvd, respectively.  If the SVD is used,
# i.e. if U, s, and V are specified, then standard-form regularization
# is applied:
#   min { || A x - b ||^2 + lambda^2 || x - x_0 ||^2 } .

#If an initial estimate x_0 is not specified, then x_0 = 0 is used.
#
# If lambda is a vector, then x_lambda is a matrix such that
#    x_lambda = [ x_lambda[1], x_lambda[2], ... ] .
#
# The solution norm (standard-form case) or seminorm (general-form
# case) and the residual norm are returned in eta and rho.

# From Matlab code of Per Christian Hansen, DTU Compute, April 14, 2003.

# Reference: A. N. Tikhonov & V. Y. Arsenin, "Solutions of Ill-Posed
# Problems", Wiley, 1977.

# Initialization.
    if (minimum(lambda)<0)
      error("Illegal regularization parameter lambda")
    end

    m = size(U,1);
    n = size(V,1);
    p = size(s,1);
    beta = U[:,1:p]' * b;
    zeta = s[:,1] .* beta;
    ll = length(lambda); 
    x_lambda = zeros(n,ll);
    rho = zeros(ll,1);
    eta = zeros(ll,1);
    omega = V'*x_0; 
    for i=1:ll
        x_lambda[:,i] = V[:,1:p]*((zeta .+ lambda[i]^2*omega)./(s.^2 .+ lambda[i]^2));
        rho[i] = lambda[i]^2*norm((beta .- s.*omega)./(s.^2 .+ lambda[i]^2));
        eta[i] = norm(x_lambda[:,i]);
    end
    if size(U,1) > p
        rho = sqrt(rho.^2 + norm(b .- U[:,1:n]*[beta;U[:,p+1:n]'*b])^2);
    end

    return x_lambda,rho,eta
end

function tsvd(U,s,V,b,k)
#SVD Truncated SVD regularization.
#
# x_k,rho,eta = tsvd(U,s,V,b,k)
#
# Computes the truncated SVD solution
#    x_k = V(:,1:k)*inv(diag(s(1:k)))*U(:,1:k)'*b .
# If k is a vector, then x_k is a matrix such that
#    x_k = [ x_k(1), x_k(2), ... ] .
# U, s, and V must be computed by the svd function.
#
# The solution and residual norms are returned in eta and rho.
    
# From Matlab code of Per Christian Hansen, DTU Compute, 12/21/97.
# Initialization.
    n,p = size(V); lk = length(k);
    if (minimum(k)<0 | maximum(k)>p)
        error("Illegal truncation parameter k")
    end
    x_k = zeros(n,lk);
    eta = zeros(lk,1); rho = zeros(lk,1);
    beta = U[:,1:p]'*b;
    xi = beta./s;
#Treat each k separately.
    for j=1:lk
        i = k[j];
        if (i>0)
            x_k[:,j] = V[:,1:i]*xi[1:i];
            eta[j] = norm(xi[1:i]);
            rho[j] = norm(beta[i+1:p]);
        end
    end

    if (size(U,1) > p)
        rho = sqrt(rho.^2 + norm(b .- U(:,1:p)*beta)^2);
    end
    return x_k,rho,eta
end



function maxk(a, k)
    b = partialsortperm(a, 1:k, rev=true)
    return [b a[b]]
end

function findrightmax(a,k)
    b=maxk(a,k)
    temp =  b[findmax(b[:,1])[2],:]
    return temp[2],Int(temp[1])
end

function detectCorner(x,y,ndiscr=100,ic=1,prop=0.02)
    lx=  log.(x[end-1:-1:2])
    ly = log.(y[end-1:-1:2])
    spl = Spline1D(lx, ly,k=2)
    dlt=(lx[end]-lx[1])/ndiscr
    xfine = lx[1]:dlt:lx[end]
    deriv1 = derivative(spl, xfine,1)
    deriv2 = derivative(spl, xfine,2)
    crit = deriv2./(abs.(deriv1) .^(3/4) .+ 0.1)
    crit=abs.(conv(ones(2ic+1),crit[ic+1:end-ic]))
    maxi,imax=findrightmax(crit,Int(round(ndiscr*prop)))
    bestRes = exp(xfine[imax])
    minsq,ibest = findmin((lx[end:-1:1].-xfine[imax]).^2)
    return bestRes,ibest[1],crit,xfine
end



function l_curve(U,s,b,method="Tikh")
#L_CURVE Plot the L-curve and find its "corner".
#
# reg_corner,rho,eta,reg_param =
#                  l_curve(U,s,b,method)

# Plots the L-shaped curve of eta, the solution norm || x || or
# semi-norm || L x ||, as a function of rho, the residual norm
# || A x - b ||, for the following methods:
#    method = 'Tikh'  : Tikhonov regularization   (solid line )
#    method = 'tsvd'  : truncated SVD or GSVD     (o markers  )
# The corresponding reg. parameters are returned in reg_param.  If no
# method is specified then 'Tikh' is default. 
#
# Note that 'Tikh', 'tsvd' require either U and s (standard-
# form regularization) computed by the function svd.
#
# The corner of the L-curve
# is identified and the corresponding reg. parameter reg_corner is
# returned.  Use routine l_corner if an upper bound on eta is required.

# Reference:  P. C. Hansen & D. P. O'Leary, "The use of the L-curve in
# the regularization of discrete ill-posed problems",  SIAM J. Sci.
# Comput. 14 (1993), pp. 1487-1503.

# From Matlab code of Per Christian Hansen, DTU Compute, October 27, 2010.
    
    eps = Base.eps()
    npoints = 200;  # Number of points on the L-curve for Tikh 
    smin_ratio = 16*eps;  # Smallest regularization paramete
    m,n = size(U); 
    p = size(s,1);
    beta = U'*b;
    beta2 = norm(b)^2 - norm(beta)^2;
    beta = beta[1:p];
    xi = beta[1:p]./s;
    xi[abs.(xi).==Inf] .= 0;
    reg_corner=0
    if method=="Tikh"
        eta = zeros(npoints,1); 
        reg_param= zeros(npoints,1); 
        rho  = zeros(npoints,1);; 
        s2 = s.^2;
        reg_param[npoints] = maximum([s[p],s[1]*smin_ratio]);
        ratio = (s[1]/reg_param[npoints])^(1/(npoints-1));
        for i=npoints-1:-1:1
            reg_param[i] = ratio*reg_param[i+1];
        end
        for i=1:npoints
            f = s2./(s2 .+ reg_param[i]^2);
            eta[i] = norm(f.*xi);
            rho[i] = norm((1. .-f).*beta[1:p]);
        end
        if (m > n) & (beta2 > 0)
            rho = sqrt(rho.^2 + beta2); 
        end
      marker = :circle; txt = "Tikh.";fmt=string("%.2e")
    end
    if method =="tsvd"
        eta = zeros(p,1); 
        rho = zeros(p,1);
        eta[1] = abs(xi[1])^2;
        for k=2:p
            eta[k] = eta[k-1] + abs(xi[k])^2; 
        end
        eta .= sqrt.(eta);
        if m > n
            if beta2 > 0
                rho[p] = beta2;
            else 
                rho[p] = eps^2; 
            end
        else
            rho[p] = eps^2;
        end
        for k=p-1:-1:1
            rho[k] = rho[k+1] + abs(beta[k+1])^2; 
        end
        rho .= sqrt.(rho);
        reg_param = (1:p)'; marker =:circle;
        U = U[:,1:p]; 
        txt = "TSVD";
        fmt=string("%.d")
    end
        
    
    np = 10;  
    n = length(rho); 
    ni = Int(round(n/np));
    bestRes,ibest,crit,xfine = detectCorner(rho,eta,100000,5) 
    a=plot(rho[2:end-1],eta[2:end-1],yaxis=:log,xaxis=:log,title="L-curve",legend=false)  
    vline!([bestRes],legend=false)
    xlabel!(L"\textrm{Residual norm } |\!|\!A x - b |\!|_2")
    ylabel!(L"\textrm{Solution norm }  |\!|\!x |\!|_2")
     
    for k = ni:ni:n
        annotate!(rho[k],eta[k], text((@eval@sprintf($fmt, $reg_param[$k])), :red, :right, 6))
    end
    display(a)
    return ibest,rho,eta,reg_param
end



function newton(lambda_0,delta,s,beta,omega,delta_0)
#NEWTON Newton iteration (utility routine for DISCREP).
#
# lambda = newton(lambda_0,delta,s,beta,omega,delta_0)
#
# Uses Newton iteration to find the solution lambda to the equation
#    || A x_lambda - b || = delta ,
# where x_lambda is the solution defined by Tikhonov regularization.
#
# The initial guess is lambda_0.
#
# The norm || A x_lambda - b || is computed via s, beta, omega and
# delta_0.  Here, s holds  the singular values of A.  Moreover,
# beta = U'*b and omega is either V'*x_0 or the first p elements of
# inv(X)*x_0.  Finally, delta_0 is the incompatibility measure.

# Reference: V. A. Morozov, "Methods for Solving Incorrectly Posed
# Problems", Springer, 1984; Chapter 26.

#From Matlab code of Per Christian Hansen, IMM, 12/29/97.

# Set defaults.
    thr = sqrt(Base.eps());  #Relative stopping criterion.
    it_max = 50;      # Max number of iterations.

# Initialization.
    if lambda_0 < 0
      print("Initial guess lambda_0 must be nonnegative")
    end
            
    p = size(s,1);
    s2 = s.^2;

# Use Newton's method to solve || b - A x ||^2 - delta^2 = 0.
# It was found experimentally, that this formulation is superior
# to the formulation || b - A x ||^(-2) - delta^(-2) = 0.

    lambda = lambda_0; 
    step = 1; 
    it = 0;
    while (abs(step) > thr*lambda) & (abs(step) > thr) & (it < it_max)
        it = it+1;
        f = s2./(s2 .+ lambda^2);
        r = (1 .-f).*(beta .- s.*omega);
        z = f.*r;
        step = (lambda/4)*((r'*r)[1] + (delta_0+delta)*(delta_0-delta))/(z'*r)[1];
        lambda = lambda - step;
  # If lambda < 0 then restart with smaller initial guess.
        if lambda < 0
            lambda = 0.5*lambda_0; 
            lambda_0 = 0.5*lambda_0;
        end
    end

# Terminate with an error if too many iterations.
    if  (abs(step) > thr*lambda) & (abs(step) > thr)
        print(["Max. number of iterations ("*string(it_max)*") reached"])
    end
    return lambda
end



function  discrep(U,s,V,b,delta,x_0)
#DISCREP Discrepancy principle criterion for choosing the reg. parameter.
#
# x_delta,lambda = discrep(U,s,V,b,delta,x_0)
#
# Least squares minimization with a quadratic inequality constraint:
#    min || x - x_0 ||       subject to   || A x - b || <= delta
#
# where x_0 is an initial guess of the solution, and delta is a
# positive constant.  Requires either the compact SVD of A saved as
# U, s, and V, 
# The regularization parameter lambda is also returned.
#
# If delta is a vector, then x_delta is a matrix such that
#    x_delta = [ x_delta(1), x_delta(2), ... ] .
#
# If x_0 is not specified, x_0 = 0 is used.
#
# Reference: V. A. Morozov, "Methods for Solving Incorrectly Posed
# Problems", Springer, 1984; Chapter 26.
#
# From Matlab code of Per Christian Hansen, IMM, August 6, 2007.
#
# Initialization.
    m = size(U,1);          n = size(V,1);
    p = size(s,1);          ld  = length(delta);
    x_delta = zeros(n,ld);  lambda = zeros(ld,1);  rho = zeros(p,1);
    if minimum(delta)<0
        print("Illegal inequality constraint delta")
    end
    x_0 = zeros(n,1)
    omega = V'*x_0

#Compute residual norms corresponding to TSVD
    beta = U'*b;

    delta_0 = norm(b - U*beta);
    rho[p] = delta_0^2;
    for i=p:-1:2
        rho[i-1] = rho[i] + (beta[i] - s[i]*omega[i])^2;
    end

# Check input.
    if (minimum(delta) < delta_0)
        print("Irrelevant delta < || (I - U*U'')*b ||")
    end

# Determine the initial guess via rho-vector, then solve the nonlinear
# equation || b - A x ||^2 - delta_0^2 = 0 via Newton's method.
    
#The standard-form case.
    s2 = s.^2;
    for k=1:ld
        if delta[k]^2 >= (norm(beta - s.*omega)^2 + delta_0^2)
            x_delta[:,k] = x_0;
        else
            dummy,kmin = findmin(abs.(rho .- delta[k]^2));
            lambda_0 = s[kmin];
            lambda[k] = newton(lambda_0,delta[k],s,beta,omega,delta_0);
            e = s./(s2 .+ lambda[k]^2); 
            f = s.*e;
            x_delta[:,k] = V[:,1:p]*(e.*beta + (1 .-f).*omega);
        end
    end
    return x_delta,lambda
end

function  gcvfun(lambda::Real,s2,beta,delta0,mn)

# Auxiliary routine for gcv.  PCH, IMM, Feb. 24, 2008.

   f = lambda^2 ./(s2 .+ lambda^2);

return (norm(f.*beta)^2 + delta0)/(mn + sum(f))^2;
        
end 

function  gcv(U,s,b,method = "Tikh")
#GCV Plot the GCV function and find its minimum.
#
# reg_min,G,reg_param = gcv(U,s,b,method)
#
#    Plots the GCV-function
#          || A*x - b ||^2
#    G = -------------------
#        (trace(I - A*A_I)^2
# as a function of the regularization parameter reg_param. Here, A_I is a
# matrix which produces the regularized solution.
#
# The following methods are allowed:
#    method = 'Tikh' : Tikhonov regularization   (solid line )
#    method = 'tsvd' : truncated SVD or GSVD     (o markers  )
#
# If method is not specified, 'Tikh' is default.  U and s,
# must be computed by the functions svd.
#
# Per Christian Hansen, DTU Compute, Dec. 16, 2003.
#
# Reference: G. Wahba, "Spline Models for Observational Data",
# SIAM, 1990.

    
    eps = Base.eps()
    npoints = 200;                      # Number of points on the curve.
    smin_ratio = 16*eps;                # Smallest regularization parameter.

# Initialization.
    m,n = size(U); 
    p  = size(s,1);
    s2 = zeros(p)
    beta = U'*b; 
    beta2 = norm(b)^2 - norm(beta)^2;

    if method == "Tikh"
      # Vector of regularization parameters.
      reg_param = zeros(npoints); 
      G = zeros(npoints); 
      s2 .= s.^2;
      reg_param[npoints] = maximum([s[p],s[1]*smin_ratio]);
      ratio = (s[1]/reg_param[npoints])^(1. /(npoints-1.));
      for i=npoints-1:-1:1
          reg_param[i] = ratio * reg_param[i+1]; 
      end

  # Intrinsic residual.
      delta0 = 0;
      if (m > n) & (beta2 > 0)
          delta0 = beta2; 
      end

  # Vector of GCV-function values.
    
      for i= 1:npoints
          G[i] = gcvfun(reg_param[i], s2, beta[1:p],delta0,m-n)
      end
      
        xlab = L"\lambda"
        ylab = L"G(\lambda)"
        fmt=string("%.2e")
  # Find minimum
     
     minG,minGi = findmin(G); # Initial guess.
     minp=reg_param[min(minGi[1]+1,npoints)]
     maxp=reg_param[max(minGi[1]-1,1)]
     gfunc(x::Real) = gcvfun(x::Real, s2, beta[1:p],delta0,m-n)
     reg_min = Optim.minimizer(optimize(gfunc,minp,maxp))
     xaxis = :log
        
    elseif method =="tsvd"
   
  # Vector of GCV-function values.
      rho2 = zeros(p-1)
      rho2[p-1] = abs(beta[p])^2;
      if (m > n) & (beta2 > 0)
          rho2[p-1] = rho2[p-1] + beta2;
      end
      for k=p-2:-1:1
          rho2[k] = rho2[k+1] + abs(beta[k+1])^2; 
      end
      G = zeros(p-1);
      for k=1:p-1
        G[k] = rho2[k]/(m - k + (n - p))^2;
      end
      reg_param = (1:p-1);
      xaxis = :none;
      minG,reg_min = findmin(G);
      xlab = L"k"
      ylab = L"G(k)"
      fmt=string("%.d")
    end
    
  # Plot GCV function.
    
  a=plot(reg_param,G,xaxis=xaxis,yaxis=:log,title = "GCV function",legend=false)
  xlabel!(xlab)
  ylabel!(ylab)
  plot!([reg_min,reg_min],[minG/2,minG],legend=false)
  annotate!(reg_min,minG, text((@eval@sprintf($fmt, $reg_min)), :red, :right, 7))
  
  display(a)
  

return reg_min,reg_param ,G
end




