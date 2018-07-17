function mixobjective(L, x, eps = 0)
  return -sum(log.(L * x + eps))
end
          
# Input argument x specifies the initial iterate of the optimization
# algorithm. When x = -1, or when any of the entries of x are
# negative, the default setting for x is used.
#
# Since the MOSEK solver for the quadratic subproblem sometimes does
# not tolerate iterates in which most of the entries are nonzero
# ("dense" vectors), the default initial estimate for qpsubprob =
# "mosek" is a sparse vector in which only the first two entries are
# nonzero.
function mixSQP(L; x = -1,
                convtol = 1e-8, pqrtol = 1e-8, eps = 1e-4, sptol = 1e-4,
                maxiter = 20, maxqpiter = 100,
                lowrank = "svd", qpsubprob = "activeset",
                nullprior = 10,
                linesearch = false,
                verbose = false)
    
# Get the number of rows (n) and columns (k) of L.
  n = size(L,1);
  k = size(L,2);

  # If any entries of input x are negative, set the initial estimate
  # to the default setting.
  #
  # When the MOSEK solver is used, the initial estimate is set to a
  # sparse vector in which all but two of the entries are zero.
  if any(x .< 0)
    if qpsubprob == "activeset"
      x = ones(k)/k;
    elseif qpsubprob == "mosek"
      x = zeros(k);
      x[1:2] = 1/2;
    else
      error("Input \"method\" should be \"activeset\" or \"mosek\"");
    end
  end
    
  # If requested (i.e., if pqrtol > 0), compute the partial QR
  # decomposition with relative precision "tol", then retrieve the
  # permutation matrix, P. For details on the partial QR, see
  # https://github.com/klho/LowRankApprox.jl.
    
  # Start timing for low-rank approximation of L.
  if lowrank != "none" && qpsubprob == "mosek"
    error("lowrank must be \"none\" when qpsubprob = \"mosek\"");
  end
  tic();
  if lowrank == "qr"
    F = pqrfact(L, rtol=pqrtol);
    P = sparse(F[:P]);
  elseif lowrank == "svd"
    F = psvdfact(L, rtol=pqrtol);
    S = Diagonal(F[:S]);
  end
  lowranktime = toq();
    
  # Summarize the analysis here.
  if verbose
    @printf("Running SQP algorithm with the following settings:\n")
    @printf("- %d x %d data matrix\n",n,k)
    @printf("- convergence tolerance = %0.2e\n",convtol)
    @printf("- zero threshold        = %0.2e\n",sptol)
    if lowrank == "qr"
      err = maximum(abs.(F[:Q]*F[:R]*P' - L));
      @printf("- partial QR tolerance  = %0.2e\n",pqrtol)
      @printf("- partial QR max. error = %0.2e\n",err)
    elseif lowrank == "svd"
      err = maximum(abs.(F[:U]*S*F[:Vt] - L));
      @printf("- partial SVD tolerance  = %0.2e\n",pqrtol)
      @printf("- partial SVD max. error = %0.2e\n",err)
    else
      @printf("- Exact derivative computation (partial QR not used).\n")
    end
  end

  # Initialize storage for the outputs obj, gmin, nnz and nqp.
  obj      = zeros(maxiter);
  gmin     = zeros(maxiter);
  nnz      = zeros(Int,maxiter);
  nqp      = zeros(Int,maxiter);
  timing   = zeros(maxiter);
  qptiming = zeros(maxiter);
    
  # Initialize loop variables used in the loops below so that they
  # are available outside the scope of the loop.
  i     = 0;
  j     = 0;
  D     = 0;
  t     = 0;
  numls = -1;
    
  # Print the column labels for reporting the algorithm's progress.
  if verbose
    @printf("iter      objective -min(g+1)  #nz #qp #ls\n")
  end

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  for i = 1:maxiter

    # Start timing the iteration.
    tic();
      
    # Compute the gradient and Hessian, optionally using the partial
    # QR decomposition to increase the speed of these computations.
    # gradient and Hessian computation -- Rank reduction method
    if lowrank == "qr"
      D = 1./(F[:Q]*(F[:R]*(P'*x)) + eps);
      g = -P * F[:R]' * (F[:Q]'*D)/n; g[1] -= nullprior/x[1]/n;
      H = P * F[:R]' * (F[:Q]'*Diagonal(D.^2)*F[:Q])*F[:R]*P'/n + eps*eye(k); H[1,1] += nullprior/x[1]^2/n;
    elseif lowrank == "svd"
      D = 1./(F[:U]*(S*(F[:Vt]*x)) + eps);
      g = -F[:Vt]'*(S * (F[:U]'*D))/n; g[1] -= nullprior/x[1]/n;
      H = (F[:V]*S*(F[:U]'*Diagonal(D.^2)*F[:U])* S*F[:Vt])/n + eps*eye(k); H[1,1] += nullprior/x[1]^2/n;
    else
      D = 1./(L*x + eps);
      g = -L'*D/n; g[1] -= nullprior/x[1]/n;
      H = L'*Diagonal(D.^2)*L/n + eps * eye(k); H[1,1] += nullprior/x[1]^2/n;
    end

    # Report on the algorithm's progress.
    if lowrank == "qr"
      obj[i] = -sum(log.(F[:Q]*(F[:R]*(P'*x)) + eps));
    elseif lowrank == "svd"
      obj[i] = -sum(log.(F[:U]*(S*(F[:Vt]*x)) + eps));
    else
      obj[i] = mixobjective(L,x,eps);
    end
    gmin[i] = minimum(g + 1);
    nnz[i]  = length(find(x .> sptol));
    nqp[i]  = j;
    if verbose
      if i == 1
          @printf("%4d %0.8e %+0.2e %4d\n",
              i,obj[i],-gmin[i],nnz[i]);
      else
          @printf("%4d %0.8e %+0.2e %4d %3d %3d\n",
              i,obj[i],-gmin[i],nnz[i],nqp[i-1],numls);
      end
    end
      
    # Check convergence of outer loop
    if minimum(g + 1) >= -convtol
      break
    end
    
    # Solve the QP subproblem using either the active-set or
    # interior-point (MOSEK) method.
    out, qptiming[i], bytes, gctime,
    memallocs = @timed if qpsubprob == "activeset"
      y,nqp[i] = qpactiveset(x,g,H,convtol = convtol,sptol = sptol,
                      maxiter = maxqpiter);
    elseif qpsubprob == "mosek"
      y = qpmosek(x,g,H);
    end
    
    # Perform backtracking line search
    if linesearch == true
        for t = 1:10
          if lowrank == "qr"
            D_new = 1./(F[:Q]*(F[:R]*(P'*y)) + eps);
          elseif lowrank == "svd"
            D_new = 1./(F[:U]*(S*(F[:Vt]*y)) + eps);
          else
            D_new = 1./(L*y + eps);
          end
          if all(D_new .> 0)
            if sum(log.(D)) - sum(log.(D_new)) > sum((x-y) .* g) / 2
              break
            end
          end
          y = (y-x)/2 + x;
        end
        numls = t;
    else
        numls = -1;
    end

    # Update the solution to the original optimization problem.
    x = copy(y);

    # Get the elapsed time for the ith iteration.
    timing[i] = toq();
  end

  # Return: (1) the solution (after zeroing out any values below the
  # tolerance); (2) the value of the objective at each iteration; (3)
  # the minimum gradient value of the modified objective at each
  # iteration; (4) the number of nonzero entries in the vector at each
  # iteration; and (5) the number of inner iterations taken to solve
  # the QP subproblem at each outer iteration.
  x[x .< sptol] = 0;
  x             = x/sum(x);
  totaltime     = lowranktime + sum(timing[1:i]);
  if verbose
    @printf("Optimization took %d iterations and %0.4f seconds.\n",i,totaltime)
  end
  return Dict([("x",full(x)), ("totaltime",totaltime),
               ("lowranktime",lowranktime), ("obj",obj[1:i]),
               ("gmin",gmin[1:i]), ("nnz",nnz[1:i]),
               ("nqp",nqp[1:i]), ("timing",timing[1:i]),
               ("qptiming",qptiming[1:i])])
end

# Solve the QP subproblem for the mix-SQP algorithm using an active
# set method.
function qpactiveset(x, g, H; convtol = 1e-8, sptol = 1e-3, maxiter = 100)

  # Get the number of degrees of freedom in the optimization problem.
  k = length(x);
    
  # Initialize the solution to the QP subproblem.
  y      = sparse(zeros(k));
  ind    = find(x .> sptol);
  y[ind] = 1/length(ind);
  i = 0;

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  for i = 1:maxiter
        
    # Define the smaller QP subproblem.
    s   = length(ind);
    H_s = H[ind,ind];
    d   = H*y + 2*g + 1;
    d_s = d[ind];

    # Solve the smaller problem.
    p      = sparse(zeros(k));
    p_s    = -H_s\d_s;
    p[ind] = p_s;

    # Check convergence using KKT
    if norm(p_s) < convtol
            
      # Compute the Lagrange multiplier.
      z = d;
      if all(z .>= -convtol)
        break
      elseif length(ind) < k
        notind  = setdiff(1:k,ind);
        ind_min = notind[findmin(z[notind])[2]];
        ind     = sort([ind; ind_min]);
      end
    else
          
      # Find a feasible step length.
      alpha     = 1;
      alpha0    = -y[ind]./p_s;
      ind_block = find(p_s .< 0);
      alpha0    = alpha0[ind_block];
      if ~isempty(ind_block)
        v, t = findmin(alpha0);
        if v < 1

          # Blocking constraint.
          ind_block = ind[ind_block[t]]; 
          alpha     = v;
              
          # Update working set if there is a blocking constraint.
          deleteat!(ind,find(ind - ind_block .== 0));
        end
      end
          
      # Move to the new "inner loop" iterate (y) along the search
      # direction.
      y = y + alpha * p;
    end
  end

  # Return the solution to the quadratic program.
  return y, Int(i)
end

# Solve the QP subproblem for the mix-SQP algorithm using MOSEK.
function qpmosek(x, g, H)
  k      = length(g);
  y      = copy(x);
  bkx    = repmat([MSK_BK_LO],k,1)[:];
  blx    = zeros(k);
  bux    = Inf * ones(k);
  numvar = length(bkx);
  c      = 2*g + 1;
    
  maketask() do task
    appendvars(task,numvar)
    putclist(task,[1:numvar;],c)
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    Q = LowerTriangular(H); ind = (Q .> 0); a, b = findn(ind);
    putqobj(task,a,b,Q[ind])
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    optimize(task)
    y = getxx(task,MSK_SOL_ITR)
  end
  return y
end

function ash2(x,s2; nv = 10, nullprior = 10)

    n = length(x); flag = 0;
    x2_div_s2 = x.^2 ./ s2;
    s2_max = 2 * maximum(x2_div_s2 - 1);
    
    flag = 0; # to check whether ash solution is a vector of zero entries
    if s2_max <= 0
        flag = 1; # flag 1 means ash will output 0 solution
        return Dict([(:flag, flag)])
    end
    
    s2_max_round = ceil(Int,log10(s2_max));
    grid = logspace(0, s2_max_round, nv) - 1;
    
    # matrix likelihood
    s2_matrix = 1 + grid' # n by m matrix of standard deviations of convolutions
    log_lik = -x2_div_s2./s2_matrix/2 .- log.(s2_matrix)/2 - log(2*pi)/2;
    offset  = maximum(log_lik,2);
    log_lik = log_lik .- offset;
    log_lik = exp.(log_lik);
    
    # fit the model
    p = mixSQP(log_lik, nullprior = nullprior)["x"];
    
    # exploit sparsity
    ind = find(p .> 0);
    ind = unique([1;ind]); # we don't need this if nullprior > 0, since the null component never gonna be 0
    ps2 = grid[ind];       # always considers a null component
    
    # check if the ash solution is a vector of 0 entries
    if length(ind) == 1
        flag = 1; # flag 1 means ash outputs 0 solution
        return Dict([(:flag, flag)])
    end
    
    # posterior calculation
    temp = 1 .+ ps2';
    comp_post_mean = (x * ps2') ./ temp;
    comp_post_sd2 = (s2 * ps2') ./ temp;
    comp_post_mean2 = comp_post_sd2 + comp_post_mean.^2;
    comp_post_prob0 = log_lik[:,ind] .* p[ind]';
    comp_post_prob = comp_post_prob0 ./ sum(comp_post_prob0,2);
    post_mean = sum(comp_post_prob .* comp_post_mean,2);
    post_mean2 = sum(comp_post_prob .* comp_post_mean2,2);
    
    # do objective calculation
    obj = sum(log.(sum(comp_post_prob0,2)) + offset);
    loglik = -0.5 * sum(log.(2*pi*s2) + (1./s2) .* (post_mean2 - 2*x.*post_mean + x.^2));
    kl = obj - loglik;
    
    
    # return posterior first/second moments
    return Dict([
                (:pm, post_mean), (:pm2, post_mean2), (:ll, log_lik), (:pp, p),
                (:cpp, comp_post_prob), (:cpm, comp_post_mean), (:cps2, comp_post_sd2),
                (:x, x), (:s2, s2), (:grid, grid), (:ps2, ps2), 
                (:obj,obj),  (:loglik, loglik), (:kl, kl),
                (:flag, flag)
                ])

end


function update_u(X, tau, v, v2, precision_type)
    if precision_type == "columnwise_constant"
        s2 = repmat(1 ./ (tau * v2), size(X,1), 1);
    elseif precision_type == "rowwise_constant"
        s2 = 1./ (tau * sum(v2,1));
    elseif precision_type == "constant"
        s2 = repmat(1./ (tau * sum(v2,1)), size(X,1), 1);
    elseif precision_type == "elementwise"
        s2 = 1./ (tau * v2);
    end
    x = ((X .* tau) * v) .* s2; x2 = x.^2;
    u = x ./ (s2/sum(x2) + 1);
    u2 = u.^2 .+ (s2 ./ (s2/sum(x2) + 1));
    log_lik = -x2./(2*s2) - log.(s2)/2 - log(2*pi)/2;
    obj = -0.5 * sum(log.(2*pi*s2) + (1./s2) .* x2);
    loglik = -0.5 * sum(log.(2*pi*s2) + (1./s2) .* (u2 - 2*x.*u + x2));
    kl = obj - loglik;
    return Dict([
                (:u,u), (:u2,u2), (:kl, kl), (:loglik, loglik), (:obj, obj)
                ])
end

function update_v(X, tau, u, u2, precision_type, nv, nullprior; alpha = 0)
    if precision_type == "columnwise_constant"
        s2 = 1./ (tau' * sum(u2,1));
    elseif precision_type == "rowwise_constant"
        s2 = repmat(1./ (tau' * u2), size(X,2), 1);
    elseif precision_type == "constant"
        s2 = repmat(1./ (tau' * sum(u2,1)), size(X,2), 1);
    elseif precision_type == "elementwise"
        s2 = 1./ (tau' * u2);
    end
    
    x = ((X .* tau)' * u) .* s2;
    if alpha == 0
        temp = ash(x,s2, nv = nv, nullprior = nullprior);
    elseif alpha == 1
        temp = ash2(x,s2, nv = nv, nullprior = nullprior);
    else
        error("Error: \"alpha\" should be 0 or 1");
    end
    return temp
end

function update_v_group(X, tau, u, u2, precision_type, nv, nullprior; alpha = 0)
    if precision_type == "columnwise_constant"
        s2 = 1./ (tau' * sum(u2,1));
    elseif precision_type == "rowwise_constant"
        s2 = repmat(1./ (tau' * u2), size(X,2), 1);
    elseif precision_type == "constant"
        s2 = repmat(1./ (tau' * sum(u2,1)), size(X,2), 1);
    elseif precision_type == "elementwise"
        s2 = s2 = 1./ (tau' * u2);
    end
    
    x = ((X .* tau)' * u) .* s2;
    if alpha == 0
        temp = ash(x[:],s2[:], nv = nv, nullprior = nullprior);
    elseif alpha == 1
        temp = ash2(x[:],s2[:], nv = nv, nullprior = nullprior);
    else
        error("Error: \"alpha\" should be 0 or 1");
    end
    return temp
end

function update_tau(R2, precision_type)
    if precision_type == "rowwise_constant"
        tau = 1./mean(R2,2);
    elseif precision_type == "columnwise_constant"
        tau = 1./mean(R2,1);
    elseif precision_type == "constant"
        tau = 1/mean(R2);
    elseif precision_type == "elementwise"
        tau = 1./R2;
    end    
end

function update_R2(X, X2, u, u2, v, v2)
    # return u2 * v2' + X2 - 2 * (u * v') .* X
    return (X - u*v').^2 + u2 * v2' - (u.^2)*(v.^2)'
end

function spca(X; iter = 100, init = "by_svd",
              precision_type = "rowwise_constant",
              return_type = "posterior_mean",
              nv = 20, nullprior = 20,
              convtol = 1e-8,
              verbose = true)
    
    # get size
    n,p = size(X);
    
    # initialize
    if init == "by_svd"
        temp = svds(X, nsv = 1)[1];
        u = temp[:U]; u2 = u.^2;
        v = temp[:V] * temp[:S]; v2 = v.^2;
    elseif init == "by_random"
        u = randn(n); u = u/norm(u); v = X'*u; u2 = u.^2; v2 = v.^2;
    else
        error("Error")
    end
    
    # save initial values
    init_u = copy(u);
    init_v = copy(v);
    
    # get R2 and tau before proceeding
    X2 = X.^2; i = 0;
    R2 = update_R2(X, X2, u, u2, v, v2);
    tau = update_tau(R2, precision_type);
    
    # auxiliary
    temp = 0; obj = 0; klu = 0; klv = 0; loglik = 0;
    
    # loop start
    for i = 1:iter
        
        u_old = copy(u);
        
        # update u
        res = update_u(X, tau, v, v2, precision_type);
        u2 = res[:u2]/norm(res[:u])^2;
        u = res[:u]/norm(res[:u]);
        
        # update v
        temp = update_v(X, tau, u, u2, precision_type, nv, nullprior; alpha = 0);
        
        if temp[:flag] == 1
            println("terminate the algorithm since ash solution is a vector of zero entries.")
            break;
        end
        
        v = temp[:pm]; v2 = temp[:pm2];
        
        # update tau
        R2 = update_R2(X, X2, u, u2, v, v2)
        tau = update_tau(R2, precision_type)
        
        loglik = -0.5 * sum(log(2*pi) - log.(tau) .+ tau .* R2);
        klu = res[:kl];
        klv = temp[:kl];
        obj = loglik + klu + klv;
        
        diff_u = norm(u - u_old);
        if verbose
            if i >= 5
                if diff_u < convtol || i == iter
                    @printf "iter %d done: norm(v) = %0.3f and diff %0.2e\n" i norm(v) diff_u;
                    break;
                end
            end
            
            if rem(i,5) == 0
                @printf "iter %d done: norm(v) = %0.3f and diff %0.2e\n" i norm(v) diff_u;
            end
        end

    end
    
    out = Dict([
                (:u, u), (:u2, u2), (:v, v), (:v2, v2), (:tau, tau),
                (:temp, temp), (:init_u, init_u), (:init_v, init_v),
                (:loglik, loglik), (:obj, obj), (:klu, klu), (:klv, klv),
                (:iter, i), (:flag, 0)
                ])
    
    if temp[:flag] == 1
        return Dict([(:flag, 1)])
    end
    
    if return_type == "posterior_median"
        v_med = compute_posterior_median(temp);
        out[:v_med] = v_med;
    end
    
    return out
end


function compute_posterior_median(temp)
    
    # redefine
    x = temp[:x]; comp_post_prob = temp[:cpp];
    comp_post_mean = temp[:cpm]; comp_post_sd2 = temp[:cps2];
    
    ind = find(comp_post_prob[:,1] .< 0.5);
    post_median = zeros(length(x));
    for i in ind
        a = comp_post_prob[i,2:end];
        b = comp_post_mean[i,2:end];
        c = comp_post_sd2[i,2:end];
        d = comp_post_prob[i,1];
        
        # post_cdf
        function post_cdf(x)
            return a' * erfc.((x - abs.(b))./sqrt.(2*c)) + d - 0.5
        end
        post_median[i] = fzero(post_cdf,-1,1e10) * sign(x[i]);
    end
    return post_median
end


function gram_schmidt(A)
    V = copy(A); k = size(V,2);
    for i = 1:k-1
        numer = sum(V[:,i] .* V[:,i+1:k],1)
        denom = norm(V[:,i])^2;
        V[:,i+1:k] = V[:,i+1:k] - repmat(V[:,i],1,k-i) .* numer / denom;
    end
    return V
end