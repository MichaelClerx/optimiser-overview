# Optimisers

The structure of each entry is as follows:

### Method name

_An attempted one-line summary_

A literature reference.
And a link.

A bit more info if needed.

- Global or local
- Derivatives required
- Constrained or unconstrained
- Smoothness requirements
- Heuristic vs comes-with-proof-of-convergence


## Proposed classes

### Direct search

_Methods that neither copmute nor approximate derivatives; they work without creating some sort of (linear) model of the objective function._

Lewis, Torczon, Trosset (2000). Direct search: then and now. 
https://doi.org/10.1016/S0377-0427(00)00423-4

Their distinguishing feature is that they compare values of f, instead of approximating them (i.e. via a Taylor expansion).
Typically proposed as heuristics without theoretical foundation, but often a theoretical reason for why they work has been found afterwards.

Aliases:
- Zero-order methods (because they don't construct a 1st or 2nd order taylor approximation around a sample of f).

Subclass of:
- Derivative-free methods.

Subclasses (from Lews et al.):
- Pattern search methods
- Simplex methods
- Adaptive sets of search directions-

Methods:
- Nelder-Mead simplex method (see Lewis et al.)

#### Pattern search

_Search using a "pattern" or lattice of points that are independent of f_

Torczon (1997). On the convergence of pattern search algorithms.
https://doi.org/10.1137/S1052623493250780

Define some kind of lattice (integer lattice but scaled to the search space).
Methods have a step size initially set by user, but can adapt this if nothing is found: crucial different w. other methods is that this happens based on f values, not on derivatives or derivative estimates.

Subclass of:
- Direct search

Methods:
- Coordinate search
- "Evolutionary operation", of G.E.P. Box
- Hooke & Jeeves' Pattern Search
- Multidirectional search algorihm of Dennis and Torczon

##### Coordinate search

_Vary one parameter at a time with some step size, and when no more improvements can be made decrease the step size_

Torczon (1997). On the convergence of pattern search algorithms.
https://doi.org/10.1137/S1052623493250780

Aliases: 
- alternating directions
- alternating variable search
- axial relaxation
- local variation

Subclass of:
- Pattern search (see Torczon 1997)

##### Hooke & Jeeves' pattern search

_Like coordinate search, but has a "pattern step": if an iteration improves the result it switches to exploring around the new best solution, but goes back again if that fails._

Torczon (1997). On the convergence of pattern search algorithms.
https://doi.org/10.1137/S1052623493250780

Hooke, Jeeves (1961). "Direct Search" solution of numerical and statistical problems.
https://doi.org/10.1145/321062.321069

##### Multidirectional search of Dennis & Torczon



#### Simplex methods
#### Adaptive sets of search directions

## List of methods

















# First attempt

## All of this stuff should be moved upwards, eventually

## Derivative-free Direct methods (aka Search methods)
_Completely ignore the idea of f having a derivative_

### Brute-force exploration
_Look at the landscape and pick the lowest point_
- Constrained
- Uniform sampling
- Latin hypercube sampling
- Should we generalise this as sampling-based optimisation?

### Random search

- Random search (RS)
  - _Sample points on a hypersphere around current, jump if better_
  - Fixed step size random search (FSSRS)
  - Optimum step size random search (OSSRS) Schumer and Steiglitz
  - Adaptive step size random search (ASSRS) Schumer and Steiglitz
  - Optimized relative step size random search (ORSSRS) Schrack and Choit
- Random optimization
  - _Sample points from a distribution, jump if better_
  - Usually a normal distribution
  - Luus-Jaakola: uniform distribution

- Simulated annealing
  - _Jump to nearby points with a P dependent on f; jump almost always at the start (high temperature) but with a likelihood that decreases every iteration (cooling)_
  - Global method
  
- Stochastic tunneling (STUN)
  - _Search using a metropolis criterion for random jumps, and a transformation of f intended to let it 'tunnel' between minima_
  - Global method

- Parallel tempering
  - _Like simulated annealing but randomly hop between temperatures_
  - Global method

### Coordinate search (aka Coordinate descent, aka Compass search)
_Evaluate f at nearby points at distance d in each direction, if better go there and end iteration, if no better point found decrease d_
- See book: "Introduction to derivative-free optimization"
- Type of "hill climbing" algorithm
- Adaptive coordinate descent
  - _Perform CS, but adapt (transform) the coordinate vectors to obtain maximum decorrelation along coordinates_
- Stochastic coordinate descent (aka stochastic hill climbing)
  - _Sometimes go the wrong way_

### Search & poll / Pattern search
- Search: Evaluate f at a finite number of points, jump to lower if found
- Poll: If no lower f found, perform local search, with some iteratively decreasing distance parameter
- Fancy versions can use e.g. surrogate models
- Should converge to _a_ mode, provided f is smooth
- Rosenbrock's method (1960) https://doi.org/10.1093/comjnl/3.3.175
- Audett & Dennis: Pattern search filter methods (2004)
  - Generalized pattern search (GPS)
  - Generating set search (GSS)
  - Mesh adaptive search (MADS)

### Simplex methods
- Simplex = extension of line segment / triangle / tetrahedron to any dimension
- Nelder-Mead (aka downhill simplex)
  - _Maintain a simplex, at each iteration replace its worst point by a value determined by the values of the other points_
  - Many extensions to avoid getting stuck
  - Unconstrained
- Subplex
  _Nelder-Mead on "a sequence of subspaces"_
  - https://dl.acm.org/citation.cfm?id=100816


### Tabu search
_Perform local search (e.g. random or coordinate search), but allow moves to worse f if at an optimum; but remember previously visited points and disallow those (mark them taboo)._

### Dividing rectangles algorithm (DIRECT)
_Partitions the (bounded) subspace into rectangles, and works out where the optimum could be for any value of the Livschitz constant (which says how fast a function varies with its input), then subdivides all rectangles where it could possible be_
- 1993
- Global method
- requires boundaries
- https://doi.org/10.1007/BF00941892
- Locally biased DIRECT (DIRECT-L)
    - https://doi.org/10.1023/A:1017930332101
    - 2001

### Powell's conjugate direction method
- Brent's method (aka Brent-Dekker method)
  - aka Principal axis (PRAXIS)
  - Book: Richard Brent, Algorithms for minimization without derivatives, 1972
- Golden-section line search
  - _Define a range around the optimum and then disect it using golden ratio cuts_




## Multi-start methods
_Start from several points to avoid local optima_

### Multi-level single-linkage (MLSL)
_Multi-start from uniform places densely in the search space should find all optima, but is expensive, so use clustering methods to try and identify clusters of starting points with the same attractor_
- https://doi.org/10.1007/BF02592070
- 1987
- Requires boundaries
- Global method
- Hybrid method: Can use any local optimiser inside!



## Basin-hopping
_Repeatedly perform and then perturb a local search_
- https://doi.org/10.1155/2012/674832




## Derivative-free evolutionary methods and metaheuristics
_Maintain some population of points that gets updated at each iteration_

### Genetic algorithms (GA)
_Generate new points using crossover and mutation, estimate the fitness of all points in the population, remove the points with the worst scores_
- Note this depends on a ranking of fitness, so it may be possible to avoid evaluating f(x) if something that gives the same rankings is available

### Differential evolution
_Loop over each point in a population, generate a new proposal for each point (as in GA), and replace point with that proposal if its f is lower_
- Differential evolution
- Self-adaptive DE (jDE, iDE, and pDE)

### Evolution strategies
_Like GA, but mutate by adding a normally distributed random value to each particle - no crossover_
- Can use ranking instead of objective function (see GA)

- Stochastic ranking for constrained evolutionary optimisation
  - https://doi.org/10.1109/4235.873238
  - 2000
  - Global method
- Improved stochastic ranking evolution strategy (ISRES)
  - https://doi.org/10.1109/TSMCC.2004.841906
  - 2005
  - Global method
- (N+1)-ES Simple Evolutionary Algorithm

### Controlled random search (CRS)
_Evolve a random population using a Heuristic that's a bit like a randomised Nelder-Mead iteration_
- Requires constraints
- https://doi.org/10.1093/comjnl/20.4.367
- 1977
- CRS2 with local mutation
    - https://doi.org/10.1007/s10957-006-9101-0
    - 2006


### Swarm algorithms

- Particle-swarm optimisation
  - _Several particles buzz around the search space randomly, but with a bias towards the previous personal and global best_
  - Unconstrained / Constrained-with-caveats
  - Particle-based method
  - Evolutionary / biologically inspired method
  - No proof of convergence
  - Global method 

- Ant colony optimization (see wikipedia)
- Artificial bee colony algorithm (see wikipedia)
- Artifical swarm intelligence
- Bees algorihtm (see wikipedia)
- Cuckoo search (see wikipedia)
- Glowworm swarm optimization (see wikipedia)
- Particle swarm optimization generational (GPSO)

### Metaheuristics

- Adaptive dimensional search (see wikipedia)
- Bat algorithm (see wikipedia)
- Colliding bodies optimisation (see wikipedia)
- Cuttlefish optimization (see wikipedia)
- Duelist algorithm (see wikipedia)
- Flower polination algorithm (see wikipedia)
- Firefly algorithm (see wikipedia)
- Gaussian adaptation (see wikipedia)
  - Apparently "based on information theory"
- Gravitational search algorithm (see wikipedia)
- Harmony search (see wikipedia)
  - Improved Harmony Search
- Hunting search (see wikipedia)
- Hydrological cycle algorithm (see wikipedia)
- Imperialist competitive algorithm (see wikipedia)
- Intelligent water drops algorithm (see wikipedia)
- Killer whale algorithm (see wikipedia)
- Mass and energy balances algorithm (see wikipedia)
- Memetic algorithm (see wikipedia)
- Rain water algorithm (see wikipedia)
- River formation dyanmics (see wikipedia)
- Spiral optiomization (see wikipedia)



## Bayesian optimistation
 _Treat f as an unknown distribution, define a prior over it, generate samples from f, combine with prior to form posterior, repeat._



## Gradient-estimating methods
_Assume f' exists, and then approximate it_

### Finite difference methods
_Approximate f' with finite differences, and find its root_

- Quasi-newton methods
  - "Iteratively guess the root is where the tangent hits the x-axis"
  - Unconstrained

### Simplex gradient methods / Implicit filtering
"Line-search using the simplex gradient"

### Natural evolution strategies (NES)
_Use a search population to estimate the "natural gradient", a gradient which takes different scaling of the coordinates into account_

- Covariance matrix adaptation evolution strategy (CMA-ES)
- Exponential natural evolution strategy (xNES)
- Separable natural evolution strategy (SNES)



## Surrogate-model methods
_Replace f by some function g for which g' is easy to find_

### Trust-region methods
_Select a region of search space, approximate f (typically with a quadratic function), and move to its minimum_

- Powell's Constrained Optimisation By Linear approximations COBYLA (constrained), linear function
- Powell's UOBYQA (unconstrained), quadratic function
- Powell's NEWUOA (unconstrained), quadratic function
  - Succeeeded by BOBYQA
- Powell's BOBYQA (constrained), quadratic function
- Powell's LINCOA (constrained), quadratic function
- Trust region reflective method (TRR)
- Dog-leg trust-region algorithm

### Data-based Online Nonlinear Extremumseeker (DONE)
_Use Fourier-based model, then optimize on that with L-BFGS_

### Root finding methods
_Find a point where f'(x) = 0 (and hope there's only one)_

- Bisection
- Newton's method
  - _Iteratively guess that the root is where the tangent hits the x-axis
  - Unconstrained
  - Truncated Newton methods
    - https://doi.org/10.1007/BF02592055
- Secant method (finite diff version of Newton's method)
- Steffensen's method
- Inverse interpolation
- Brent's method

### Line search based methods

    1. Esimate a descent direction `p_k`
    2. Esimate a step size `alpha_k`
    3. Update `x_{k+1} = x_k + alpha_k * p_k`

Choosing an `alpha_k` then requires solving a 1-dimensional sub-optimisation, leading to a split into _exact_ and _inexact_ methods.
An example of an exact method is the conjugate gradient method.
Inexact methods include "backtracking line search" and methods that try to choose an alpha that satisfies the "Wolfe conditions"

- Broyden-Fletcher-Goldfarb-Shanno (BFGS)
  - _Line-search-based method that chooses search function using an estimated Hessian, which gets updated at each step_
  - Unconstrained
  - BFGS-B is constrained
  - Limited-memory BFGS (L-BFGS or LM-BFGS)
    - _Like BFGS but with linear memory requirement, so good for high numbers of variables_


### Gradient descent (aka Steepest descent)

- Unconstrained / Constrained-with-projection
- Basic gradient descent
  - _Walk down the steepest direction_
  - Step size is fraction of gradient: sensitive parameter!
- Conjugate gradient descent
- Stochastic gradient descent (SGD)
  - Approximate or partially evaluate the gradient
  - For example, calculate the gradient in each direction separately, and perform a separate step in each direction
  - Implicit updates (ISGD)
    - Use gradient information _after_ the step
  - Adaptive subgradient methods (AdaGrad)
    - http://www.jmlr.org/papers/v12/duchi11a.html
  - RMSProp
    - Adaptive Moment Estimation (Adam) https://arxiv.org/abs/1412.6980
  - Natural gradient stochastic descent
    - Kalman-based stochastic gradient descent https://arxiv.org/abs/1512.01139
  - Momentum
    - Linear combination of last update + gradient update
- Levenberg-Marquardt algorithm
    - For nonlinear least-squares problems
    - Alternates between Gauss-Newton root-finding method, and gradient descent



### Continuation
_Define a class of functions F(x,k), so that f'(x)=F(x,1), and some solution F(x*, 0) = 0 is known, then 'continue' k to find the x where f'(x)=0_

### Others
- Conservative convex separable approximation (CCSA)
  - https://doi.org/10.1137/S1052623499362822
  - 2006
  - Must have inequality constraints
  - Can deal with very high-dimensional problems
  - Some relation to Sequential Quadratic Programming (SQP)
  - Based on Method of Moving Asymptotes (MMA)

- Shifted limited-memory variable metrics methods
  - https://doi.org/10.1016/j.cam.2005.02.010
  - 2006


## Methods requiring the 2nd-order gradient

?





## Good resources for optimisers:

- [NLopt](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/)
- [Pagmo/Pygmo](https://esa.github.io/pagmo2/docs/algorithm_list.html)
- [Scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize)
