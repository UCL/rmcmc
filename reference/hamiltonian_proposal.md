# Create a new Hamiltonian proposal object.

The Hamiltonian proposal augments the target distribution with normally
distributed auxiliary momenta variables and simulates the dynamics for a
Hamiltonian function corresponding to the negative logarithm of the
density of the resulting joint target distribution using a leapfrog
integrator, with the proposed new state being the forward integrate
state with momenta negated to ensure reversibility.

## Usage

``` r
hamiltonian_proposal(
  n_step,
  scale = NULL,
  shape = NULL,
  sample_auxiliary = function(state) stats::rnorm(state$dimension()),
  sample_n_step = NULL
)
```

## Arguments

- n_step:

  Number of leapfrog steps to simulate Hamiltonian dynamics for in each
  proposed move, or parameter passed to function specified by
  `sample_n_step` argument if not `NULL`.

- scale:

  Scale parameter of proposal distribution. A non-negative scalar value
  determining scale of steps proposed.

- shape:

  Shape parameter of proposal distribution. Either a vector
  corresponding to a diagonal shape matrix with per-dimension scaling
  factors, or a matrix allowing arbitrary linear transformations.

- sample_auxiliary:

  A function which samples new values for auxiliary variables
  (corresponding to a linear transform of momentum) given current chain
  state, leaving their standard normal target distribution invariant.
  Defaults to a function sampling independent standard normal random
  variates but can be used to implement alternative updates such as
  partial momentum refreshment. Function should accept a single argument
  which is passed the current chain state.

- sample_n_step:

  Optionally a function which randomly samples number of leapfrog steps
  to simulate in each proposed move from some integer-valued
  distribution, or `NULL` (the default) to use a fixed deterministic
  number of steps as specified by `n_step` argument. If a function it
  should accept a single argument which will be passed the value of
  `n_step` which can be used to specify parameter(s) of distribution to
  sample number of steps from.

## Value

Proposal object. A list with entries

- `sample`: a function to generate sample from proposal distribution
  given current chain state,

- `log_density_ratio`: a function to compute log density ratio for
  proposal for a given pair of current and proposed chain states,

- `update`: a function to update parameters of proposal,

- `parameters`: a function to return list of current parameter values.

- `default_target_accept_prob`: a function returning the default target
  acceptance rate to use for any scale adaptation.

- `default_initial_scale`: a function which given a dimension gives a
  default value to use for the initial proposal scale parameter.

## References

Duane, S., Kennedy, A. D., Pendleton, B. J., & Roweth, D. (1987). Hybrid
Monte Carlo. *Physics Letters B*, 195(2), 216-222.

Neal, R. M. (2011). MCMC Using Hamiltonian Dynamics. In *Handbook of
Markov Chain Monte Carlo* (pp. 113-162). Chapman and Hall/CRC.

## Examples

``` r
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  gradient_log_density = function(x) -x
)

# Proposal with fixed number of leapfrog steps
proposal <- hamiltonian_proposal(scale = 1., n_step = 5)
state <- chain_state(c(0., 0.))
withr::with_seed(
  876287L, proposed_state <- proposal$sample(state, target_distribution)
)
log_density_ratio <- proposal$log_density_ratio(
  state, proposed_state, target_distribution
)
proposal$update(scale = 0.5)

# Proposal with number of steps randomly sampled uniformly from 5:10
sample_uniform_int <- function(lower, upper) {
  lower + sample.int(upper - lower + 1, 1) - 1
}
proposal <- hamiltonian_proposal(
  scale = 1.,
  n_step = c(5, 10),
  sample_n_step = function(n_step) sample_uniform_int(n_step[1], n_step[2])
)
withr::with_seed(
  876287L, proposed_state <- proposal$sample(state, target_distribution)
)

# Proposal with partial momentum refreshment
partial_momentum_update <- function(state, phi = pi / 4) {
  momentum <- state$momentum()
  if (is.null(momentum)) {
    stats::rnorm(state$dimension())
  } else {
    cos(phi) * momentum + sin(phi) * stats::rnorm(length(momentum))
  }
}
proposal <- hamiltonian_proposal(
  scale = 1.,
  n_step = 1,
  sample_auxiliary = partial_momentum_update
)
withr::with_seed(
  876287L, proposed_state <- proposal$sample(state, target_distribution)
)
```
