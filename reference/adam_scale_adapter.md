# Create object to adapt proposal scale to coerce average acceptance rate using the Adam optimizer of Kingma and Ba (2014).

Applies the Adam stochastic-gradient optimizer to the log of the
proposal scale, using the acceptance-rate residual
`target_accept_prob - accept_prob` as the (stochastic) gradient. This
corresponds to treating `-0.5 * (accept_prob - target_accept_prob)^2` as
the objective and is the same gradient signal used by
[`dual_averaging_scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/dual_averaging_scale_adapter.md)
and
[`stochastic_approximation_scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/stochastic_approximation_scale_adapter.md),
but plugged into the Adam update rule instead of a dual-averaging or
Robbins-Monro schedule. Follows the reference implementation in the
`walnuts` library (see references below).

## Usage

``` r
adam_scale_adapter(
  initial_scale = NULL,
  target_accept_prob = NULL,
  learning_rate = 0.05,
  beta_1 = 0.9,
  beta_2 = 0.999,
  epsilon = 1e-08,
  learn_rate_decay = 0.5
)
```

## Arguments

- initial_scale:

  Initial value to use for scale parameter. If not set explicitly a
  proposal and dimension dependent default will be used.

- target_accept_prob:

  Target value for average accept probability for chain. If not set a
  proposal dependent default will be used.

- learning_rate:

  Learning rate for Adam optimizer (the `alpha` hyperparameter in the
  original Adam paper). Controls the magnitude of the update applied to
  the log-scale on each iteration. Should be positive. Defaults to
  `0.05`, which is a practical setting that gives fast convergence of
  the acceptance-rate coercion within typical MCMC warm-up lengths; the
  original Adam paper uses `1e-3`.

- beta_1:

  Exponential decay rate for the first-moment estimate of the gradient
  (the `beta_1` hyperparameter in the original Adam paper). Should be in
  `[0, 1)`. Defaults to `0.9`.

- beta_2:

  Exponential decay rate for the second-moment (squared gradient)
  estimate (the `beta_2` hyperparameter in the original Adam paper).
  Should be in `[0, 1)`. Defaults to `0.999`.

- epsilon:

  Small positive constant added to the square root of the second-moment
  estimate in the denominator of the Adam update for numerical
  stability. Defaults to `1e-8`.

- learn_rate_decay:

  Exponent controlling the decay of the effective learning rate across
  iterations: on iteration `t` the effective learning rate is
  `learning_rate / t^learn_rate_decay`. Should be in `[0, 1]`. A value
  of `0` recovers standard Adam; the default of `0.5` matches the value
  recommended in the `walnuts` implementation.

## Value

List of functions with entries

- `initialize`, a function for initializing adapter state and proposal
  parameters at beginning of chain,

- `update` a function for updating adapter state and proposal parameters
  on each chain iteration,

- `finalize` a function for performing any final updates to adapter
  state and proposal parameters on completion of chain sampling (may be
  `NULL` if unused).

- `state` a zero-argument function for accessing current values of
  adapter state variables.

## Details

To match the stability provided by dual averaging, a learning-rate decay
`learning_rate / t^learn_rate_decay` is applied to the per-iteration
Adam step; a value of `learn_rate_decay = 0` recovers standard Adam,
while `learn_rate_decay = 0.5` (the default here) is recommended by the
`walnuts` authors for stable convergence.

## References

Kingma, D. P., & Ba, J. (2014). Adam: A method for stochastic
optimization. *arXiv preprint* arXiv:1412.6980.

Reference implementation in `walnuts`:
<https://github.com/flatironinstitute/walnuts/blob/main/include/walnuts/adam.hpp>

## Examples

``` r
proposal <- barker_proposal()
adapter <- adam_scale_adapter(
  initial_scale = 1., target_accept_prob = 0.4
)
adapter$initialize(proposal, chain_state(c(0, 0)))
```
