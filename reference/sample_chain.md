# Sample a Markov chain

Sample a Markov chain using Metropolis-Hastings kernel with a
user-specified target distribution and proposal (defaulting to Barker
proposal), optionally adapting proposal parameters in a warm-up stage.

## Usage

``` r
sample_chain(
  target_distribution,
  initial_state,
  n_warm_up_iteration,
  n_main_iteration,
  proposal = barker_proposal(),
  adapters = list(scale_adapter(), shape_adapter()),
  show_progress_bar = TRUE,
  trace_warm_up = FALSE
)
```

## Arguments

- target_distribution:

  Target stationary distribution for chain. One of:

  - A one-sided formula specifying expression for log density of target
    distribution which will be passed to
    [`target_distribution_from_log_density_formula()`](http://github-pages.ucl.ac.uk/rmcmc/reference/target_distribution_from_log_density_formula.md)
    to construct functions to evaluate log density and its gradient
    using [`deriv()`](https://rdrr.io/r/stats/deriv.html).

  - A
    [`bridgestan::StanModel`](https://rdrr.io/pkg/bridgestan/man/StanModel.html)
    instance (requires `bridgestan` to be installed) specifying target
    model and data. Will be passed to
    [`target_distribution_from_stan_model()`](http://github-pages.ucl.ac.uk/rmcmc/reference/target_distribution_from_stan_model.md)
    using default values for optional arguments - to override call
    [`target_distribution_from_stan_model()`](http://github-pages.ucl.ac.uk/rmcmc/reference/target_distribution_from_stan_model.md)
    directly and pass the returned list as the `target_distribution`
    argument here.

  - A list with named entries `log_density` and `gradient_log_density`
    corresponding to respectively functions for evaluating the logarithm
    of the (potentially unnormalized) density of the target distribution
    and its gradient (only required for gradient-based proposals). As an
    alternative to `gradient_log_density` an entry
    `value_and_gradient_log_density` may instead be provided which is a
    function returning both the value and gradient of the logarithm of
    the (unnormalized) density of the target distribution as a list
    under the names `value` and `gradient` respectively. The list may
    also contain a named entry `trace_function`, correspond to a
    function which given current chain state outputs a named vector or
    list of variables to trace on each main (non-adaptive) chain
    iteration. If a `trace_function` entry is not specified, then the
    default behaviour is to trace the position component of the chain
    state along with the log density of the target distribution.

- initial_state:

  Initial chain state. Either a vector specifying just the position
  component of the chain state or a list output by `chain_state`
  specifying the full chain state.

- n_warm_up_iteration:

  Number of warm-up (adaptive) chain iterations to run.

- n_main_iteration:

  Number of main (non-adaptive) chain iterations to run.

- proposal:

  Proposal distribution object. Defaults to Barker proposal, that is the
  output of
  [`barker_proposal()`](http://github-pages.ucl.ac.uk/rmcmc/reference/barker_proposal.md).
  Proposal objects are lists which must minimally define entries
  `sample`, a function to generate sample from proposal distribution
  given current chain state and `log_density_ratio`, a function to
  compute log density ratio for proposal for a given pair of current and
  proposed chain states. If adapters are being used to adaptively tune
  the proposal scale and shape parameters, which is the default
  behaviour of `sample_chain`, then additionally the list must also
  define entries: `update` a function for updating parameters of
  proposal, `parameters` a function for getting current proposal
  parameter values, `default_target_accept_prob` a function for getting
  proposal specific default target acceptance probability for scale
  adaptation and `default_initial_scale` a function for getting proposal
  and dimension dependent default initial value for scale parameter.

- adapters:

  List of adapters to tune proposal parameters during warm-up. Defaults
  to using list with instances of
  [`scale_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/scale_adapter.md)
  and
  [`shape_adapter()`](http://github-pages.ucl.ac.uk/rmcmc/reference/shape_adapter.md),
  corresponding to respectively, adapting the scale to coerce the
  average acceptance rate to a target value using a dual-averaging
  algorithm, and adapting the shape to an estimate of the covariance of
  the target distribution.

- show_progress_bar:

  Whether to show progress bars during sampling. Requires `progress`
  package to be installed to have an effect.

- trace_warm_up:

  Whether to record chain traces and adaptation / transition statistics
  during (adaptive) warm-up iterations in addition to (non-adaptive)
  main chain iterations.

## Value

A list with entries

- `final_state`: the final chain state,

- `traces`: a matrix with named columns contained traced variables for
  each main chain iteration, with variables along columns and iterations
  along rows.

- `statistics`: a matrix with named columns containing transition
  statistics for each main chain iteration, with statistics along
  columns and iterations along rows.

- `warm_up_traces`: a matrix with named columns contained traced
  variables for each warm-up chain iteration, with variables along
  columns and iterations along rows. Only present if
  `trace_warm_up = TRUE`.

- `warm_up_statistics`: a matrix with named columns containing
  adaptation and transition statistics for each warm-up chain iteration,
  with statistics along columns and iterations along rows. Only present
  if `trace_warm_up = TRUE`.

## Examples

``` r
target_distribution <- list(
  log_density = function(x) -sum(x^2) / 2,
  gradient_log_density = function(x) -x
)
withr::with_seed(876287L, {
  results <- sample_chain(
    target_distribution,
    initial_state = stats::rnorm(2),
    n_warm_up_iteration = 1000,
    n_main_iteration = 1000
  )
})
```
