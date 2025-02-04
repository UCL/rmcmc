## Resubmission

* Enclosed reference to library name, _BridgeStan_, in quotes in description
  field in DESCRIPTION file.
  
* Updated plotting function in `barker-proposal` vignette to restore values of
  graphical parameters set using `par()` during function to previous values on
  exit.

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.

* Suggests or Enhances not in mainstream repositories:
    bridgestan

  The package suggests the non-CRAN package 
  [bridgestan](https://roualdes.github.io/bridgestan/latest/languages/r.html)
  as an optional dependency. Links to the bridgestan installation instructions
  and associated GitHub repository are provided in the DESCRIPTION file. 
  Functions using the bridgestan API check for presence of the package using 
  `rlang::check_installed()`, the tests for  the related functionality are set
  to be skipped when running on CRAN and the  corresponding examples are only
  run if bridgestan is installed. bridgestan is installed and the corresponding
  tests and examples run, as part of the continuous integration runs on the
  package's GitHub repository.
