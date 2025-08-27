If you are interested in using this model, please feel free to contact us (Alasdair Clarke and/or Anna Hughes) for help and advice. We're happy to help. 


# Model Versions

More details are available in the pre-print.

| Version Number             | Overview                           | item class | proximity ($\delta$) | rel dir ($\psi$)  | abs dir ($\phi$) |
|----------------------------|------------------------------------|-----------------|--------------------|-----------|-----------|
| v0.0                       | prior predictions for 1.0          | $b_a$ and $b_s$ | abs | x | | 
| v1.0                       | Clarke, Hunt & Hughes (2022)       | $b_a$ and $b_s$ | abs | x |  |
| v1.1                       | Does relative prox help?    | $b_a$ and $b_s$ | rel | x | | 
| v1.2                       | Is rel direction doing anything? | $b_a$ and $b_s$ | abs |       |   |
| v1.3                       | Four comp directions | $b_a$ and $b_s$ | abs | x|  4 | 
| v1.4                       | Eight comp directions | $b_a$ and $b_s$ | abs | x |  8 | 
| v1.5                       | Is rel direction doing anything? | $b_a$ and $b_s$ | abs |  |  4 |


# Import Data

Data should be imported using the `import_data()` function: if you use `fomo_preprocess()` (see below) it will attempt to use this function to import the data. If you have a new dataset in a new format, you will need to write your own code to get it into the correct format.

The output of this import function should return a list of tibbles: d$stim and d$found.

(x, y) coordinates are in Euclidean coordinates, with (0, 0) indicating the *bottom left* corner of the display. 
Coordinates are scaled such that $x \in (0, 1)$ and $y \in (0, a)$ where $a$ is the aspect ratio.



# Running The Model

## Simulation Tests

See the `tests` folder for a range of options for simulating foraging data (`1_test_data_simulation`), testing simulated data using models without a random effects structure (`2_test_simple_models`) and testing simulated data using models with a random effects structure (`3_test_multilevel_models`). The model fitting procedure is the same as for real data, below.

## Real Data


Some preprocessing is required. `fomo_preprocess(dataset)` in `functions/prep_data.R` computes list required for passing data to cmdStan and saves as `.rds`. (This function also takes a `d0` input which is a scaling factor for `delta` to keep it on approximately the same scale as the other parameters: this has a default value of 20 which should work well for most screen-based tasks, but get in touch if you need further assistance with your particular paradigm). See `examples/1_fit_models/1_do_all_preprocess.R` for an example script.

`fit_model(dataset, fomo_ver, mode = "traintest", iter = 500)` in the `functions/fit_model.R` script is a wrapper function that carries out all the required steps to fit a model. Setting `mode = "traintest"` will mean that training and testing take part on separate partitions of the data: `mode = "all"` means that training and testing will be done on all the data. (This function also takes a `kappa` input which determines the value of the hyperparameter `kappa`: this has a default value of 20 which we anticipate will work well in most cases). See `examples/1_fit_models/2_fit_all_models.R` for an example script.

In most cases, you will also want to generate predictions from the model. This uses the `gen_quant()` function which takes the same arguments as `fit_model()`. See `examples/1_fit_models/3_calc_gen_quant.R` for an example script.

Model diagnostics are a work in progress: you can see a current example script at `examples/1_fit_models/4_check_diagnostics.R`.


# Tools for Working with Fitted Models

## Posterior

`extract_post(model, dataset)` will extract posterior density estimates. We provide a number of functions that then use the output of this function for plotting: `plot_model_fixed()` to plot fixed effects (group averages), `plot_model_random()` to plot random effects and `plot_model_theta()` to plot theta values (where relevant).

## Accuracies

`extract_pred(dataset, fomo_ver, folder, mode = "traintest")` computes the item-by-item accuracy from the model's generated quantities and joins it with the empirical data so that we can measure accuracy. See `examples/1_fit_models/5_extract_pred.R` for an example script.

There are a number of plotting functions to look at accuracies: `plot_model_accuracy()` (to plot model accuracy throughout a trial), `plot_models_accuracy()` (to plot average accuracy, for all model versions) and `plot_model_accuracy_comparison()` (to compare per person accuracies for two different model versions).

## Posterior predictions

You can calculate a range of summary statistics (either on empirical data, or model predictions). `get_run_info_over_trials()` calculates maximum run length, the number of runs, best r, PAO and the number of intersections. `get_iisv_over_trials()` returns inter-item selection vectors: delta (distance), theta (direction) and psi (relative direction).  See `examples/1_fit_models/6_compute_summary_stats.R` for an example script.

As an extra step, you can use the IISV statistics to calculate a levy flight statistic: see `examples/1_fit_models/7_compute_levy_from_iisv.R` for an example script.
