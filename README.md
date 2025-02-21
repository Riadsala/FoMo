# Model Versions

More details will be made available in the pre-print (once we finish it)

| Version Number             | Overview                           |
|----------------------------|------------------------------------|
| v0.0                       | prior predictions for 1.0          |
| v1.0                       | Clarke, Hunt & Hughes (2022)       |
| v1.1                       | As v1.0, but with rel distance     |
| v1.2                       | As v1.0, but removes rel direction |
| v1.3                       | As v1.2, but with abs direction <br> (four components) |
| v1.4                | As v1.2, but with abs direction <br> (eight components) |
| v1.5               | As v1.4, but adds rel direction back in                          |




# Import Data

Data should be imported using the import_data() function. You will need to write your own code if you have a new dataset in a new format 

The output of this import function should return a list of tibbles: d$stim and d$found

(x, y) coordinates are in Euclidean coordinates, with (0, 0) indicating the *bottom left* corner of the display. 
Coordinates are scaled such that $x \in (0, 1)$ and $y \in (0, a)$ where $a$ is the aspect ratio.



# Running The Model

## Simulation Tests


## Real Data

Some preprocessing is required. `fomo_preprocess(dataset, model_components = c("spatial", "item_class"))` in `functions/prep_data.R` computes list required for passing data to cmdStan and saves as `.rds`.

`fit_model(dataset, fomo_ver, mode = "traintest", iter = 500)` in the `functions/fit_model.R` script is a wrapper function that carries out all the required steps to fit a model.

MORE DETAILS on the inputs etc

# Tools for Working with Fitted Models

## Posterior

## Accuracies

## Predictive Checks
