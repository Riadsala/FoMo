*Model Versions*

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




*Import Data*

Data should be imported using the import_data() function. You will need to write your own code if you have a new dataset in a new format 

The output of this import function should return a list of tibbles: d$stim and d$found

(x, y) coordinates are in Euclidean coordinates, with (0, 0) indicating the *bottom left* corner of the display. 
Coordinates are scaled such that $x \in (0, 1)$ and $y \in (0, a)$ where $a$ is the aspect ratio.



*Running Simulation Tests*


*Running on Real Data*
