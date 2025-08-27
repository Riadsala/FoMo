library(quarto)
library(this.path)

# set the model version & dataset you want to generate graphs for
model_ver <- "1_3"
dataset <- "kristjansson2014plos"
free_or_fixed <- "f" # "v" for free, "f" for fixed

# some naming of things
filename <- paste0("model_plots_", dataset, "_", free_or_fixed, model_ver, ".html")


# render output file
quarto_render("helper/make_model_plots.qmd", output_file = filename, execute_params = list(model_ver = model_ver, 
                                                                                      dataset = dataset,
                                                                                      free_or_fixed = free_or_fixed),
              execute_dir = getwd())

# move generated file to a better folder
setwd('helper')
stopifnot(file.exists(filename))
# move filename to path.out
path.out <- '../reports';
file.rename(filename, file.path(path.out, filename))
# set working directory back to here, in case you want to repeat
setwd(this.path::here())
