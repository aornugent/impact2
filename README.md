# impact2

Repository for reproducible analysis of the Pinnacle vegetation data using joint species tobit models. 

Pre-print available on EcoEvoRxiv: 

**O'Reilly-Nugent, A.**, Wandrag, E., Catford, J., Gruber, B., Driscoll, D., & Duncan, R. (2018, December 18). *Measuring impact: joint-species modelling of invaded plant communities.* https://doi.org/10.32942/osf.io/rcvt4



### Installation

There are two options:

1. Run `devtools::install_github("aornugent/impact2")` to install the development version. There are several dependencies, but sometimes it has trouble updating loaded packages. If you don't have the patience to resolve this, feel free to ignore package updates. You can then load the package with `library(impact2)`.

2. If you'd like to edit or work with the code more closely, you can clone or download the repository (top right of this page), set the working directory to `../impact2` and run `devtools::load_all('.')` to build and load the package.
(or open `impact2.Rproj` in RStudio and press `Ctrl + Shift + L`)

### Analysis

Data wrangling steps are wrapped up in `format_data()`. See the documentation `?format_data` for the possible knobs and dials.
While these steps are specific to our dataset, the Stan files in `models/` are generic. 

``` r
## load data
data(cover)

## format data (specific to Pinnacle dataset)
data_list <- format_data(model = "m1",
                         years = c(2013:2016),
                         threshold = 0.2,
                         subset = "presence",
                         lkj_prior = 25)

## examine data format
str(data_list)
                         
## run model
mod <- rstan::stan(
  file = "models/m1_interacting_species_tobit_varying_intercepts_slopes.stan",
  data = data_list,
  iter = 2000,
  chains = 3,
  init_r = 0.5,
  save_warmup = F,
  control = list(max_treedepth = 15,
                 adapt_delta = 0.8))

## summarise model
summary <- rstan::summary(mod)$summary %>%
  as.data.frame()
```

### Alternatively, these are bundled up in some higher level functions

However because the model outputs are large we do not distribute them as part of the package.

```r
## run models and save output
run_models("m1")

## previous runs can be loaded
m <- load_model("m1")

## but most functions will do this
check_models(model = "m1")

## figure numbers may not match manuscript
create_figures(model = "m1", figs = 3)
```

