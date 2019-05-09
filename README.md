# impact2

Repository for reproducible analysis of the Pinnacle vegetation
data using joint species tobit models.

### Installation

Download code, set the working directory to `../impact2` and run `devtools::load_all('.')` to build and load the package.
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

### Alternatively, these are bundled up in some higher level functoins

```r
## run models
run_models()

## check models
check_models()

## create figures
create_figures()
```

