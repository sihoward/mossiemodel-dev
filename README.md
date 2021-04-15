# README #

### What is this repository for? ###

This repository contains code for an in-development R package recreating the mosquito model of - Niebuhr, 2016. Avian malaria transmission dynamics in New Zealand: investigating host and vector relationships along an elevational gradient. PhD Thesis, University of Otago.

The model is updated to use climate station data accessed from NIWA's CliFlo climate information database. Currently only one location, Musselburgh, Dunedin (CliFlo station ID  - 15752), is accessed and all models are run for that station.

CliFlo requests require a user login (see <https://cliflo.niwa.co.nz>) and users need to add these as **cliflo_usrid** and **cliflo_pwd** environment variables on their local machine (see [this page](https://support.rstudio.com/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf)).


### How do I get set up? ###

Install the appropriate version of Rtools and install the 'devtools' R package.

To install directly from the github repository run `devtools::install_github(repo = "sihoward/mossiemodel-dev", ref = "master")` in the R console. 
The model is structured as an R package.

R package dependencies are flagged during install, and are deSolve, clifro, shiny, dplyr, ggplot2, tidyr, rlang.

After installation run `example(mosqmod::runModel)` to see a worked example.

To run the included shiny app locally run `mosqmod::runMosqModApp()`.

Users on networks with restrictions may have to edit the permissions file 'cacert.pem' at 'C:/Program Files/R/R-X.0.0/library/RCurl/CurlSSL/cacert.pem' to send CliFlo requests. Alternatively, use runMosqModApp(cliflo_requests = FALSE) to supress requests when running locally.

### TO DO

* Add documentation to functions
* Vignette?
* Add datasets?
* Additional climate stations
* Manual upload of temperature time series
* Usage tracking?


