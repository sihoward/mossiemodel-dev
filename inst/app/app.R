#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(clifro)
library(readr)
library(ggplot2)
# install mosqmod package from repo before pushing to shinyapps
# devtools::install_github(repo = "sihoward/mossiemodel-dev", ref = "master")
library(mosqmod)

# set cliflo_requests to TRUE if missing
if(!exists("cliflo_requests", envir = globalenv())) cliflo_requests <- TRUE


# UI ----------------------------------------------------------------------

## UI: sidebar ----
ui.sidebar <-
    list(selectInput(inputId = "selectStation", label = "Select site",
                     choices = c("Dunedin city (Musselburgh)", "Catlins (Nugget Point)")),
         dateRangeInput("runDates", label = "Run model over date range", start = as.Date("2020-07-01")),
         numericInput(inputId = "Mfloor",label = "Minimum number of adult mosquitos (M)", value = 100),
         numericInput(inputId = "extend_days",label = "Project temperature by n days", value = 30),
         actionButton(inputId = "runModel", label = "Run model", width = '100%'),
         conditionalPanel('true',
                          numericInput(inputId = "MTD", label = "Minimum temperature for mosquito development", value = 7.783),
                          numericInput(inputId = "degdays.plasmod.devel", label = "Degree days for plasmodium development", value = 86.21),
                          numericInput(inputId = "timeout.plasmod.devel", label = "Days below MTT that reset plasmodium development", value = 30),
                          numericInput(inputId = "MTT", label = "Minimum temperature for plasmodium development", value = 12.97),
                          wellPanel(numericInput(inputId = "burnin.reps", label = "Repeat burnin sequence n times", value = 100),
                                    numericInput(inputId = "yrng", label = "set plot y scale maximum", value = NULL)))
         )

## UI: main panel ----
ui.main <-
    list(# wellPanel(
        fluidRow(# column(4, actionButton(inputId = "runModel", label = "Run model", width = '100%')),
                 column(4, conditionalPanel('false', checkboxGroupInput(inputId = "selectPopn", label = "Select population to plot",
                                              choiceNames = list("Adults","Larvae")[1],
                                              choiceValues = list("M","L")[1],
                                              selected = "M"))),
                 conditionalPanel(condition = "input.runModel > 0",
                                  column(4,
                                         downloadButton("downloadData", "Download results"),
                                         downloadButton("downloadReport", "Download report"),
                                         conditionalPanel('false',
                                         selectInput("dowloadFormat",
                                                     label = "Report format",
                                                     choices = c("pdf", "word")[1]))
                                  )
                 )

        # )
    ),
    plotOutput("popnplot"),
    plotOutput("compareyearplot")
    )

ui <-
    fluidPage(title = "MossieModel",
              # Application title
              fluidRow(img(src = "MW_LR_Landscape_lge_blk.png", width = '200px', style="float:right;width:200px"),
                       h1("MossieModel"), style = "padding-left:35px; padding-right:35px;"),
              fluidRow(
                  column(wellPanel(ui.sidebar),
                         h6(paste0("Version ", packageVersion("mosqmod"), " (", packageDate("mosqmod"), ")")),
                         h6("Developed by Manaaki Whenua - Landcare Research."),
                         h6("email: ", tags$a(href = "https://www.landcareresearch.co.nz/about-us/our-people/simon-howard", target="_blank", "Simon Howard"), "or",
                            tags$a(href = "https://www.landcareresearch.co.nz/about-us/our-people/chris-niebuhr", target="_blank", "Chris Niebuhr")),
                         width = 3),
                  column(ui.main, width = 9)
              )
    )

# Define server logic required to draw a histogram
server <- function(session, input, output) {


    # server: temperature times series ----

    ## load saved temperature series
    saved_station_temps <- mosqmod::saved_station_temps

    ## select station from saved temperature series
    temp_past <- reactive({

        # convert saved temperature series to dataframe
        dat <- data.frame(saved_station_temps)
        # subset selected temperature series
        subdat <- switch(input$selectStation,
                         "Dunedin city (Musselburgh)" = subset(dat, Station == "Dunedin, Musselburgh Ews"),
                         "Catlins (Nugget Point)" = subset(dat, Station == "Nugget Point Aws"))

        ## request CliFlo temperatures ----
        ## append selected temperature series (can fail if no rows retreived)
        try({
            showNotification(id = "cliflo_note", ui = "Making request to CliFlo climate server")

            subdat <- mosqmod::append_TempSeq(# temp_stored = read.csv("www/temp_data/Musselburgh_15752_2000-2021.csv",
                # stringsAsFactors = FALSE, na.strings = ""),
                temp_stored = subdat,
                request = cliflo_requests,
                date_append = as.Date(Sys.Date()),
                username = Sys.getenv("cliflo_usrid"), password = Sys.getenv("cliflo_pwd"))

            showNotification("CliFlo request successful", duration = 1)
            removeNotification(id = "cliflo_note")
        })

        # format sequence (fill gaps, format dates etc. )
        mosqmod::formatTempSeq(d = subdat)

    })

    # server: projected degree days and temperatures ----
    temp_seq <- reactive({

        # checks for input$extend_days
        validate(
            need(!is.na(input$extend_days), "Enter a number of days to project temperature data ahead")
        )
        if(input$extend_days < 0){
            updateNumericInput(session, inputId = "extend_days", value = 0)
        }

        ## get past temperature series ----
        dat <- subset(temp_past(), Date <= input$runDates[2])

        # calculate calendar day means
        calday_means <- mosqmod::getCalendarDayMeans(temp_hist = dat)

        # project temperature ahead using calendar day means ----
        if(input$extend_days > 0) {
            # add projected temperatures
            temp_projected <-
                mosqmod::project_TempSeq(temp_seq = dat,
                                         extend_days = input$extend_days, lookback_days = 90,
                                         calday_means = calday_means)
            temp_seq <- dplyr::bind_rows(dat, temp_projected)
        } else {
            temp_seq <- dat
        }

        return(temp_seq)
    })

    ## calendar degree day means for plasmod dev ----
    calDegDay_means_plasmod <-
      reactive({
        getCalendarDegreeDays(temp_hist = temp_past(), lookback = "-10 year",
                              developtemp = input$MTT, timeout = input$timeout.plasmod.devel)
      })

    ## plasmodium development requirements ----
    plasmod_devel_seq <-
      reactive({
      # browser()
      mosqmod::plasmod_devel(temp_seq = temp_seq(), MTT = input$MTT,
                    devel_degdays = input$degdays.plasmod.devel,
                    timeout = input$timeout.plasmod.devel, extend_days = input$extend_days,
                    cal_degday_means = calDegDay_means_plasmod())
    })


    # server: update projected days if end date > projected days --------------
    observe({
      latest_record <- max(temp_seq()$Date[temp_seq()$source %in% c("stored", "retreived")])

      # update date range if runDates doesn't equal latest retreived record
      if(latest_record != input$runDates[2]){
        updateDateRangeInput(session, inputId = "runDates", end = latest_record)
      }

    }, priority = 1)


    # server: runModel() ------------------------------------------------------

    res <- eventReactive({ input$runModel | input$MTD }, {

        validate(
            need(!is.na(input$Mfloor) & input$Mfloor > 0, "Enter a starting number of adult mosquitos > 0"),
            need(input$Mfloor <= 2710200, "Starting number of adults must be < 2710200"),
            need(!is.na(input$MTD), "Enter a minimum development temperature"),
            need(input$runDates[1] >= min(temp_seq()$Date)+365,
                     sprintf("Select a start date after %s", min(temp_seq()$Date)+364)),
            need(grepl("07-01", format(input$runDates[1], "%m-%d")),
                 sprintf("Select a 1st July start date after %s", min(temp_seq()$Date)+364))
        )


        showNotification(id = "model_update", ui = "Running model ...", duration = NULL)

        # set burn-in dates to be preceeding year
        burnin.start <- seq(input$runDates[1], length.out = 2, by = "-1 year")[2]
        burnin.end <- input$runDates[1] - 1

        # run model function
        model_results <- mosqmod::runModel(temp_seq = temp_seq(),
                          burnin.dates = seq(burnin.start, burnin.end, 1),
                          run.dates = seq(from = input$runDates[1],
                                          to = max(temp_seq()$Date),
                                          by = "1 day"),
                          M = input$Mfloor, Mfloor = input$Mfloor,
                          MTD = input$MTD,
                          burnin.reps = input$burnin.reps)

        showNotification("Model completed", duration = 1)
        removeNotification(id = "model_update")

        return(model_results)

    }, ignoreInit = TRUE)


    # server: download results .csv -------------------------------------------
    output$downloadData <- downloadHandler(
        filename = function() {
            format(Sys.time(), "model_results_%Y%m%d_%H%M%S.csv")
        },
        content = function(file) {

            selectCols <- c("Date", "Station", "Tmean", "L", "M", "source",
                            "datetime", "interpolated", "Tmean_calday")

            write.csv(res()[selectCols], file, row.names = FALSE)
        }
    )

    # server: download report -------------------------------------------
    output$downloadReport <- downloadHandler(
        filename = function() {
            ext <- c(pdf = ".pdf", word = ".docx")[input$dowloadFormat]

            paste0(format(Sys.time(), "model_results_%Y%m%d_%H%M%S"),
                   ext)
        },
        content = function(file) {

            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).

            showNotification(id = "report_update", ui = "Generating report ...", duration = NULL)

            tempReport <- file.path(tempdir(), "mm_report_template.Rmd")
            file.copy(system.file("report/mm_report_template.Rmd", package = "mosqmod"), tempReport, overwrite = TRUE)

            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the document
            # from the code in this app).

            output_format <- c(pdf = "pdf_document", word = "word_document")[input$dowloadFormat]

            rmarkdown::render(tempReport, output_file = file,
                              params = list(res = res()),
                              output_format = output_format,
                              envir = new.env(parent = globalenv()))

            showNotification("Report completed", duration = 1)
            removeNotification(id = "report_update")

        }
    )

    # server: tempplot --------------------------------------------------------
    output$tempplot <- renderPlot({

        validate(
            need(input$runDates[1] >= min(temp_seq()$Date)+365,
                 sprintf("Select a start date after %s", min(temp_seq()$Date)+364))
        )

        ggplot2::ggplot(dplyr::filter(temp_seq(), Date >= input$runDates[1] & Date <= input$runDates[2]),
                                      ggplot2::aes(y = Tmean, x = Date)) +
                            ggplot2::geom_line() +
            ggplot2::geom_point(aes(color = source))

    })

    # server: popnplot --------------------------------------------------------
    output$popnplot <- renderPlot({

      validate(
        need(input$runModel > 0, "Press 'Run model' to display model results")
      )

      req(res())

      mosqmod::plot_popn(resdf = res(),
                         selectPopn = input$selectPopn, include_temp = TRUE,
                         plasmod_devel = plasmod_devel_seq(), MTT = 12.97, MTD = input$MTD) +
        ggplot2::coord_cartesian(ylim = c(0, input$yrng))
    })


    # server: compare years ---------------------------------------------------
    output$compareyearplot <- renderPlot({
        validate(
            need(!is.null(input$selectPopn), "Select checkbox for plotting adults, larvae or both")
        )

        mosqmod::plot_popn_years(resdf = res(), selectPopn = input$selectPopn)
    })

    # server: run lines -------------------------------------------------------

    printRes <- eventReactive(eventExpr = input$runLine, {
        eval(str2expression(input$consoleIn))
    })

    output$consoleOut <- renderPrint({
        printRes()
    })
}

# Run the application
shinyApp(ui = ui, server = server)
