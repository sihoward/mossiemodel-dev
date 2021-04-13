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


# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("MossieModel - dev"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(dateRangeInput("runDates", label = "Run model over date range", start = as.Date("2020-07-01")),
                     numericInput(inputId = "M",label = "Starting adult mosquitos (M)", value = 100),
                     numericInput(inputId = "extend_days",label = "Project temperature by n days", value = 30),
                     actionButton(inputId = "runModel", label = "Run model"),
                     checkboxGroupInput(inputId = "selectPopn", label = "Select population to plot",
                                        choiceNames = list("Adults","Larvae"),
                                        choiceValues = list("M","L"),
                                        selected = "M")),

        # Show a plot of the generated distribution
        mainPanel(
            list(plotOutput("popnplot"),
                 plotOutput("tempplot") #,
                 # wellPanel(
                 #     h4("Run R commands"),
                 #     fluidRow(column(textInput(inputId = "consoleIn", label = "consoleIn", value = "getwd()"), width = 6),
                 #              column(actionButton(inputId = "runLine", label = "runLine"), width = 6)),
                 #     verbatimTextOutput("consoleOut"),
                 #     verbatimTextOutput("reactOut")
                 # )

            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(session, input, output) {


    # load stored temperatures

    # progress <- shiny::Progress$new(min = 0, max = 5)
    temp_stored <- mosqmod::append_TempSeq(temp_stored = read.csv("www/temp_data/Musselburgh_15752_2000-2021.csv",
                                                                  stringsAsFactors = FALSE, na.strings = ""),
                                           date_append = as.Date(Sys.Date()),
                                           username = Sys.getenv("cliflo_usrid"), password = Sys.getenv("cliflo_pwd"))
    # progress$inc(amount = 0.1)
    # progress$close()


    # format sequence (fill gaps, format dates etc. )
    temp_stored <- mosqmod::formatTempSeq(d = temp_stored)[c("Date", "Tmean", "calday", "source")]

    print(tail(temp_stored))

    # get projected temperatures from calendar day mean temperatures
    temp_seq <- reactive({

        updateDateRangeInput(session, inputId = "runDates", end = max(temp_seq$Date))
        # project temperature ahead using calendar day means
        if(input$extend_days > 0) {
            # add projected temperatures
            temp_projected <-
                mosqmod::project_TempSeq(temp_seq = temp_stored,
                                         extend_days = input$extend_days, lookback_days = 90,
                                         calday_means = read.csv("www/temp_data/calday_means.csv"))
            temp_seq <- dplyr::bind_rows(temp_stored, temp_projected)
        } else {
            temp_seq <- temp_stored
        }

        return(temp_seq)
    })



    # server: update projected days if end date > projected days --------------
    observeEvent(input$runDates, {
        if(input$runDates[2] > max(temp_seq()$Date)){
            updateNumericInput(session, inputId = "extend_days", value = input$extend_days + as.numeric(input$runDates[2] - max(temp_seq()$Date)))
        }
    })

    observeEvent(input$extend_days, {
        updateDateRangeInput(session, inputId = "runDates", end = max(temp_seq()$Date))
    })




    res <- eventReactive(input$runModel, {

        validate(
            need(!is.na(input$M) & input$M > 0, "Enter a starting number of adult mosquitos > 0"),
            need(input$M <= 2710200, "Starting number of adults must be < 2710200")
        )

        mosqmod::runModel(temp_seq = temp_seq(),
                          burnin.dates = seq(input$runDates[1] - 365, input$runDates[1] - 1, 1),
                          run.dates = seq(input$runDates[1],
                                          input$runDates[2], 1),
                          M = input$M)
        })

    # server: popnplot --------------------------------------------------------
    output$tempplot <- renderPlot({
        print(tail(temp_seq()))

        validate(
            need(input$runDates[1] >= min(temp_seq()$Date)+365,
                 sprintf("Select a start date after %s", min(temp_seq()$Date)+364))
        )

        ggplot2::ggplot(dplyr::filter(temp_seq(), Date >= input$runDates[1] & Date <= input$runDates[2]),
                                      ggplot2::aes(y = Tmean, x = Date)) +
                            ggplot2::geom_line() +
            ggplot2::geom_point(aes(color = source))

    })

    output$popnplot <- renderPlot({
        validate(
            need(!is.null(input$selectPopn), "Select checkbox for plotting adults, larvae or both")
        )

        mosqmod::plotModOut(resdf = res(), selectPopn = input$selectPopn)
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
