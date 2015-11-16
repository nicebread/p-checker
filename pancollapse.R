library(shiny)
library(xtable)
library(htmltools)
library(utils)

pancollapse <- function(){
  tags$head(tags$link(rel="stylesheet", type="text/css", href="pancollapse.css"),
            tags$script(src="pancollapse.js"))
}

pancollapse.create <- function(title, content){
  return(
    div(class="panel panel-primary", 
      div(class="panel-heading panel-heading-collapse",
        h3(class="panel-title", 
          title, 
          tags$span(class="pull-right", 
            tags$i(class="glyphicon glyphicon-chevron-up")
          )
        )
      ),
      div(class="panel-body panel-collapse collapse", content)
    )
  )
}

getTable <- function(data) {
  return(
    div(class="table-responsive",
      HTML(
        paste(
          utils::capture.output(
            print(
              xtable(data), 
              type = "html", 
              html.table.attributes = paste("style=\"font-size:80%;\" class=\"", htmlEscape("data table table-bordered table-condensed", TRUE), "\"", sep = "")
            )
          ), 
          collapse = "\n"
        )
      )
    )
  )
}

loadHTML <- function(filename) {
  fileConnection = file(filename, encoding="UTF-8")
  text <- readChar(fileConnection, file.info(filename)$size, useBytes = TRUE)
  Encoding(text) <- "UTF-8"
  close(fileConnection)
  
  return(HTML(text))
}