library(shiny)
library(xtable)
library(htmltools)
library(utils)

pancollapse <- function(){
  tags$head(tags$link(rel="stylesheet", type="text/css", href="pancollapse.css"),
            tags$script(src="pancollapse.js"))
}

alert.create <- function(content, style="info") {
  HTML(paste0('<div class="alert alert-', style, ' alert-dismissible" role="alert">'),
    '<button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>',
    content, 
    '</div>'
  )
}

pancollapse.create <- function(title, content, class="panel-default"){
  HTML('<div class="panel ', class, '"><div class="panel-heading panel-heading-collapse"><h4 class="panel-title">',
       title,
       '<span class="pull-right"><i class="glyphicon glyphicon-chevron-up"></i></span></h4></div><div class="panel-body panel-collapse collapse">',
       content,
       '</div></div>'
  )
}

panel.create <- function(title, content) {
  HTML('<div class="panel panel-default"><div class="panel-heading"><h4 class="panel-title">',
       title,
       '</h4></div><div class="panel-body">',
       content,
       '</div></div>'
  )
}

getTable <- function(df, cbGetClass = NULL) {
  thead <- paste0('<th>', htmlEscape(names(df)), '</th>', collapse='')
  
  tbody <- rep("",nrow(df))
  
  for(row in 1:nrow(df)) {
    format <- '<td>%s</td>'

    tbody[row] <- paste0(sapply(df[row,], function(x){ sprintf('<td>%s</td>', htmlEscape(as.character(x))) }), collapse='');
    
    cls <- NULL
    if( is.function(cbGetClass) ) {
      cls <- cbGetClass(df[row,])
    } 
    
    if(is.character(cls)) {
      tbody[row] <- sprintf('<tr class="%s">%s</tr>', htmlEscape(cls), tbody[row])
    } else {
      tbody[row] <- sprintf('<tr>%s</tr>', tbody[row])
    }
  }
  
  tbody2 <- paste0(tbody, collapse='')
  
  HTML(
    '<div class="table-responsive"><table style="font-size:80%;" class="data table table-condensed"><thead><tr>', 
    thead, 
    '</tr></thead><tbody>',
    tbody2,
    '</tbody></table></div>'
  )
}

classWhenValueFun <- function(col.name, value, cls) {
  return(function(named.row){
    tmp <- named.row[[col.name]]
    if(!is.null(tmp) && !is.na(tmp) && tmp == value) cls
  })
}

readFile <- function(filename) {
  fileConnection <- file(filename, encoding="UTF-8")
  text <- readChar(fileConnection, file.info(filename)$size, useBytes = TRUE)
  Encoding(text) <- "UTF-8"
  close(fileConnection)
  text
}

loadHTML <- function(filename) {
  HTML(readFile(filename))
}