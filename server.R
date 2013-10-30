library(shiny)
library(shinyAce)
library(psych)
library(ltm)
library(CTT)



shinyServer(function(input, output) {


    options(warn=-1)
    
    
    check <- reactive({
        if (input$colname == 0) {
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            dat <- score(x, ans, output.scored=TRUE)$scored
            
        } else {
            
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            dat <- score(x, ans, output.scored=TRUE)$scored
        }
    })
    
    
    
    bs <- reactive({
        if (input$colname == 0) {
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            dat <- score(x, ans, output.scored=TRUE)$scored
            
            total <- rowSums(dat, na.rm=T)
            result <- describe(total)[2:13]
            row.names(result) <- "Total   "
            result
        
        } else {
            
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            dat <- score(x, ans, output.scored=TRUE)$scored
  
            total <- rowSums(dat, na.rm=T)
            result <- describe(total)[2:13]
            row.names(result) <- "Total   "
            result
        }
    })
    
    
    
    alpha.result <- reactive({
        if (input$colname == 0) {
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            dat <- score(x, ans, output.scored=TRUE)$scored

            result1 <- cronbach.alpha(dat)
            result2 <- alpha(dat, check.keys=F)
            list(result1, result2)
        
        } else {
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            dat <- score(x, ans, output.scored=TRUE)$scored

            result1 <- cronbach.alpha(dat)
            result2 <- alpha(dat, check.keys=F)
            list(result1, result2)
        
        }
    })
    
    
    
    distractor <- reactive({
        
      if (input$type == "frequency") {
        
        if (input$colname == 0) {
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            distractor.analysis(x, ans)


        } else {
            
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            dat <- score(x, ans, output.scored=TRUE)$scored
            
            distractor.analysis(x, ans)
        }
        
      } else {
          
          if (input$colname == 0) {
              x <- read.table(text=input$text1, sep="\t")
              x <- as.matrix(x)
              
              ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
              ans <- as.character(ans)
              
              distractor.analysis(x, ans, p.table = T)
              
              
          } else {
              
              x <- read.csv(text=input$text1, sep="\t")
              
              ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
              ans <- as.character(ans)
              
              dat <- score(x, ans, output.scored=TRUE)$scored
              
              distractor.analysis(x, ans, p.table = T)
          }
      }
    })
   
    
    
    output$distPlot <- renderPlot({
        if (input$colname == 0) {
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            x <- score(x, ans, output.scored=TRUE)$scored
            
            x <- rowMeans(x, na.rm=T)


        } else {
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            x <- score(x, ans, output.scored=TRUE)$scored

            x <- rowMeans(x, na.rm=T)
        }

            simple.bincount <- function(x, breaks) {
                nx <- length(x)
                nbreaks <- length(breaks)
                counts <- integer(nbreaks - 1)
                    for (i in 1:nx) {
                        lo <- 1
                        hi <- nbreaks
                        if (breaks[lo] <= x[i] && x[i] <= breaks[hi]) {
                                while (hi - lo >= 2) {
                                new <- (hi + lo) %/% 2
                                if(x[i] > breaks[new])
                                lo <- new
                                else
                                hi <- new
                                }
                            counts[lo] <- counts[lo] + 1
                        }
                    }
            return(counts)
            }
        
        nclass <- nclass.FD(x)
        breaks <- pretty(x, nclass)
        counts <- simple.bincount(x, breaks)
        counts.max <- max(counts)
        
        h <- hist(x, na.rm= T, las=1, breaks="FD", xlab= "Red vertical line shows the mean.",
        ylim=c(0, counts.max*1.2), main="", col = "cyan")
        rug(x)
        abline(v = mean(x, na.rm=T), col = "red", lwd = 2)
        xfit <- seq(min(x, na.rm=T), max(x, na.rm=T))
        yfit <- dnorm(xfit, mean = mean(x, na.rm=T), sd = sd(x, na.rm=T))
        yfit <- yfit * diff(h$mids[1:2]) * length(x)
        lines(xfit, yfit, col = "blue", lwd = 2)

    })
    
    
    
    output$boxPlot <- renderPlot({
        if (input$colname == 0) {
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            x <- score(x, ans, output.scored=TRUE)$scored
            
            x <- rowMeans(x, na.rm=T)
            
            
        } else {
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            x <- score(x, ans, output.scored=TRUE)$scored
            
            x <- rowMeans(x, na.rm=T)
        }
        
        boxplot(x, horizontal=TRUE, xlab= "Mean and +/-1 SD are displayed in red.")
        stripchart(x, pch = 16, add = TRUE)
        points(mean(x, na.rm=T), 0.9, pch = 18, col = "red", cex = 2)
        arrows(mean(x, na.rm=T), 0.9, mean(x, na.rm=T) + sd(x, na.rm=T), length = 0.1, angle = 45, col = "red")
        arrows(mean(x, na.rm=T), 0.9, mean(x, na.rm=T) - sd(x, na.rm=T), length = 0.1, angle = 45, col = "red")
    })
    
    
    testnorm <- reactive({
        if (input$colname == 0) {
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            x <- score(x, ans, output.scored=TRUE)$scored
            
            x <- rowMeans(x, na.rm=T)
            
            
        } else {
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            x <- score(x, ans, output.scored=TRUE)$scored
            
            x <- rowMeans(x, na.rm=T)
        }

        list(ks.test(scale(x), "pnorm"), shapiro.test(x))
    })
    
    
    output$qqPlot <- renderPlot({
        if (input$colname == 0) {
            
            x <- read.table(text=input$text1, sep="\t")
            x <- as.matrix(x)
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            x <- score(x, ans, output.scored=TRUE)$scored
            
            x <- rowMeans(x, na.rm=T)
        
        } else {
            x <- read.csv(text=input$text1, sep="\t")
            
            ans <- read.delim(text=input$text2, sep="\t", fill=TRUE, header=FALSE, stringsAsFactors=FALSE)
            ans <- as.character(ans)
            
            x <- score(x, ans, output.scored=TRUE)$scored
            
            x <- rowMeans(x, na.rm=T)

        }

        qqnorm(x, las=1)
        qqline(x, col=2)
    })




    output$check <- renderTable({
        head(check(), n = 10)
    }, digits = 0)

    output$downloadData <- downloadHandler(
        filename = function() {
            paste('Data-', Sys.Date(), '.csv', sep='')
        },
        content = function(file) {
            write.csv(check(), file)
        }
    )

    output$textarea.out <- renderPrint({
        bs()
    })
    
    output$alpha.result.out <- renderPrint({
        alpha.result()
    })
    
    output$distractor.out <- renderPrint({
        distractor()
    })

    output$testnorm.out <- renderPrint({
        testnorm()
    })
    

})
