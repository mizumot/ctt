library(shiny)
library(shinyAce)


shinyUI(bootstrapPage(

    headerPanel("Classical Test Theory (Item Analysis)"),

    mainPanel(
        tabsetPanel(

        tabPanel("Main",

#h3("Executing the analysis"),
#p('If you input a new dataset, click on this button to execute the analysis.'),

#submitButton("Compute"),

#br(),
#br(),

            strong('Option:'),

            checkboxInput("colname", label = "The input data includes variable names (the header) in the first row.", TRUE),

            br(),

            p('Note: Input values (either numeric or character) must be separated by tabs. Copy and paste from Excel/Numbers.'),

            aceEditor("text1", value="i01\ti02\ti03\ti04\ti05\ti06\ti07\ti08\ti09\ti10\ti11\ti12\ti13\ti14\ti15\ti16\ti17\ti18\ti19\ti20\nA\tB\tB\tB\tB\tC\tB\tC\tB\tD\tD\tC\tA\tB\tA\tD\tB\tD\tA\tC\nC\tD\tA\tD\tC\tB\tD\tB\tD\tA\tD\tD\tA\tB\tC\tC\tC\tA\tD\tC\nB\tD\tC\tD\tA\tB\tA\tC\tB\tD\tB\tA\tA\tD\tD\tA\tB\tC\tB\tB\nC\tC\tD\tD\tD\tA\tA\tD\tD\tD\tA\tB\tC\tB\tD\tB\tC\tB\tC\tA\nA\tA\tA\tD\tA\tA\tD\tB\tA\tC\tA\tD\tC\tC\tC\tC\tA\tA\tA\tB\nA\tA\tB\tC\tC\tA\tA\tA\tA\tA\tB\tC\tC\tC\tC\tB\tD\tC\tD\tD\nA\tA\tB\tA\tA\tA\tD\tB\tC\tC\tB\tC\tD\tA\tB\tD\tB\tB\tB\tD\nA\tC\tA\tD\tC\tA\tD\tA\tA\tA\tD\tD\tC\tC\tB\tA\tD\tC\tA\tD\nD\tB\tA\tD\tD\tA\tD\tB\tB\tA\tB\tB\tB\tC\tA\tA\tD\tA\tC\tB\nC\tC\tA\tC\tB\tC\tD\tC\tA\tA\tD\tD\tA\tA\tB\tC\tB\tB\tC\tC\nD\tA\tC\tB\tD\tA\tD\tB\tD\tA\tA\tD\tD\tC\tA\tC\tD\tC\tA\tD\nA\tD\tA\tD\tC\tA\tA\tA\tC\tD\tB\tB\tB\tA\tA\tC\tC\tD\tC\tC\nD\tC\tB\tA\tD\tA\tD\tB\tB\tA\tB\tD\tC\tC\tC\tC\tD\tA\tB\tC\nD\tC\tC\tD\tA\tA\tD\tB\tD\tB\tD\tD\tC\tB\tB\tB\tD\tB\tC\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tC\tC\tC\tB\tC\tB\tA\tA\tB\nB\tC\tB\tC\tB\tB\tD\tD\tA\tA\tB\tA\tC\tD\tD\tB\tA\tA\tB\tC\nD\tA\tA\tA\tD\tD\tB\tC\tB\tA\tB\tA\tD\tB\tD\tC\tC\tD\tD\tC\nC\tD\tB\tC\tC\tA\tD\tC\tC\tB\tC\tA\tC\tA\tA\tC\tC\tB\tA\tC\nD\tB\tA\tC\tD\tB\tD\tB\tD\tC\tD\tD\tA\tC\tA\tD\tC\tD\tD\tD\nA\tC\tA\tD\tB\tC\tC\tD\tD\tC\tA\tB\tC\tC\tA\tD\tB\tC\tA\tB\nD\tC\tA\tD\tB\tA\tD\tA\tA\tA\tA\tD\tC\tC\tC\tC\tB\tB\tD\tD\nB\tA\tB\tB\tB\tA\tD\tD\tD\tD\tD\tD\tB\tB\tB\tD\tB\tD\tC\tC\nB\tC\tC\tA\tC\tC\tD\tB\tD\tA\tC\tD\tD\tA\tA\tA\tC\tB\tD\tC\nA\tC\tA\tA\tA\tB\tD\tB\tA\tA\tC\tC\tC\tC\tA\tB\tA\tB\tB\tB\nD\tD\tD\tB\tC\tA\tD\tB\tC\tA\tD\tD\tB\tD\tB\tC\tB\tA\tB\tA\nA\tD\tA\tB\tA\tA\tB\tC\tB\tB\tA\tA\tB\tA\tC\tA\tD\tB\tD\tB\nB\tD\tD\tA\tC\tB\tD\tD\tA\tA\tB\tD\tC\tA\tA\tD\tD\tA\tD\tD\nD\tC\tD\tB\tA\tA\tB\tC\tC\tB\tC\tD\tC\tC\tB\tD\tA\tB\tC\tB\nD\tB\tD\tB\tD\tA\tD\tC\tD\tD\tC\tD\tC\tD\tD\tB\tC\tC\tB\tD\nB\tB\tD\tC\tD\tA\tB\tB\tB\tB\tB\tC\tA\tC\tC\tA\tB\tA\tB\tA\nD\tC\tC\tC\tC\tA\tD\tA\tB\tA\tC\tD\tC\tC\tB\tA\tA\tC\tB\tB\nB\tA\tA\tB\tD\tA\tA\tD\tD\tA\tA\tC\tD\tA\tD\tA\tA\tC\tB\tA\nB\tA\tA\tB\tB\tC\tB\tA\tC\tA\tA\tC\tD\tD\tB\tD\tA\tB\tB\tB\nD\tA\tC\tB\tB\tA\tC\tB\tB\tD\tB\tD\tC\tA\tA\tA\tC\tD\tD\tD\nB\tC\tB\tA\tD\tA\tD\tB\tA\tA\tD\tB\tC\tA\tC\tB\tA\tB\tC\tC\nA\tA\tC\tB\tA\tA\tC\tB\tB\tC\tD\tA\tB\tC\tA\tD\tA\tB\tB\tD\nD\tC\tA\tC\tA\tA\tD\tB\tC\tB\tA\tD\tC\tB\tD\tC\tD\tD\tC\tA\nD\tC\tC\tB\tA\tA\tD\tB\tC\tB\tD\tD\tB\tA\tD\tA\tD\tC\tD\tD\nB\tC\tB\tC\tA\tA\tD\tB\tA\tA\tA\tC\tD\tB\tD\tC\tB\tD\tC\tC\nD\tC\tA\tB\tD\tA\tD\tA\tB\tB\tC\tC\tC\tA\tA\tA\tC\tC\tA\tC\nC\tC\tD\tC\tB\tD\tA\tC\tD\tA\tC\tB\tA\tD\tD\tA\tB\tD\tA\tC\nC\tD\tB\tD\tA\tA\tC\tB\tC\tA\tB\tD\tC\tD\tC\tA\tD\tC\tB\tA\nC\tD\tC\tB\tB\tA\tD\tB\tD\tA\tD\tD\tC\tB\tD\tB\tB\tA\tA\tB\nB\tB\tA\tB\tB\tA\tD\tB\tC\tC\tD\tB\tC\tC\tD\tA\tD\tC\tD\tC\nD\tB\tD\tC\tC\tA\tD\tD\tB\tA\tC\tA\tC\tC\tA\tB\tC\tC\tA\tA\nD\tA\tD\tA\tA\tA\tD\tB\tD\tB\tB\tB\tA\tB\tB\tC\tC\tB\tB\tA\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tB\tD\tC\tC\tB\tA\tD\tA\tB\tB\nD\tC\tA\tD\tB\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tB\tD\tA\tC\tB\nD\tC\tD\tA\tB\tA\tD\tB\tD\tD\tA\tD\tC\tD\tB\tC\tD\tA\tA\tB\nA\tB\tA\tB\tA\tA\tD\tB\tA\tC\tC\tB\tC\tD\tB\tD\tC\tA\tD\tA\nB\tB\tA\tB\tC\tA\tD\tA\tB\tA\tC\tD\tC\tC\tC\tB\tB\tA\tD\tD\nD\tC\tC\tD\tD\tA\tD\tB\tD\tA\tC\tD\tC\tD\tD\tA\tB\tB\tA\tC\nD\tC\tD\tD\tD\tA\tD\tC\tB\tA\tA\tD\tC\tD\tB\tD\tD\tC\tB\tA\nA\tB\tA\tD\tC\tA\tB\tB\tC\tC\tC\tB\tD\tA\tA\tA\tB\tB\tC\tD\nB\tC\tA\tC\tA\tA\tD\tA\tB\tC\tD\tA\tC\tB\tC\tB\tC\tC\tB\tC\nC\tA\tB\tD\tC\tA\tD\tA\tA\tD\tC\tB\tA\tA\tB\tB\tB\tA\tC\tD\nA\tD\tA\tB\tC\tB\tD\tB\tB\tC\tB\tA\tC\tC\tA\tA\tA\tC\tC\tD\nB\tA\tA\tB\tD\tA\tD\tC\tD\tA\tA\tD\tC\tC\tD\tA\tD\tA\tA\tB\nD\tC\tA\tB\tC\tA\tD\tD\tA\tA\tD\tA\tC\tB\tC\tB\tA\tA\tA\tC\nA\tB\tD\tC\tC\tA\tD\tA\tD\tA\tD\tD\tC\tC\tA\tB\tC\tB\tB\tD\nD\tC\tA\tD\tD\tA\tD\tB\tC\tA\tA\tD\tA\tC\tB\tC\tD\tA\tD\tB\nB\tA\tA\tC\tC\tC\tD\tB\tB\tC\tA\tA\tA\tD\tB\tB\tD\tD\tD\tC\nD\tC\tA\tB\tB\tA\tD\tA\tB\tB\tA\tD\tC\tC\tA\tC\tD\tB\tA\tB\nD\tC\tA\tB\tC\tA\tD\tB\tA\tA\tA\tD\tC\tC\tB\tC\tC\tA\tD\tB\nC\tB\tA\tC\tC\tA\tB\tA\tA\tD\tB\tA\tC\tA\tA\tA\tA\tB\tC\tB\nC\tC\tD\tA\tD\tC\tB\tB\tA\tB\tC\tD\tC\tB\tC\tC\tD\tA\tA\tA\nB\tA\tB\tA\tB\tA\tA\tA\tC\tC\tC\tB\tD\tD\tB\tB\tB\tB\tA\tA\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tB\tD\tB\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tC\tD\tA\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tB\tD\tA\tA\tB\nA\tC\tB\tA\tD\tA\tA\tB\tB\tB\tA\tD\tC\tB\tB\tA\tD\tD\tA\tD\nA\tC\tA\tB\tC\tB\tD\tB\tC\tD\tA\tD\tC\tC\tA\tD\tD\tD\tB\tA\nD\tA\tA\tA\tA\tA\tD\tB\tB\tB\tA\tA\tA\tC\tB\tC\tD\tD\tA\tB\nA\tC\tC\tD\tB\tA\tC\tD\tA\tA\tA\tA\tC\tC\tB\tA\tD\tC\tB\tB\nD\tC\tA\tD\tD\tA\tD\tA\tA\tD\tD\tD\tC\tC\tC\tD\tA\tA\tA\tB\nD\tC\tD\tB\tD\tA\tC\tA\tC\tA\tA\tD\tC\tC\tD\tD\tC\tC\tB\tD\nD\tC\tA\tD\tA\tA\tD\tB\tA\tB\tA\tD\tD\tC\tB\tB\tD\tD\tD\tC\nB\tC\tA\tB\tC\tA\tD\tB\tA\tD\tA\tD\tC\tC\tA\tC\tB\tB\tA\tB\nD\tC\tB\tD\tB\tA\tD\tC\tD\tA\tB\tD\tC\tC\tB\tC\tD\tD\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tC\tD\tC\tC\tD\tA\tD\tA\tB\tD\nD\tC\tA\tB\tD\tA\tD\tB\tD\tB\tA\tD\tC\tC\tB\tD\tA\tB\tA\tD\nC\tC\tA\tC\tD\tA\tD\tB\tD\tB\tA\tD\tA\tC\tD\tC\tC\tD\tA\tB\nB\tC\tA\tD\tC\tB\tB\tD\tB\tA\tB\tA\tD\tC\tD\tC\tA\tD\tD\tB\nC\tA\tD\tA\tA\tA\tB\tC\tA\tC\tB\tC\tC\tD\tB\tB\tC\tC\tD\tB\nB\tB\tA\tC\tA\tA\tD\tB\tC\tC\tD\tA\tB\tA\tC\tD\tD\tB\tC\tA\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tD\tD\tA\tA\tC\nC\tC\tA\tB\tC\tA\tD\tB\tD\tD\tA\tD\tA\tC\tD\tA\tC\tC\tB\tB\nA\tC\tB\tD\tD\tA\tD\tB\tA\tA\tA\tD\tC\tC\tD\tC\tA\tB\tA\tD\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tA\tD\tA\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tB\tD\tA\tD\tA\tC\tC\tC\tD\tA\tC\tB\nD\tB\tA\tB\tD\tA\tA\tC\tC\tD\tC\tC\tB\tD\tC\tC\tB\tB\tD\tB\nD\tA\tD\tB\tA\tA\tA\tC\tB\tA\tC\tA\tC\tD\tB\tC\tB\tC\tD\tD\nB\tD\tA\tB\tB\tA\tC\tD\tC\tB\tB\tD\tD\tA\tB\tA\tB\tD\tC\tC\nD\tA\tB\tB\tC\tA\tD\tB\tC\tA\tA\tD\tC\tC\tB\tC\tD\tC\tD\tB\nA\tA\tB\tB\tD\tA\tA\tD\tB\tA\tB\tD\tC\tC\tD\tB\tB\tA\tB\tD\nD\tC\tA\tD\tD\tA\tD\tB\tB\tC\tB\tD\tA\tA\tC\tB\tD\tA\tD\tB\nD\tC\tC\tD\tD\tA\tD\tB\tC\tA\tA\tD\tC\tC\tB\tB\tD\tD\tB\tB\nD\tC\tA\tD\tB\tA\tD\tB\tA\tA\tA\tD\tC\tC\tB\tB\tD\tA\tA\tB\nB\tB\tA\tC\tD\tA\tD\tB\tC\tA\tD\tB\tC\tC\tD\tC\tD\tC\tB\tD\nD\tA\tA\tD\tD\tA\tD\tB\tB\tA\tA\tD\tB\tC\tC\tC\tD\tA\tA\tB",
                mode="r", theme="cobalt", height="400px"),

            p("Input answer keys (Either numeric or character, separated by tabs.):"),

            aceEditor("text2", value="D\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tC\tD\tA\tA\tB", mode="r", theme="chrome", height="50px"),

            br(),

            h3("Checking the 1-0 converted data"),
            p('Only the first 10 observations are displayed.'),
            p('If you want to download the converted data, use',
            a('Binary (1-0) Data Converter', href='https://langtest4.shinyapps.io/biconv/', target="_blank"), '.'),

            tableOutput("check"),

            br(),

            h3("Basic statistics"),
            verbatimTextOutput("textarea.out"),

            br(),

            h3("Cronbach's coefficient alpha"),
            verbatimTextOutput("alpha.result.out"),

            br(),

            h3("Item analysis"),
            verbatimTextOutput("item.analysis.out"),
            p('Item_Mean: item facility (IF)', br(),
            'I-R_Correl: Item-Remainder score correlation or "corrected item-total correlation"', br(),
            'r_pbi: Point-biserial correlation or "item-total correlation"', br(),
            'U-L_DISC: item discrimination (upper 1/3 - lower 1/3)', br(),
            'AENO: actual equivalant number of options (out of the total number of options)'),

            br(),

            h3("Distractor analysis"),

            radioButtons("type", "",
                        list("Frequency" = "frequency", "Proportion" = "proportion"), selected = "frequency"),

            verbatimTextOutput("distractor.out"),

            br(),

            h3("Histogram of the total score"),

            plotOutput("distPlot"),

            br(),

            h3("Box plot with individual data points"),

            plotOutput("boxPlot"),

            br(),

            h3("Test of normality"),
            verbatimTextOutput("testnorm.out"),

            br(),

            h3("Q-Q plot"),
 
            plotOutput("qqPlot", width="70%"),

            br(),
            br(),

            strong('R session info'),
            verbatimTextOutput("info.out")

            ),



        tabPanel("About",

            strong('Note'),
            p('This web application is developed with',
            a("Shiny.", href="http://www.rstudio.com/shiny/", target="_blank"),
            ''),

            br(),

            strong('List of Packages Used'), br(),
            code('library(shiny)'),br(),
            code('library(shinyAce)'),br(),
            code('library(psych)'),br(),
            code('library(ltm)'),br(),
            code('library(CTT)'),br(),
            code('library(beeswarm)'),br(),
            br(),

            strong('Code'),
            p('Source code for this application is based on',
            a('"The handbook of Research in Foreign Language Learning and Teaching" (Takeuchi & Mizumoto, 2012).', href='http://mizumot.com/handbook/', target="_blank"),'I also referred to the code used in', a("MacR.", href="https://sites.google.com/site/casualmacr/", target="_blank")),

            p('The code for this web application is available at',
            a('GitHub.', href='https://github.com/mizumot/ctt', target="_blank")),

            p('If you want to run this code on your computer (in a local R session), run the code below:',br(),

            code('library(shiny)'),br(),
            code('runGitHub("ctt","mizumot")')
            ),

            br(),

            strong('Citation in Publications'),
            p('Mizumoto, A. (2015). Langtest (Version 1.0) [Web application]. Retrieved from http://langtest.jp'),

            br(),

            strong('Article'),
            p('Mizumoto, A., & Plonsky, L. (2015).', a("R as a lingua franca: Advantages of using R for quantitative research in applied linguistics.", href='http://applij.oxfordjournals.org/content/early/2015/06/24/applin.amv025.abstract', target="_blank"), em('Applied Linguistics,'), 'Advance online publication. doi:10.1093/applin/amv025'),

            br(),

            strong('Recommended'),
            p('To learn more about R, I suggest this excellent and free e-book (pdf),',
            a("A Guide to Doing Statistics in Second Language Research Using R,", href="http://cw.routledge.com/textbooks/9780805861853/guide-to-R.asp", target="_blank"),
            'written by Dr. Jenifer Larson-Hall.'),

            p('Also, if you are a cool Mac user and want to use R with GUI,',
            a("MacR", href="https://sites.google.com/site/casualmacr/", target="_blank"),
            'is defenitely the way to go!'),

            br(),

            strong('Author'),
            p(a("Atsushi MIZUMOTO,", href="http://mizumot.com", target="_blank"),' Ph.D.',br(),
            'Professor of Applied Linguistics',br(),
            'Faculty of Foreign Language Studies /',br(),
            'Graduate School of Foreign Language Education and Research,',br(),
            'Kansai University, Osaka, Japan'),

            br(),

            a(img(src="http://i.creativecommons.org/p/mark/1.0/80x15.png"), target="_blank", href="http://creativecommons.org/publicdomain/mark/1.0/"),

            p(br())

            )

))
))
