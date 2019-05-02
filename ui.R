#####################################
## VARIOUS PACKAGES REQUIRED
#####################################

  ## Front End 
    # Doc Viewer
      library(markdown)
    # Shiny Components
      library(shiny)
      library(shinyRGL)
      library(shinydashboard)
    # RGL components
      library(rgl)
      library(rglwidget)
    # Graphical Components
      library(ggplot2)
      library(plot3D)
      library(plotly)
    # Graphical Arrangement Components
      library(grid)
      library(gridBase)
      library(gridExtra)
      library(lattice)

#####################################
## USER INTERFACE CODE
#####################################

shinyUI(fluidPage(
  dashboardPage(skin = "blue",
                
  dashboardHeader(title = span(tagList("Minimum Density Hyperplane Clustering",icon=icon("fab fa-connectdevelop"))),
                  titleWidth=500),
  dashboardSidebar( 
    
  ## Sidebar with a slider input for number of bins
  sidebarPanel(

    tags$span(style="color:navy; font-family:Times; font-weight:bold; font-size=40px;",tags$strong("Options")),
               width=300, tags$style("body {background-color: #242d32;}"),
    
  selectInput('image',
               tags$span(style="color:navy; font-family: Times;",tags$strong('Image select')),
               choices = c('Dog', 'Robot','Tamarin', 'Dinner','Circle'),
               selectize=FALSE),
  radioButtons('DCS',
                tags$span(style="color:navy; font-family: Times;",tags$strong('Decorrelate and stretch image?')),  
                choiceNames = list(tags$span(style = "color:navy", "Yes"),
                tags$span(style = "color:navy", "No")),
                choiceValues = c(TRUE, FALSE), 
                selected=FALSE),  
  
  tags$style(".irs-bar, .irs-bar-edge, .irs-single{background: #222d32; border-color: #5a79ff;}
                .irs-grid-pol {display: none;} .irs-max {background: #222d32;}
                .irs-min {font-family: 'arial'; background: #222d32;}
                .irs-grid-text {display: none;}"),                             
  
  sliderInput("SPROP",
               tags$span(style="color:navy; font-family: Times;",tags$strong("Sample proportion:")),
               min=0.01, max=0.25, step=0.01, value=0.1),
  sliderInput("GAMMA",
               tags$span(style="color:navy; font-family: Times;",tags$strong("Reassignment \u0393 region:")),
               min=0.01, max=0.5, step=0.01, value=0.15),     
  sliderInput("ALPHA",
               tags$span(style="color:navy; font-family: Times;",tags$strong("Enforced solution plane distance from mean (\u03b1):")),
               min=0.1, max=3, step=0.1, value=0.5),

  radioButtons('ASS',
                tags$span(style="color:navy; font-family: Times;",tags$strong('Reassignment method:')),
                choiceNames = list(tags$span(style = "color:navy", "Right of plane"),
                                   tags$span(style = "color:navy", "Left of plane"),
                                   tags$span(style = "color:navy", "Heuristic shift")),
                choiceValues = c(1, 2, 3),
                selected=1),
  
tags$div(
tags$p(tags$span(style="color:navy; font-family: Times;float:left; padding:12px",
tags$strong("Submit changes:"))),
  submitButton(icon("paper-plane"),"Update Analysis", width='100%'),  

tags$p(tags$span(style="color:navy;font-family: Times; float:left; padding:12px",
"Results will vary when adjusting sampling proportions.")),

tags$p(tags$span(style="color:navy;font-family: Times; float:left; padding:12px",
"Information on minimum density hyperplane clustering can be found under the 'Documentation' tab.")),

tags$p(tags$span(style="color:navy; font-family: Times; padding:12px","Enjoy!"))))),
                
                
dashboardBody(   
  
  # Show a plot of the generated distribution
  tags$head(tags$style(HTML('.main-header .logo{
             font-family: "Times New Roman";
             font-weight: bold;
             font-size: 24px;}'))),
                  
  shinyUI(navbarPage("",
    tabPanel(title="Results",icon=icon("fas fa-camera-retro"),
      plotOutput("ImgPlot", width='100%', height = 500),

tags$div(tags$p(tags$span(style="color:black; padding:5px; float:left; font-family: courier new",
tags$strong("Top left:"),br(),
"Original Image.", br(),br(),
tags$strong("Top right image:"),br(),
"Minimum Density Solution, colour is representative of average RGB values within cluster.", br(), br(),
tags$strong("Bottom left image:"),br(),
"\u0393 region considered for reassignment coloured in cyan", br(),br(),
tags$strong("Bottom right image:"),br(),
"Final enhanced version of MDH solution.", br(),br(),br(),
tags$strong("NOTE: "),'Pixels are reassigned to the cluster class selected within "Reassignment method" options.',
"Reassignment method using the heuristic approach entails calculating first order kernel derivatives",
"for all \u03B3 points relative to the sample space.  This process may take up to an additional 15 seconds",
"for larger samples.")))),         
  
                                                                        
                                     
navbarMenu("Manual Separation",icon=icon("fas fa-filter"),
  tabPanel(title="Density Adjustment",icon=icon("fas fa-adjust"),
  br(), br(), br(),
    plotOutput("ManPlotDensity", width='100%'),br(),br(),br(),
    tags$style(".irs-bar, .irs-bar-edge, .irs-single{background: #222d32; border-color: #5a79ff;}
                .irs-grid-pol {display: none;} .irs-max {background: white;}
                .irs-min {font-family: 'arial'; background: white;}
                .irs-grid-text {display: none;}"),
      sliderInput("DD", 
                  tags$span(style="color:black; padding:5px; font-family: courier new; float:left; ","Manual Separation:"),
                  min = -3, max = 3, step=0.01, value=0, width='100%'),

          
tags$div(tags$p(tags$span(style="color:black; padding:5px; font-family: courier new; float:left; ",
'Choose a point along the density on which to separate the image into two ideal',
            'clusters and submit changes by clicking "Update Analysis".')),

tags$p(tags$span(style="color:black; padding:5px; float:left; font-family: courier new",
'Increasing x-axis values results in a larger cluster represented by red pixels in the image',
            'while decreasing the cluster size represented in black and vice-versa.')),


tags$p(tags$span(style="color:black; padding:5px; float:right; font-family: courier new",
"Density clustering seeks to find the minimum density",
"along which to segment an image into clusters which preserves as much information from the original image.",
"The challenge here is to maintain as much detail in the image with only two colours.",
"It may be necessary to increase sample proportion to attain a density estimation more similar to the original image,",
"thus possibly indicating a low density plane to optimally separate upon.",
"Some instances require a solution plane which may lay far from the mean and thus the user must indicate a larger \u03b1 value.",
"But care must be taken when adjusting \u03b1, large enough and the solution plane will minimise a density where few to no values are located,",
"resulting in a single clustered image."))))),

  navbarMenu("Image Densities",icon=icon("fab fa-stack-overflow"),
    tabPanel("2-Dimensional",icon=icon("fas fa-adjust"),
      plotOutput("DenPlot", width='100%', height = 700),
      
      
      
tags$div(tags$p(tags$span(style="color:black; font-family: courier new; padding:5px; float:left; ",
tags$strong("Top Left:"),br(),
"The density of each colour channel within the image.",
"The density of a colour at a certain intensity indicates how frequent that intensity is observed relative",
"to all other intensities of said colour.",br(),br(),
tags$strong("Top Right/Bottom:"),br(),
"The minimum density hyperplane (MDH) univarite projection solution with \u0393 region indicated in cyan(top).",
"The dashed line represents the solution plane while colours on either side are average channel colours within that cluster.",
"The penalised estimated density with final projected position of pixels (bottom).",
"The red line represents the MDH solution plane",
"Sampling the data will result in a different estimated kernel density and therefore a different solution.",
"Interestingly, some instances require relatively few observations for an acceptable MDH result.",
"Adjusting sample proportion down to 1% of the image in some cases segments the image better than clustering on a full image data set.",
br(),br(),
"Note: More insight can be gleamed by viewing the solution plane and pixel dispersion in their natural domain",
"within the 3-dimensional plot tabs.")))),

  
  tabPanel("3-Dimensional",icon=icon("fas fa-adjust"),
      plotlyOutput("plot"),
      

tags$div(tags$p(tags$span(style="color:red; font-family: courier new; float:left; padding:12px",
tags$strong("Please be patient while the interactive 3-dimensional scatter plot renders.")))),
tags$div(tags$p(tags$span(style="color:black; font-family: courier new; float:left; padding:12px",
tags$strong("Points represent red, green, and blue intensities drawn from a sample within the image.",
"Each point are coloured as they are found within the image for clarity of association.")))))),

navbarMenu("Adjusted Densities",icon=icon("fab fa-stack-overflow"),
    tabPanel("3-Dimensional \u0393 region",icon=icon("fas fa-adjust"),
      plotlyOutput("GAMplot"),
      
tags$div(tags$p(tags$span(style="color:red; font-family: courier new; float:left; padding:12px",
tags$strong("Please be patient while the interactive 3-dimensional scatter plot renders.")))),
tags$div(tags$p(tags$span(style="color:black; padding:5px; float:left; font-family: courier new",
"Each sampled pixel is coloured according to the mean colour of cluster with \u0393 region highlighted in cyan")))),

  
  tabPanel("3-Dimensional final assignment",icon=icon("fas fa-adjust"),
      plotlyOutput("Finalplot"),
      

tags$div(tags$p(tags$span(style="color:red; font-family: courier new; float:left; padding:12px",
tags$strong("Please be patient while the interactive 3-dimensional scatter plot renders.")))),
tags$div(tags$p(tags$span(style="color:black; font-family: courier new; float:left; padding:12px",
tags$strong("Sampled points in the plot represent the final adjusted segmented image.",
"All \u0393 region points are reassigned via 'Reassign cluster to?' option.",
"The mean shift heuristic enhancement will be introduced as a method of reassignment in future versions.")))),


tags$div(tags$p(tags$span(style="color:black; font-family: courier new; float:left; padding:12px",
tags$strong("\u0393 Region:"),br(),
"Increase or decrease the percentage of values around the MDH solution plane to reassign to pre-defined cluster and",
"then head over to results to see the final adjusted image."))))),

  navbarMenu("Documentation",icon=icon("far fa-file"),                              
  tabPanel("R Program", icon=icon("fas fa-adjust"),
    htmlOutput("RESMDH")),
   tabPanel("Thesis Abstract", icon=icon("fas fa-adjust"),
    mainPanel(tags$head(tags$style(HTML(".main{overflow-y: auto;}"))),
      includeHTML("./www/MDH Abstract.html"))))

   ))
   )
   )
  ))
