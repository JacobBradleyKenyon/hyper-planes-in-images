#####################################
## VARIOUS PACKAGES REQUIRED
#####################################

  ## Back End
    # Project Pursuit Package for Hyperplane SOlution
      library(PPCI)
    # JPEG package to read in Image
      library(jpeg)
    # Package to convert Data into Long Format
      library(reshape2)
    # Various other packages for mathematical calculations etc.
      library(akima) 
      library(grpss)
      library(ks)
      library(Matrix)


#####################################
## VARIOUS FUNCTIONS 
#####################################

  ## Function to read in jpeg with indices
    Read.jpeg <- function(FILE, Image_Name, Url=TRUE){
  ## Function to download jpeg images from the internet using full html address
  ## Function requires jpeg package


    ## Empty list to collect URLs of image/images
      Images<-list()
    if(Url==TRUE){
    ## Download image/images from url FILE
      for(i in 1:length(FILE)){
              z <- tempfile()
              download.file(FILE[i], z, mode = 'wb')
              Images[[i]] <- jpeg::readJPEG(z)
              file.remove(z)
      }
    }else{
      for(i in 1:length(FILE)){
      Images[[i]] <- jpeg::readJPEG(FILE[i])
      }
    }
    ## Indexing images with x and y coordinates
      ## Empty set to index images
        Indexed <- list()

      ## indexing image/images
        for(i in 1:length(Images)){
            Indexed[[i]] <- data.frame( x = rep(1:dim(Images[[i]])[2], each = dim(Images[[i]])[1]),
                                        y = rep(dim(Images[[i]])[1]:1, dim(Images[[i]])[2]),
                                        R = as.vector(Images[[i]][,,1]),
                                        G = as.vector(Images[[i]][,,2]),
                                        B = as.vector(Images[[i]][,,3]))
        }

  ## Output final indexed images with coordinates and R, G, B values measuring intensities of each value [0,1]
    names(Indexed) <- Image_Name
    Indexed
  }

  ## Contrasting image function, a.k.a. DeCorrelation Stretch
    DC.Stretch<-function(input, Up = 1, Low = 0){
  ## A function to decorrelate and stretch range of input according to upper limit "Up" and lower limit "Low"
    
    ## Decorrelating input
      DC <- as.matrix(input) %*% eigen(cov(input))$vector
    
    ## Stretching input according to Up and Low
      Stretch <- function(input, Up=Up, Low=Low){
        Max<-max(input)
        Min<-min(input)
        ((Up-Low)/(Max-Min))*(input-Min)+Low
      }
     
    ## Applying Stretching function to decorrelated input 
      DCS <- apply(DC, 2, function(x) Stretch(x, Up=Up, Low=Low))
  
  ## final decorrelated stretched output from input  
  DCS
  }
  
  ## Function to plot images
    Plot.img <- function(data = DCS.img, cluster = 0, Title = '', Gamma.region.color = 'black', preprocessed = FALSE, original.data = Original.img){
    ## Function to plot image data containing x y R G B 
    ## Option to plot image according to  mean RGB values of each cluster
    ## If region around mdh solution is included within cluster assignment then 
    ## Gamma.region.color will be colored according to selection

    ## Choosing an CEX size according to image size 
      NR <- nrow(data)
      CEX <-ifelse(NR <= 2500, 4.51, 
                 ifelse(NR <= 10000, 2.1,
                        ifelse( NR <= 40000, 0.85,0.225)))
    
    ## Setting theme for pictures
    THEME.Pic<-theme_bw() + theme(text=element_text(family='sans'),
                                  panel.border     = element_blank(), 
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  plot.title       = element_text(size=12, hjust=0.5),
                                  axis.line        = element_blank(),
                                  axis.text        = element_blank(),
                                  axis.ticks       = element_blank(),
                                  axis.title       = element_blank(),
                                  plot.margin      = rep(unit(0,"cm"),4),
                                  panel.spacing    = unit(0,"null"),
                                  axis.ticks.length = unit(0,"null"))
                                  
    
    ## Generating SuperPixel color based on cluster
    if(cluster[1]!=0){
      
      ## To maintain original mean color clusters while using preprocessed images,
      ## indicate preprocessed==TRUE & include original image data (original.data)
      if(preprocessed==TRUE){
        Data <- cbind(original.data, cluster)
      }
      else{
        Data <- cbind(data, cluster)
      }
      
      RGB.mean <- Data %>% group_by(cluster) %>% dplyr::summarize(R = mean(R), G = mean(G), B =mean(B))
      
      ## Mean R G B values to use within plot
      if(nrow(RGB.mean) == 1){
        levels(RGB.mean$cluster)<-c(1,2,3)
        COLOR <- rbind(RGB.mean, RGB.mean, RGB.mean)
        COLOR[3,1]<-3
        COLOR[2,1]<-2
      }
      else
       COLOR <- data.frame(RGB.mean)

      ## Plotting image according to cluster mean color
      plt<-ggplot(data=data, aes(x=x, y=y, col=ifelse(cluster == 1, 'Cluster 1',
                                                 ifelse(cluster == 2, 'Cluster 2', 'Gamma')))) +
        geom_point(size=as.numeric(CEX), shape=15) + 
        scale_color_manual(values=c(rgb(COLOR[1,2:4]), rgb(COLOR[2,2:4]), Gamma.region.color)) +
        coord_fixed() +
        ggtitle(Title) +
        THEME.Pic +
        theme(legend.position = '')   
    }
    
    else{
      ## Original R G B values to use within plot  
      COLOR <- data[, c('R', 'G', 'B')]  
      
      ## Plotting image according to original R G B values
      plt<-ggplot(data=data, aes(x=x, y=y)) +
        geom_point(size=as.numeric(CEX), shape=15, col=rgb(COLOR[,c('R', 'G', 'B')])) + 
        scale_color_manual() +
        coord_fixed() +
        ggtitle(Title) +
        THEME.Pic 
    }
    list(Colors=COLOR, Plot=plt)
  }  

  ## Function to pull values within gamma of the hyperplane
    GammaHype <- function(gamma = 0.1, B = sol.mdh[[1]]$b, V = sol.mdh[[1]]$v, CL = sol.mdh[[1]]$cluster, input = DCS.img, Picture = FALSE){
    ## Function that takes an input and outputs a subset of pixels related to a gamma region (i.e. half-width interval around mdh hyperplane solution)
    ## The gamma parameter captures an interval that contains a % of the data around B
    
    if(!Picture){
    ## Checks to ensure user has correct format
      if(length(V) != ncol(input))
        stop('Error in input %*% V : non-conformable arguments.\nInput must have columns equal to rows of projection vector V.\n')
    
      ## Linear transformation according to MDH solution vector V
        Xv <- as.matrix(input) %*% as.matrix(V)
    }
    
    if(Picture){
    ## Checks to ensure user has correct format
      CN <- colnames(input)
      FN <- sum(CN=='R' | CN=='B' | CN=='G')
      if(FN!=3)
        stop('Incorrect data input format.\nPlease rename red, green, and blue column names to "R", "G" , and "B".')
      
      ## Linear transformation according to MDH solution vector V
      Xv <- as.matrix(input[,c('R','G','B')]) %*% as.matrix(V)
    }
     
    ## Creating half width length around MDH hyperplane solution (B)
      hw <- diff(range(Xv))*(gamma/2)
    
    ## Collecting index of observations around hyperplane 
      Index <- which(Xv >= (B-hw) & Xv <= (B+hw))
    
    ## Combining current mdh cluster assignment with data
      data <- cbind(input, CL)
        colnames(data)[(ncol(input)+1)] <- 'cluster'
    
    ## Subsetting input via index created according to gamma half width to mdh hyperplane solution
      Output <- data[Index,]  
    
  ## Output half-width, mdh projected solution, index of gamma values, and subset of points within gamma of hyperplane
    list(hw=hw, sol.Xv=Xv, Index=Index, Points=Output)
  }

    
#####################################
## SERVER SIDE CODE
#####################################

shinyServer(function(input, output) {

  
  ## URLs of images 
  URLs<-c("https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcT2w5H7Jpr0_DJkpLpCZE80T7aiRzeJtyjhH6wQpPwwUfg8XRbW-w",
          "http://i42.photobucket.com/albums/e315/tntoxfox/wall-esmall.jpg",
          "https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcQD7G1JmkGeFg63G99upZycjQhq_9VZN9V25Vqx3tK9Loe2MNQgkQ",
          "http://static.newworldencyclopedia.org/thumb/6/62/CreoleFood.jpg/200px-CreoleFood.jpg",
          "http://wiki.metin2.co.uk/images/thumb/0/06/Ds_attack_round.jpg/200px-Ds_attack_round.jpg")
  
  ## Importing images using function from line 5
  img<-Read.jpeg(Url=FALSE,
                 FILE = c("./www/dog.jpg",
                              "./www/robot.jpg",
                              "./www/tamarin.jpg",
                              "./www/dinner.jpg",
                              "./www/circle.jpg"),
                 Image_Name = c('Dog', 'Robot', 'Tamarin', 'Dinner','Circle'))
 

  ## Theme for plots
  THEME<-theme_bw() + theme(text=element_text(family='serif')) +
    theme(panel.border=element_blank(), panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), axis.line=element_line(colour='gray'),
          plot.title=element_text(face='bold', size=15, hjust=0.5))
  
  
  datasetInput <- reactive({
    switch(input$image,
           "Dog" = "Dog",
           'Robot' = 'Robot',
           'Tamarin'='Tamarin',
           'Dinner'='Dinner',
           'Circle'='Circle') 
  })
  
  ContrastInput<-reactive({
    DTA      <- datasetInput()
    pic      <- img[[DTA]]
    Contrast <- cbind(pic[,1:2], DC.Stretch(pic[,3:5]))
    
    colnames(Contrast) <- c('x', 'y', 'R', 'G', 'B')
    
    list(Contrast=Contrast, pic=pic)
  })
  
  DCSInput <- reactive({
     input$DCS
  }) 
  
  S.Prop <- reactive({
    input$SPROP
  })
  
  Indy<-reactive({
    Contrast <- ContrastInput()
    S.prop   <- S.Prop()
    Indy     <- sample(x=1:nrow(Contrast$Contrast), size=ceiling(nrow(Contrast$Contrast)*S.prop), replace = FALSE)
    
    return(list(Indy=Indy))
    })
  
  MethodInput<-reactive({
    Contrast <- ContrastInput()
    Alpha    <- alphainput()
    Reassign <- AssInput()
    Gamma    <- gammainput()
    Indy     <- Indy()$Indy
    BOUND    <- 0 
    DCS      <- DCSInput()
        
      if(DCS==TRUE)
        Contrast<-Contrast$Contrast
      else
        Contrast<-Contrast$pic
    
    
       USE  <- as.matrix(Contrast[Indy,c('R','G','B')])
       sol.mdh <- mdh(X = USE, alphamin = Alpha, alphamax = Alpha)
       

      V       <- sol.mdh[[1]]$v
      B       <- sol.mdh[[1]]$b
      BW      <- sol.mdh[[1]]$params$h
      oim     <- as.matrix(Contrast[,3:5]) %*% V 
      Clust   <- ifelse(oim > B, 1, 2)
        rownames(Clust) <- seq(1:nrow(Clust))

    O.lab <- Clust
    
    ## Gamma Region 
    Gamma.region <- GammaHype(gamma=Gamma, input = Contrast, B = B, V = V, CL = Clust, Picture = TRUE)
    
    ## Gamma Region RGB values
    Gam.points <- Gamma.region$Points
    G.index    <- Gamma.region$Index
    
    if(Reassign==3){
    Heuristic <- kdde(x=Contrast[Indy,3:5], H=diag(rep(BW,3)), eval.points=Gam.points[,3:5], deriv.order = 1)
      Uni.side<- Heuristic$estimate %*% V 
      NewLabel  <- ifelse(Uni.side > B, 1, 2)
      O.lab[Gamma.region$Index]<-NewLabel
    }else
      O.lab[Gamma.region$Index]<-Reassign
    
      return(list(V=V, B=B, Clust=Clust, BW=BW, sol=sol.mdh, Gamma.region=Gamma.region, G.index=G.index, Gamma.Labels=O.lab))
  })
  
  Meltinput<-reactive({
    Contrast <- ContrastInput()
    DCS      <- DCSInput()
    
      if(DCS==TRUE)
        Contrast <- Contrast$Contrast
      else
        Contrast <- Contrast$pic
    
    melt(Contrast[,3:5])

  })
  
  getPage <- function(X){
    return(includeHTML(X))}
    
  DensityInput<-reactive({
    Contrast <- ContrastInput()
    DCS      <- DCSInput()
    
      if(DCS==TRUE)
        Contrast<-Contrast$Contrast
      else
        Contrast<-Contrast$pic
    
    Clust  <- MethodInput()
    B      <- Clust$B
    V      <- Clust$V
    Rezone <- as.matrix(Contrast[,c('R','G','B')]) %*% V
    dens   <- density(Rezone, n=nrow(Contrast)) 
    DEN    <- as.data.frame(cbind(x=dens$x, y=dens$y))
    return(list(Rezone=Rezone, DEN=DEN, B=B))
  })
  
  cexInput<-reactive({
    Contrast <- ContrastInput()
    DCS      <- DCSInput()
    
      if(DCS==TRUE)
        Contrast <- Contrast$Contrast
      else
        Contrast <- Contrast$pic
    
    NR  <- nrow(Contrast)
    CEX <- ifelse(NR <= 2500, 4.51, 
                 ifelse(NR <= 10000, 2.1,
                        ifelse( NR <= 40000, 0.85,0.225)))
    return(CEX)
  })
  
  DimInput<-reactive({
    Contrast <- ContrastInput()
    DCS      <- DCSInput()
    
      if(DCS==TRUE)
        Contrast <- Contrast$Contrast
      else
        Contrast <- Contrast$pic
    
    fhat <- kde(Contrast[,3:5])
    return(fhat)
  })
  
  inputDD<-reactive({
    input$DD})
  
  gammainput<-reactive({
  input$GAMMA
    })
  
  alphainput<-reactive({
    input$ALPHA
  })
  
  AssInput <- reactive({input$ASS})
  
## Output image according to ImgPlot
  output$ImgPlot <- renderPlot({

  ## Collecting Reactives
  Contrast  <- ContrastInput()
  O.pic     <- Contrast$pic
  Contrast  <- Contrast$Contrast

  Outsidein <- MethodInput()
  G.labs    <- Outsidein$Gamma.Labels
  G.points  <- Outsidein$Gamma.region
  CL        <- Outsidein$Clust
  
  DCS       <- DCSInput()

  CEX       <- cexInput()

  ## Updating Cluster List (CL) to account for 3rd cluster, Gamma cluster
  CL.gamma <- CL
  CL.gamma[G.points$Index] <- 3

  ## Illustrate DCS image or not
  if(DCS==TRUE){

  p1 <- Plot.img(data = Contrast, Title = 'Decorrelated & Stretched Image')

  p2 <- Plot.img(data = Contrast, cluster = CL, Title = 'MDH Solution',
                 preprocessed = FALSE, original.data = O.pic)

  p3 <- Plot.img(data = Contrast, preprocessed = FALSE, Title = '\u0393 Region',
                 cluster = CL.gamma, Gamma.region.color = 'cyan')

  p4 <- Plot.img(data = Contrast, cluster = G.labs, Title = 'Adjusted Solution',
                 preprocessed = FALSE, original.data = O.pic)
  }else{

  p1 <- Plot.img(data = O.pic, Title = 'Original Image')

  p2 <- Plot.img(data = Contrast, cluster = CL, Title = 'MDH Solution',
                 preprocessed = TRUE, original.data = O.pic)

  p3 <- Plot.img(data = O.pic, cluster = CL.gamma, Title = '\u0393 Region', Gamma.region.color = 'cyan')

  p4 <- Plot.img(data = Contrast, cluster = G.labs, Title = 'Adjusted Solution',
                 preprocessed = TRUE, original.data = O.pic)
  }




grid.arrange(p1$Plot, p2$Plot, p3$Plot, p4$Plot, ncol=2, nrow=2)

})

  output$DenPlot <- renderPlot({
    ## Collecting Reactives of image and cluster by method
    Contrast <- ContrastInput()
    O.pic    <- Contrast$pic
    Contrast <- Contrast$Contrast 
    G.dex    <- gammainput()
    Clust    <- MethodInput()
    Sol      <- Clust$sol
    Indy     <- Indy() 
    USE      <- Indy$Indy      
    B        <- Clust$B
    V        <- Clust$V
    Clust    <-Clust$Clust
    
    ## Using reactive to melt RGB for density plots
    RGB <- Meltinput()
    
    ## Using reactive to select Rezoned values and density   
    DEN<-DensityInput()
    Rezone <- DEN$Rezone
    DEN<-DEN$DEN
    
    dens <- density(Rezone, n=nrow(Contrast)) 
    DEN <- as.data.frame(cbind(x=dens$x, y=dens$y))
    
    
    p2<-ggplot(data=RGB, aes(x=value, fill=variable, color=variable)) +
          geom_density(alpha=0.35) +
            labs(x="Color Intensity", y = "Density", title="Initial RGB Density") +
              THEME +
                theme(legend.position = c(0.99,0.8),
                      legend.title = element_text(colour = 'white'),
                      legend.spacing.y = unit(4,'cm'))
    
    
  ## Illustrate DCS colours or original in density plot
  DCS<-DCSInput()
 
  if(DCS==TRUE){
    Colors <- Plot.img(data = Contrast, cluster = Clust, Title = 'Clustered Image',
                 preprocessed = FALSE, original.data = O.pic)$Colors
    DATA<-as.matrix(Contrast[,3:5])
    }
    
  if(DCS==FALSE){
    Colors <- Plot.img(data = Contrast, cluster = Clust, Title = 'Clustered Image',
                 preprocessed = TRUE,
                 original.data = O.pic)$Colors
    DATA<-as.matrix(O.pic[,3:5])
  }
   
  ## Gamma region colour
  G.low<-B-(G.dex/2)
  G.hig<-B+(G.dex/2)
  ecdf_fun <- function(x,perc) ecdf(x)(perc)
  Lower<-ecdf_fun(DEN$x, G.low)
  Upper<-ecdf_fun(DEN$x, G.hig)

  
  p4<-ggplot(data=DEN, aes(x, y)) +
    geom_ribbon(data=subset(DEN, x <= quantile(DEN$x,Lower)[1]),
                aes(x=x, ymax=y), ymin=0, fill=rgb(Colors[2,2:4]), alpha=0.4)+
    geom_ribbon(data=subset(DEN, x > quantile(DEN$x,Lower)[1] & x <= quantile(DEN$x,Upper)[1]),
                aes(x=x, ymax=y), ymin=0, fill=rgb(0, 1, 1), alpha=0.4)+
    geom_ribbon(data=subset(DEN, x > quantile(DEN$x,Upper)[1]),
                aes(x=x, ymax=y), ymin=0, fill=rgb(Colors[1,2:4]), alpha=0.4) +
    geom_vline(xintercept = B, linetype='dashed') +
    ggtitle("Minimum Density Hyperplane") +
    labs(x="Optimized values for separation", y = "") +
    THEME +
    scale_x_continuous(breaks = round(seq(min(DEN$x), max(DEN$x), by = 0.2),1))
    
  plot.new()
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))

  #Draw ggplot
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p4, newpage = FALSE)
  popViewport()
  
  #Draw ggplot
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, newpage = FALSE)
  popViewport()
  
  #Draw base plot
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = c(1:2)))
  par(fig = gridFIG(), new = TRUE)
  hp_plot(Sol, DATA)
  popViewport()
  
#grid.arrange(p2, p4, ncol=2)
})
  
  output$ManPlotDensity<- renderPlot({
    ## Manual Adjustments
    Contrast <- ContrastInput()
    O.pic    <- Contrast$pic
    Contrast <- as.data.frame(Contrast$Contrast)
    
    DEN      <- DensityInput()
    Rezone   <- DEN$Rezone
    DEN      <- DEN$DEN
    MD       <- input$DD
    CEX      <- cexInput()

    Contrast$MCG <- ifelse(Rezone <= MD, 'red', 'black')

    p5 <- ggplot(data=Contrast, aes(x=x, y=y, col=MCG)) + 
            geom_point(shape=15)+ 
              scale_color_identity() +
                coord_fixed() +
                  ggtitle("Manual Clustered Image") +
                    THEME +
                      theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    
    p6 <- ggplot(data=DEN, aes(x, y)) +
            geom_ribbon(data=subset(DEN, x <= MD), aes(x=x, ymax=y), ymin=0, fill='red', alpha=0.4)+
            geom_ribbon(data=subset(DEN, x > MD), aes(x=x, ymax=y), ymin=0, fill='black', alpha=0.4) +
              geom_vline(xintercept = MD, linetype='dashed') +
                labs(x="Solution values", y = "Density", title="Kernel Density of RGB Values") +
                  THEME +
                    scale_x_continuous(breaks = round(seq(min(DEN$x), max(DEN$x), by = 0.2),1))

    grid.arrange(p5, p6, nrow = 1)
  })
  
  output$ManPlotImage<- renderPlot({
    ## Manual Adjustments
    Contrast<-ContrastInput()
    O.pic <- Contrast$pic
    Contrast<-Contrast$Contrast
    
    DEN<-DensityInput()
    Rezone<-DEN$Rezone
    DEN<- DEN$DEN
    MD<-input$DD 
    
    CEX<-cexInput()
    
    Contrast$MCG <- ifelse(Rezone <= MD, 'red', 'black')
    
    p5 <-ggplot(data=Contrast, aes(x=x, y=y, col=MCG)) + 
          geom_point(size=as.numeric(CEX), shape=15)+ 
            scale_color_identity() +
              coord_fixed() +
                ggtitle("Manual Clustered Image") +
                  THEME +
                    theme(axis.line=element_blank(),
                          axis.text.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
    
    print(p5)
  })
  
  output$plot<-renderPlotly({
    
    ## Reactive elements
    Contrast <- ContrastInput()
    USE      <- Indy()
    USE      <- USE$Indy
    DCS      <- DCSInput()
 
      if(DCS==TRUE)
        Contrast <- Contrast$Contrast[USE,]
      else
        Contrast <- Contrast$pic[USE,] 
    
    ax <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)    
    
    plot_ly(Contrast, x = ~R, y = ~G, z = ~B,
            hoverinfo="none", 
            marker = list(size = 3,
                          color = rgb(Contrast[,3:5]))) %>%
            add_markers()%>%
            layout(xaxis = ax, yaxis=ax, showlegend = FALSE, plot_bgcolor = '#222d32',
                   paper_bgcolor = '#222d32') 

    })
  
  output$GAMplot<-renderPlotly({
    GWIDTH   <- gammainput()
    Contrast <- ContrastInput()
    DCS      <- DCSInput()
    
    Outsidein<- MethodInput()
    G.index  <- Outsidein$G.index
    B        <- Outsidein$B
    V        <- Outsidein$V
    BW       <- Outsidein$BW
    CL       <- Outsidein$Clust
    Indy     <- Indy()
    USE      <- Indy$Indy
    
    if(DCS==TRUE)
      Contrast <- Contrast$Contrast[USE,]
    else
      Contrast <- Contrast$pic[USE,]
    
    
    ## Updating Cluster List (CL) to account for 3rd cluster, Gamma cluster
      CL.gamma <- CL
     # CL.gamma[Gamma.region$Index] <- 3
      CL.gamma[G.index] <- 3
    
    p3 <- Plot.img(data = Contrast, preprocessed = FALSE, Title = '\u0393 Region',
                   cluster = CL[USE])$Color


    ax <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)    

    COLORS<-c(rgb(p3[1,2:4], maxColorValue = 1),
              rgb(p3[2,2:4], maxColorValue = 1), rgb(red=0, blue = 1, green = 1, maxColorValue = 1))
    
    Contrast<-cbind(Contrast, CL.gamma[USE])

    Contrast<-as.data.frame(Contrast)
    Contrast$CL.gamma <- ifelse(Contrast$CL.gamma == 1, COLORS[1], 
                         ifelse(Contrast$CL.gamma == 2, COLORS[2], COLORS[3]))

    plot_ly(Contrast, x = ~R, y = ~G, z = ~B, type="scatter3d", hoverinfo="none",
            marker = list(size=3,   
                          color = ~CL.gamma,
                          colors= c(COLORS[1], COLORS[2], COLORS[3]))) %>%  
      
      layout(xaxis = ax, yaxis=ax, showlegend = FALSE, plot_bgcolor = '#222d32',
             paper_bgcolor = '#222d32') 
  })
 
  output$Finalplot<-renderPlotly({
    Gamma.lab <- MethodInput()
    Indy      <- Indy()
    USE       <- Indy$Indy
    CL        <- Gamma.lab$Clust
    G.index   <- Gamma.lab$G.index
    Gamma.lab <- Gamma.lab$Gamma.Labels

    Contrast <- ContrastInput()
    
    DCS      <- DCSInput()

    if(DCS==TRUE) 
      Contrast <- Contrast$Contrast[USE,]
    else
      Contrast <- Contrast$pic[USE,]
    

    p3 <- Plot.img(data = Contrast, preprocessed = FALSE, Title = '\u0393 Region',
                   cluster = Gamma.lab[USE])$Color

    ax <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)    

    COLORS<-c(rgb(p3[1,2:4], maxColorValue = 1), rgb(p3[2,2:4], maxColorValue = 1), rgb(red=0, blue = 1, green = 1, maxColorValue = 1))
    
    Contrast<-cbind(Contrast, Gamma.lab[USE])

    Contrast<-as.data.frame(Contrast)
    Contrast$CL.gamma <- ifelse(Contrast$Gamma.lab == 1, COLORS[1], 
                         ifelse(Contrast$Gamma.lab == 2, COLORS[2], COLORS[3]))

      plot_ly(Contrast, x = ~R, y = ~G, z = ~B, type="scatter3d", hoverinfo="none",
        marker = list(size=3, color=~CL.gamma, colors=c(COLORS[1], COLORS[2]))) %>%  
           layout(xaxis=ax, yaxis=ax, showlegend=FALSE, plot_bgcolor='#222d32', paper_bgcolor='#222d32') 
    })
  
output$RESMDH   <- renderUI({getPage("https://cran.r-project.org/web/packages/PPCI/index.html")})

})


