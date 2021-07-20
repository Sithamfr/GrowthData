# SHINY
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)

# IMAGES
library(ijtiff)
library(raster)
library(imager)

# VISUALISATION
library(ggplot2)
library(hrbrthemes)
library(plotly)
library(scales)

# DATA
library(data.table)


options(shiny.maxRequestSize=2000*1024^2) # Import file max size : 2GB

dic_col = list("Red"=1,"Green"=2,"Blue"=3)
dic_col_inv = c("red","green","blue")

RVB = function(v){
  return(rgb(v[1],v[2],v[3]))
}

signe = function(x){
  s = sign(x)
  if (s==1){return("+")}
  return("-")
}

extract_var = function(mat3d){
  v = apply(mat3d,3,as.numeric)
  return(sum(apply(v,2,var)))
}

notna = function(value){
  if (is.na(value)) {return(0)}
  return(value)
}

decompimgNG = function(mat,col_couche){
  df = data.frame(X=as.numeric(mat))
  colnames(df) = "X"
  p = ggplot(df,aes(x=X,fill=..x..)) + geom_histogram(bins=50) + theme_classic() +
    scale_fill_gradient(low="black",high=col_couche) + 
    theme(legend.position="none") + labs(x="Pixel intensities", y="Occurences")
  return(p)
}

decompimgRGB = function(im,threshold){
  df = as.data.frame(cbind(as.numeric(im[,,1]),as.numeric(im[,,2]),as.numeric(im[,,3])))
  colnames(df) = c("R","G","B")
  df = setDT(df)[, .N, by = c("R","G","B")]
  df = df[df$N>threshold,]
  col = apply(df,1,RVB)
  df = as.data.frame(cbind(df,col))
  p = plot_ly(df, x=~R, y=~G, z=~B,
          text=~paste("Color : ",col,"<br> R:",R,",G:",G,",B:",B,"<br> Number of occurences : ",N),
          hoverinfo = 'text',
          marker = list(size = ~2*log(N+1),
                        color = ~col,
                        line = list(width = 0))) %>%
      add_markers() %>%
      layout(title="3D Histogram",
           scene = list(xaxis = list(title = 'Red'),
                        yaxis = list(title = 'Green'),
                        zaxis = list(title = 'Blue'))
    )
  return(p)
}

segmentation = function(im,bin=TRUE){
  if (bin) {return(im>0)}
  else {
    return(im)
  }
}

harmonique = function(X,period,n){ #X un vecteur (serie chrono), n la n-ieme periode
  t = 2*pi*(n*period+1):((n+1)*period)/T
  a = 2/period*sum(X[n,]%*%cos(t))
  b = 2/period*sum(X[n,]%*%sin(t))
  c = sqrt(a^2+b^2)
  return(c)
}

energie = function(X,period,n){  #X un vecteur (serie chrono), n la n-ieme periode
  return(1/period*sum(X[n,]^2))
}

hdr = function(E,c){  #E energie d'un signal et c sa premiere harmonique
  return(100*sqrt(abs(2*E/c^2-1)))
}

pente = function(X){ #inutilise pour l'instant
  df = data.frame(X,Tps=1:length(X))
  trend = lm(X~Tps,data=df)
  return(as.numeric(trend$coefficients[2]))
}

evol3d_feats = function(signal,period,h,couleur="blue"){ # h est le décalage pour se situer en minuit (debut de periode)
  #X = 100*signal[(h+1):length(signal)]
  X = 2000*signal[(h+1):length(signal)]
  n = length(X)
  nT = n%/%period
  if (n%%period!=0){
    N = n%%period
    X = X[1:(length(X)-N)]
    n = length(X)
  }
  Tps = 1:n
  t = 1:nT
  Z = matrix(0,nrow=nT,ncol=period)
  pente = rep(0,nT)
  for (i in t){
    tmp = 1:period
    X_tmp = X[(period*(i-1)+1):(period*i)]
    trend = lm(X_tmp ~ tmp)
    pente[i] = trend$coefficients[2]
    trend = trend$coefficients[1] + trend$coefficients[2]*tmp
    Z[i,] = X_tmp - trend
  }
  c = rep(0,nT)
  E = rep(0,nT)
  HDR = rep(0,nT)
  for (i in t){ #Pour chaque periode
    c[i] = harmonique(Z,period,i)
    E[i] = energie(Z,period,i)
    HDR[i] = hdr(E[i],c[i])
    #pente[i] = pente(Z[i,])
  }
  M = as.data.frame(cbind(c,E,HDR,rep(couleur,nT),t,pente))
  colnames(M) = c("Harmonic","Energy","HDR","Couleur","N_period","Pente")
  plot_ly(M,
          x=~Pente,
          y=~Energy,
          z=~HDR,
          line=list(color=~Harmonic),
          text=~paste("Harmonic : ",Harmonic,
                      "<br>PERIOD ",N_period),
          type = 'scatter3d',
          mode = 'lines') %>%
  layout(scene=list(xaxis=list(title="Slope",zeroline=FALSE),
  			yaxis=list(zeroline=FALSE),
  			zaxis=list(zeroline=FALSE)))
}



# SHINY APP



ui <- dashboardPage(
  skin = "blue",
  title="GrowthData",
  dashboardHeader(tags$li(class = "dropdown",
                          tags$style(".main-header {max-height: 0px}"),
                          tags$style(".main-header .logo {height: 0px;}"),
                          tags$style(".sidebar-toggle {height: 0px; padding-top: 0px !important;}"),
                          tags$style(".navbar {min-height: 0px !important}")
                  ),
                  disable = TRUE
  ),
  dashboardSidebar(
    width = 200,
    tags$style(".main-sidebar {padding : 10px}"),
    tags$style(".small-box {opacity:0.7;border-radius:4px}"),
    sidebarMenu(
      id = "tabs",
      div(align="center",img(src='logo_inUI.png', width="40%")),
      br(),
      menuItem("HOME", tabName = "home", icon = icon("home"), badgeLabel="●", badgeColor="blue"),
      br(),br(),
      HTML("<h4 align='center'>GROWTH ANALYSIS<br>OF PLANTS</h4>"),
      br(),br(),br(),
      HTML("<div align='center'><font color='D5CBBB'>Color images</font></div>"),
      HTML("<div align='right'><font size='1.5'>
           <font color='#B65151'>R</font>
           <font color='#7FBC95'>G</font>
           <font color='#4B91C8'>B</font>
           : <i><font color='#B4A68F'>3D</font>
           </font></i></div>"),
      menuItem("Upload", tabName = "up", icon = icon("chevron-down")),
      menuItem(HTML("&nbsp;&nbsp;Image analysis"), tabName="im_ana", icon = icon("image")),
      menuItem("Sequence analysis", tabName="seq_ana", icon = icon("images")),
      br(),br(),
      HTML("<div align='center'><font color='D5CBBB'>Depth images</font></div>"),
      HTML("<div align='right'><font size='1.5'>
           <font color='#A0A0A0'>G</font>
           <font color='#AAAAAA'>r</font>
           <font color='#B4B4B4'>a</font>
           <font color='#BEBEBE'>y</font>
           <font color='#C8C8C8'>s</font>
           <font color='#D2D2D2'>c</font>
           <font color='#DCDCDC'>a</font>
           <font color='#E6E6E6'>l</font>
           <font color='#F0F0F0'>e</font>
           : <i><font color='#B4A68F'>1D</font>
           </font></i></div>"),
      menuItem("Upload", tabName = "depth_up", icon = icon("chevron-down")),
      menuItem("Sequence analysis", tabName="depth_ana", icon = icon("images"))
    )
  ),
  dashboardBody(
    tags$head(tags$link(rel="shortcut icon", href="logo.png")),
    setBackgroundImage(src="bg.jpg",shinydashboard=TRUE),
    tabItems(
      tabItem(tabName="home",
              column(1),
              column(5,
                     div(style="border:solid;border-color:white;
                         border-top-left-radius:8px;border-top-right-radius:8px;
                         padding-left:20px;padding-right:20px;
                         margin-right:20px;margin-top:20px",
                     div(align="center",
                         style="background-color:white; opacity:0.8;
                         border-bottom-left-radius:5px;
                         border-bottom-right-radius:5px",
                         HTML("<font color='#667266' size=3><b>COLOR IMAGES</b></font>")
                         ),
                     br(),
                     valueBox("Upload",
                              HTML("- Load the TIFF file to analyse from your local disk<br>
                                - Resize, crop or rotate your images<br>
                                     - Watch images of the sequence"),
                              icon=icon("chevron-down"), color="olive", width="60%"),
                     valueBox("Image analysis",
                              HTML("- Select an image of the sequence<br>
                                - Decompose the 3 color layers of the image<br>
                                - Plot the 3D histogram of the image<br>"),
                              icon=icon("image"), color="purple", width="60%"),
                     valueBox("Sequence analysis",
                              HTML("- Evolution of the area covered by your objects of interest<br>
                                  - 3D Evolution of the average color<br>
                                   - Donwload the evolution of the area covered by your objects of interest"),
                              icon=icon("images"), color="blue", width="60%")
                     ),
                     div(align="center",
                         style="margin-right:20px",
                         actionButton('go_color', icon("chevron-down"),
                                      style="color:#395339; border-color:white;
                                      border-bottom-left-radius:8px; border-bottom-right-radius:8px;
                                      border-top-left-radius:0px; border-top-right-radius:0px;
                                      border-width:1px;opacity:0.9;
                                      background: linear-gradient(to top, #BECDBE 0%, white 100%);",
                                      width="100%")
                     )
                     ),
              column(5,
                     div(style="border:solid;border-color:white;
                         border-top-left-radius:8px;border-top-right-radius:8px;
                         padding-left:20px;padding-right:20px;
                         margin-left:20px;margin-top:20px",
                         div(align="center",
                             style="background-color:white; opacity:0.8;
                         border-bottom-left-radius:5px;
                         border-bottom-right-radius:5px",
                             HTML("<font color='#667266' size=3><b>DEPTH IMAGES</b></font>")
                         ),
                         br(),
                     valueBox("Upload",
                              HTML("- Load the TIFF file to analyse from your local disk<br>
                                - Resize, crop or rotate your images<br>
                                     - Watch images of the sequence"),
                              icon=icon("chevron-down"), color="olive", width="50%"),
                     valueBox("Sequence analysis",
                              HTML("- Evolution of the distance from the camera to the scene<br>
                                  - Plot the 3D representation of frequency features<br>
                                   - Donwload the evolution of the average depth"),
                              icon=icon("images"), color="blue", width="50%")
                     ),
                     div(align="center",
                         style="margin-left:20px",
                         actionButton('depth_go_color', icon("chevron-down"),
                                      style="color:#395339; border-color:white;
                                      border-bottom-left-radius:8px; border-bottom-right-radius:8px;
                                      border-top-left-radius:0px; border-top-right-radius:0px;
                                      border-width:1px;opacity:0.9;
                                      background: linear-gradient(to top, #BECDBE 0%, white 100%);",
                                      width="100%")
                     )
              ),
              column(1)
      ),
      tabItem(tabName="up",
              fluidPage(
                fluidRow(
                      column(4,
                             box(title="Import the image sequence to analyse",width="100%",status="primary",solidHeader=TRUE,
                              fileInput("seq", NULL, accept = 'tif',
                                buttonLabel = HTML("<img src='https://cdn2.iconfinder.com/data/icons/circle-icons-1/64/folder-256.png' height=18></img>")),
                              htmlOutput("infos_file",
                                         style="background-color: #337ab7; border-color: #2e6da4; padding: 10px;
                                         border-radius:5px"),
                              br(),
                              numericInput("factor_resizing","Factor of resizing : ",value=1),
                              br(),
                              HTML("<font color='#001A42'>Please wait after uploading the data, until the visualisation is ready and 
                                    verify your data is well uploaded, using the visualisation of the sequence.<br>
                                   Be sure your images are encoded on 8 bits.</font>"),
                             br(),br(),br(),
                             HTML('<font color="#B94646"><div align="center"><i class="fa fa-shield-alt fa-2x"></i></div>'),
                             HTML("<p>Your data is loaded on the server when you're using the app 
                                    and will disappear when you'll leave it. We don't keep any data.</p></font>")),
                             actionButton('switchtab_up', ' Next stage', icon("chevron-circle-right"),
                                          style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                            ),
                      column(8,
                             box(width="60%",status="primary",
                                 plotOutput(outputId = "image_output",width="100%"),
                                 column(1,
                                        br(),br(),
                                        actionButton("zoom",icon("crop-alt"))
                                  ),
                                 column(11,
                                        sliderInput("tps_import","",min=1,max=10,value=1,ticks=FALSE)
                                 )
                                )
                            )
                      )
              )
      ),
      tabItem(tabName="im_ana",
              fluidRow(
                column(6,
                       box(width="100%",title="Extensive tools",status="primary",solidHeader=TRUE,
                           sliderInput("tps_decomp","",min=1,max=10,value=1,ticks=FALSE),
                           br(),
                           pickerInput("choice_seg_or_decomp",
                                       choices=c("-  Select your visualisation tool  -","3D Histogram")),
                           uiOutput("seg_or_decomp")
                           ),
                       actionButton('switchtab_ana', ' Next stage', icon("chevron-circle-right"),
                                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                       ),
                column(6,
                       box(width="100%",
                           pickerInput(
                             inputId = "choix_couche",
                             choices = c("Red", "Green", "Blue"),
                             choicesOpt = list(
                               content = sprintf("<span class='label label-%s'>%s</span>", 
                                                 c("danger","success","primary"),
                                                 c("Red","Green","Blue")))
                           ),
                           uiOutput("plots_couche"))
                       )  
              )
      ),
      tabItem(tabName="seq_ana",
              box(
                title="Plant surface evolution over time",
                uiOutput("evol2d_area_ui")
                ),
              box(
                title="3D evolution of the average color over time",
                uiOutput("ana_seq_depth")
              )
      ),
      tabItem(tabName="depth_up",
              fluidPage(
                fluidRow(
                  column(4,
                         box(title="Import the depth sequence to analyse",width="100%",status="primary",solidHeader=TRUE,
                             fileInput("depth_seq", NULL, accept = 'tif',
                                       buttonLabel = HTML("<img src='https://cdn2.iconfinder.com/data/icons/circle-icons-1/64/folder-256.png' height=18></img>")),
                             htmlOutput("depth_infos_file",
                                        style="background-color: #337ab7; border-color: #2e6da4; padding: 10px;
                                         border-radius:5px"),
                             br(),
                             numericInput("depth_factor_resizing","Factor of resizing : ",value=1),
                             br(),
                             HTML("<font color='#001A42'>Please wait after uploading the data, until the visualisation is ready and 
                                    verify your data is well uploaded, using the visualisation of the sequence.</font>"),
                             br(),br(),br(),
                             HTML('<font color="#B94646"><div align="center"><i class="fa fa-shield-alt fa-2x"></i></div>'),
                             HTML("<p>Your data is loaded on the server when you're using the app 
                                    and will disappear when you'll leave it. We don't keep any data.</p></font>")),
                         actionButton('switchtab_depth_up', ' Next stage', icon("chevron-circle-right"),
                                      style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                  ),
                  column(8,
                         box(width="60%",status="primary",
                             plotOutput(outputId = "depth_image_output",width="100%"),
                             column(1,
                                    br(),br(),
                                    actionButton("depth_zoom",icon("crop-alt"))
                             ),
                             column(11,
                                    sliderInput("depth_tps_import","",min=1,max=10,value=1,ticks=FALSE)
                             )
                         )
                  )
                )
            )
      ),
      tabItem(tabName="depth_ana",
              box(
                title="Average of the canopy height over time",
                uiOutput("mean_depth_ui") %>% withSpinner(color="#5A6C80", type=4)
              ),
              box(
                title="Plot 3D HDR",
                uiOutput("feats_3d") %>% withSpinner(color="#5A6C80", type=4)
              )
      )
    )
  )
)

server = function(input, output, session) {
  
  
  # UPLOAD ---------------------------------------------------------------------------------------------------------
  
  values <- reactiveValues(sequence=NULL,depth=NULL)
  
  observe({
    req(input$seq)
    file = input$seq
    seqim = read_tif(file$datapath)
    dimensions = dim(seqim) #H, L, couches, images
    H = dimensions[1]/input$factor_resizing
    L = dimensions[2]/input$factor_resizing
    seqim = resize(seqim,H,L)
    seqim = seqim/255
    values$sequence = round(seqim,2)
    updateSliderInput(session,"tps_import",min=1,max=dimensions[4])
    updateSliderInput(session,"tps_decomp",min=1,max=dimensions[4])
  })
  
  
  output$infos_file = renderUI({
    if (is.null(input$seq)){
      HTML("<font color=white><div align='center'>Please select the factor before uploading</div></font>")
    }
    else {
      str1 = input$seq$name
      str2 = paste("Original :",
                   format(input$seq$size,big.mark=" ",trim=FALSE),"Bytes",sep=" ")
      str3 = paste("Loaded :",
                   format(as.numeric(object.size(values$sequence))/8,big.mark=" ",trim=FALSE),"Bytes",sep=" ")
      str4 = paste("Sequence of ",dim(values$sequence)[4]," images with resolution : ",dim(values$sequence)[1],"x",dim(values$sequence)[2],sep="")
      HTML("<font color=white>",
           paste(str1,str2,str3,str4, sep="<br>"),
           "</font>")
    }
  })
  
  output$image_output = renderPlot({
    if (is.null(input$seq)){
      wait = load.image("www/waitdata.png")
      plot(wait, axes=FALSE)
    }
    else {
      par(pty="m")
      plot(as.raster(values$sequence[,,,input$tps_import]))
    }
  })
  
  observeEvent(input$zoom,{
    showModal(modalDialog(
      size="l",
      fluidRow(
        column(3, numericInput("angle","Angle rotation",value=0)),
        column(3, numericRangeInput("xaxis", "X-axis interval", c(0,ncol(values$sequence)), separator = "-")),
        column(3, numericRangeInput("yaxis", "Y-axis interval", c(0,nrow(values$sequence)), separator = "-")),
        column(3,actionButton("go_crop",HTML("CROP SEQUENCE &<br>ANALYSE"),
                              style="color: #fff; background-color: #133B5F; border-color: #2e6da4"))
      ),
      plotOutput("image_output_zoom",height=600) %>% withSpinner(color="#5A6C80", type=4),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  output$image_output_zoom = renderPlot({
    if (is.null(input$seq)){
      wait = load.image("www/waitdata.png")
      plot(wait, axes=FALSE)
    }
    else {
      par(pty="m")
      seqim_tmp = round(imrotate(values$sequence,notna(input$angle)),2)
      nr = nrow(seqim_tmp)
      nc = ncol(seqim_tmp)
      plot(as.raster(seqim_tmp[,,,notna(input$tps_import)]))
      rect(0,0,notna(input$xaxis[1]),nr,col=alpha("red",0.4),lwd=0)
      rect(0,0,nc,notna(nr-input$yaxis[2]),col=alpha("red",0.4),lwd=0)
      rect(0,notna(nr-input$yaxis[1]),nc,nr,col=alpha("red",0.4),lwd=0)
      rect(notna(input$xaxis[2]),0,nc,nr,col=alpha("red",0.4),lwd=0)
    }
  })
  
  observeEvent(input$go_crop,{
    seqim = round(imrotate(values$sequence,input$angle),2)
    D = dim(seqim)[1]
    values$sequence = seqim[input$yaxis[1]:input$yaxis[2],
                            input$xaxis[1]:input$xaxis[2],,]
    removeModal()
    updateTabItems(session,"tabs","im_ana")
  })
  # ----------------------------------------------------------------------------------------------------------------------
  
  
  # IMAGE ANALYSIS -------------------------------------------------------------------------------------------------------
  output$seg_or_decomp = renderUI({
    if (is.null(input$seq)){
      HTML('<font color="#B94646"><div align="center"><i class="fa fa-tint-slash fa-2x"></i></div>')}
    else if (input$choice_seg_or_decomp=="3D Histogram"){
      tagList(
        fluidRow(
          sliderTextInput(
            inputId = "threshold_hist3d",
            label = "Color occurences threshold applied to images :", 
            choices = 0:20,
            selected = 10,
            grid = TRUE
          ),
          style = "background-color:#E3E9F0;color:#486586;margin:10px;padding:5px;"
        ),
        fluidRow(
          plotlyOutput("decomp",height=400) %>% withSpinner(color="#5A6C80", type=4),
          style = "background-color:#E3E9F0;color:#486586;margin:10px;"
        )
      )
    }
    else if (input$choice_seg_or_decomp=="Plants segmentation") {
      tagList(
        plotOutput("segmentation",height=400),
        downloadButton("download_mask", "Download image mask",
                       style="color: #fff; background-color: #588E74; border-color: #2e6da4")
      )}
  })
  
  output$decomp = renderPlotly({
    decompimgRGB(values$sequence[,,,input$tps_decomp],input$threshold_hist3d)
  })
  
  output$plots_couche = renderUI({
    if (is.null(input$seq)){
      HTML('<font color="#B94646"><div align="center"><i class="fa fa-tint-slash fa-2x"></i></div>')}
    else{
      tagList(
        plotOutput("couche_col"),
        plotOutput("hist_couche_col",height=250)
      )}
  })
  
  output$couche_col = renderPlot({
    plot(as.raster(values$sequence[,,as.integer(dic_col[input$choix_couche]),input$tps_decomp]))
  })
  
  output$hist_couche_col = renderPlot({
    decompimgNG(values$sequence[,,as.integer(dic_col[input$choix_couche]),input$tps_decomp],input$choix_couche)
  })
  
  output$segmentation = renderPlot({
    plot(as.raster(segmentation(values$sequence[,,,input$tps_decomp],FALSE)))
  })
  
  output$download_mask = downloadHandler(
    filename = paste("my_mask_for_image",input$tps_decomp,".tif",sep=""),
    content = function(file) {
      seg = segmentation(values$sequence[,,,input$tps_decomp])
      seg[seg] = 255
      write_tif(seg,
                file)
    }
  )
  
  # ---------------------------------------------------------------------------------------------------------------------
  
  
  # SEQUENCE ANALYSIS ---------------------------------------------------------------------------------------------------
  output$evol2d_area_ui = renderUI({
    if (is.null(input$seq)){
      HTML('<font color="#B94646"><div align="center"><i class="fa fa-tint-slash fa-2x"></i></div>')}
    else{
      tagList(
        plotOutput("evol2d_area")  %>% withSpinner(color="#5A6C80", type=4),
        downloadButton("download_mean_evol","Download evolution vector",
                       style="color: #fff; background-color: #588E74; border-color: #2e6da4")
      )
    }
  })
  
  output$evol2d_area = renderPlot({
    N = dim(values$sequence)[4]
    list_seg = apply(values$sequence,4,segmentation)
    evol = apply(list_seg,2,sum)
    D = data.frame(Time=1:N,evol)
    ggplot(D,aes(x=Time)) + geom_line(aes(y=evol,col=evol),alpha=0.6,lwd=1.5) +
      labs(y="Estimated plant surface (pixels)")+ guides(color=FALSE) +
      scale_color_gradient(low="#6A93B7",high="#216199") +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))# + theme_ipsum()  # + theme_ipsum() 
  })
  
  output$ana_seq_depth = renderUI({
    if (is.null(input$seq)){
      HTML('<font color="#B94646"><div align="center"><i class="fa fa-tint-slash fa-2x"></i></div>')}
    else {
      tagList(
        pickerInput("depth_seq_choice",
                    choices=c("3D evolution of the average color","Evolution of the average color intensity","Difference of growth between 2 times")),
        uiOutput("ana_seq_depth_plot") %>% withSpinner(color="#5A6C80", type=4)
      )
    }
  })
  
  
output$ana_seq_depth_plot = renderUI({
  if (input$depth_seq_choice=="3D evolution of the average color"){
    plotlyOutput("evol3d_meancol")
  }
  else if (input$depth_seq_choice=="Difference of growth between 2 times"){
    tagList(
      sliderInput("times_seg", label = "Times to compare", min = 1, ticks=FALSE, step=1,
                  max = dim(values$sequence)[4], value = c(1,dim(values$sequence)[4])),
      plotOutput("segmentation_twice"),
      fluidRow(
        column(1),
        column(3,
               valueBoxOutput("seg_twice_stats_from",width="100%")),
        column(4,
               infoBoxOutput("seg_twice_stats",width="100%")),
        column(3,
               valueBoxOutput("seg_twice_stats_to",width="100%"))
      )
    )
  }
  else {
	plotOutput("evol2d_color")
  }
})

  output$evol2d_color = renderPlot({
    N = dim(values$sequence)[4]
    D = data.frame(Time=1:N)
    for (n_col in 1:3){
    	seq = values$sequence[,,n_col,]
    	evol = apply(seq,3,mean)
    	D = cbind(D,evol)
	}
	colnames(D) = c("Time","evolR","evolG","evolB")
    ggplot(D,aes(x=Time)) +
    	geom_line(aes(y=evolR),alpha=0.6,lwd=1,col="red") +
    	geom_line(aes(y=evolG),alpha=0.6,lwd=1,col="green") +
    	geom_line(aes(y=evolB),alpha=0.6,lwd=1,col="blue") +
    	ylim(0,1) +
      	labs(y="Average intensity of color")+ guides(color=FALSE) +
      	theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))# + theme_ipsum()  # + theme_ipsum() 
  })

  output$seg_twice_stats_from <- renderValueBox({
    t1 = input$times_seg[1]
    X1 = sum(segmentation(values$sequence[,,,t1],TRUE))
    valueBox(HTML("<font color='red'>Before</font>"),
             HTML(paste("<b>",X1,"pixels</b>")),
             color = "light-blue")
  })
  
  output$seg_twice_stats_to <- renderValueBox({
    t1 = input$times_seg[2]
    X1 = sum(segmentation(values$sequence[,,,t1],TRUE))
    valueBox(HTML("<font color='darkblue'>After</font>"),
             HTML(paste("<b>",X1,"pixels</b>")),
             color = "light-blue")
  })

  output$seg_twice_stats <- renderInfoBox({
    t1 = input$times_seg[1]
    t2 = input$times_seg[2]
    X1 = sum(segmentation(values$sequence[,,,t1],TRUE))
    X2 = sum(segmentation(values$sequence[,,,t2],TRUE))
    taux = round(100*(X2-X1)/X1)
    if (taux<(-25)){views = c("angle-double-down","red")}
    else if (taux<0){views = c("angle-down","orange")}
    else if (taux<25){views = c("angle-up","yellow")}
    else {views = c("angle-double-up","green")}
    infoBox("Progress", paste(signe(taux),abs(taux),"%"), icon = icon(views[1]), color = views[2], fill=TRUE)
  })
  
  output$segmentation_twice = renderPlot({ # PENSER A METTRE UNE LEGENDE
    t1 = input$times_seg[1]
    t2 = input$times_seg[2]
    X1 = segmentation(apply(values$sequence[,,,t1],c(1,2),sum),TRUE)
    X2 = segmentation(apply(values$sequence[,,,t2],c(1,2),sum),TRUE)
    M = array(0,dim(values$sequence)[1:3])
    M[,,1][X1] = 1
    M[,,3][X2] = 1
    plot(as.raster(M))
  })
  
  output$segmentation_twice_explain = renderUI({
    space = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
    str1 = "<font color='red'>Before</font> : &nbsp;&nbsp;&nbsp; Time"
    str2 = "<font color='blue'>After</font> : &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Time"
    HTML(paste(space,str1,input$times_seg[1],"<br>",space,str2,input$times_seg[2]))
  })
  
  output$evol3d_meancol = renderPlotly({
    all = data.frame()
    seg = data.frame()
    N = dim(values$sequence)[4]
    for (i in 1:N){
      all = rbind(all,apply(values$sequence[,,,i],3,mean))
      tmp = apply(values$sequence[,,,i],c(1,2),sum)
      mask = segmentation(tmp,TRUE)
      M = c()
      for (j in 1:3){
        M = c(M,mean(values$sequence[,,j,i][mask]))
      }
      seg = rbind(seg,M)
    }
    all = cbind(round(all,2),0)
    seg = cbind(round(seg,2),1)
    colnames(all) = c("R","G","B","type")
    colnames(seg) = c("R","G","B","type")
    X = rbind(all,seg)
    col = apply(X[,1:3],1,RVB)
    X = cbind(rep(1:N,2),X,col)
    m = min(X[2:4])
    M = max(X[2:4])
    colnames(X) = c("Time","R","G","B","type","col")
    X$type = as.factor(X$type)
    levels(X$type) = c('All','Segmented')
    plot_ly(X,
            x=~R,y=~G,z=~B,
            color = ~type, colors = c('coral3', 'darkblue'),
            type = 'scatter3d', mode = 'lines',
            text=~paste("Time",Time,"-",type),
            hoverinfo = 'text',
            marker = list(size = 4,
                          color = ~col,
                          opacity = 0.5)
            )  %>%
      layout(scene = list(xaxis = list(title = 'Red', range=c(m,M)),
                          yaxis = list(title = 'Green', range=c(m,M)),
                          zaxis = list(title = 'Blue',range=c(m,M)))
             )
  })
  
  output$download_mean_evol = downloadHandler(
    filename = "area_evolution.csv",
    content = function(file) {
      write.csv(100*apply(apply(values$sequence,4,segmentation),2,mean),
                file,
                row.names=FALSE)
    }
  )
  # ---------------------------------------------------------------------------------------------------------------------
  
  # DEPTH UPLOAD --------------------------------------------------------------------------------------------------------
  observe({
    req(input$depth_seq)
    file = input$depth_seq
    seqim = read_tif(file$datapath)
    dimensions = dim(seqim) #H, L, couches, images
    H = dimensions[1]/input$depth_factor_resizing
    L = dimensions[2]/input$depth_factor_resizing
    seqim = resize(seqim,H,L)
    seqim = seqim/255
    values$depth = round(seqim,2)
    updateSliderInput(session,"depth_tps_import",min=1,max=dimensions[4])
    updateSliderInput(session,"depth_tps_decomp",min=1,max=dimensions[4])
  })
  
  
  output$depth_infos_file = renderUI({
    if (is.null(input$depth_seq)){
      HTML("<font color=white><div align='center'>Please select the factor before uploading</div></font>")
    }
    else {
      str1 = input$depth_seq$name
      str2 = paste("Original :",
                   format(input$depth_seq$size,big.mark=" ",trim=FALSE),"Bytes",sep=" ")
      str3 = paste("Loaded :",
                   format(as.numeric(object.size(values$depth))/8,big.mark=" ",trim=FALSE),"Bytes",sep=" ")
      str4 = paste("Sequence of ",dim(values$depth)[4]," images with resolution : ",dim(values$depth)[1],"x",dim(values$depth)[2],sep="")
      HTML("<font color=white>",
           paste(str1,str2,str3,str4, sep="<br>"),
           "</font>")
    }
  })
  
  output$depth_image_output = renderPlot({
    if (is.null(input$depth_seq)){
      wait = load.image("www/waitdata.png")
      plot(wait, axes=FALSE)
    }
    else {
      par(pty="m")
      plot(as.raster(values$depth[,,,input$tps_import]))
    }
  })
  
  observeEvent(input$depth_zoom,{
    showModal(modalDialog(
      size="l",
      fluidRow(
        column(3, numericInput("depth_angle","Angle rotation",value=0)),
        column(3, numericRangeInput("depth_xaxis", "X-axis interval", c(0,ncol(values$depth)), separator = "-")),
        column(3, numericRangeInput("depth_yaxis", "Y-axis interval", c(0,nrow(values$depth)), separator = "-")),
        column(3, actionButton("depth_go_crop",HTML("CROP SEQUENCE &<br>ANALYSE"),
                               style="color: #fff; background-color: #133B5F; border-color: #2e6da4"))
      ),
      plotOutput("depth_image_output_zoom",height=600) %>% withSpinner(color="#5A6C80", type=4),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  output$depth_image_output_zoom = renderPlot({
    if (is.null(input$depth_seq)){
      wait = load.image("www/waitdata.png")
      plot(wait, axes=FALSE)
    }
    else {
      par(pty="m")
      seqim_tmp = imrotate(values$depth,input$depth_angle)
      nr = nrow(seqim_tmp)
      nc = ncol(seqim_tmp)
      plot(as.raster(seqim_tmp[,,,input$depth_tps_import]))
      rect(0,0,notna(input$depth_xaxis[1]),nr,col=alpha("red",0.4),lwd=0)
      rect(0,0,nc,notna(nr-input$depth_yaxis[2]),col=alpha("red",0.4),lwd=0)
      rect(0,notna(nr-input$depth_yaxis[1]),nc,nr,col=alpha("red",0.4),lwd=0)
      rect(notna(input$depth_xaxis[2]),0,nc,nr,col=alpha("red",0.4),lwd=0)
    }
  })
  
  observeEvent(input$depth_go_crop,{
    seqim = imrotate(values$depth,input$depth_angle)
    values$depth = as_ijtiff_img(seqim[input$depth_yaxis[1]:input$depth_yaxis[2],
                                       input$depth_xaxis[1]:input$depth_xaxis[2],,])
    removeModal()
    updateTabItems(session,"tabs",selected="depth_ana")
  })
  
  # ------------------------------------------------------------------------------------------------------------------
  
  # DEPTH ANALYSIS ---------------------------------------------------------------------------------------------------
  output$mean_depth_ui = renderUI({
    if (is.null(input$depth_seq)){
      HTML('<font color="#B94646"><div align="center"><i class="fa fa-tint-slash fa-2x"></i></div>')}
    else{
      tagList(
        plotOutput("evol_depth")  %>% withSpinner(color="#5A6C80", type=4),
        downloadButton("download_mean_depth","Download evolution vector",
                       style="color: #fff; background-color: #588E74; border-color: #2e6da4")
      )
    }
  })
  
  output$evol_depth = renderPlot({
    N = dim(values$depth)[4]
    evol = apply(values$depth,4,mean)
    v_lines = seq(1,N+1,by=input$periode_feats)
    D = data.frame(Time=1:N,evol)
    ggplot(D,aes(x=Time)) + geom_line(aes(y=evol,col=evol),alpha=0.6,lwd=1.5) +
      labs(y="Average pixels intensity (distance to camera)") + guides(color=FALSE) +
      scale_color_gradient(low="#6A93B7",high="#216199") +
      geom_vline(xintercept = v_lines, col="royalblue", alpha=0.75, linetype="dashed") +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))# + theme_ipsum() 
  })
  
  output$download_mean_depth = downloadHandler(
    filename = "depth_sequence_evolution.csv",
    content = function(file) {
      write.csv(apply(values$depth,4,mean),
                file,
                row.names=FALSE)
    }
  )
  
  output$feats_3d = renderUI({
    if (is.null(input$depth_seq)){
      HTML('<font color="#B94646"><div align="center"><i class="fa fa-tint-slash fa-2x"></i></div>')}
    else{
      tagList(
        numericInput("periode_feats","Indicate the periodicity",96),
        uiOutput("feats_3d_rpz")  %>% withSpinner(color="#5A6C80", type=4),
        htmlOutput("explain_feats")
      )
    }
  })
  
  output$feats_3d_rpz = renderUI({
    if (dim(values$depth)[4]<input$periode_feats){
      htmlOutput("arg_feats")
    }
    else {
      plotlyOutput("feats_3d_plot")
    }
  })
  
  output$feats_3d_plot = renderPlotly({
    N = dim(values$depth)[4]
    evol = apply(values$depth,4,mean)
    evol3d_feats(evol,input$periode_feats,0,"blue")
  })
  
  output$explain_feats = renderUI({
    HTML("Color lines depend on the harmonic for each period<br>
         <b>Average harmonic on periods :</b>
         <font color='Red'>red</font> 
         is maximum")
  })
  
  output$arg_feats = renderUI({
    div(align="center",HTML("<font color='red'>Period lower than the number of values</font>"))
  })
  
  
  # ------------------------------------------------------------------------------------------------------------------
  
  # SWITCH TAB
  observeEvent(input$go_color, {
    updateTabItems(session, "tabs", "up")
  })
  observeEvent(input$depth_go_color, {
    updateTabItems(session, "tabs", "depth_up")
  })
  observeEvent(input$switchtab_up, {
    updateTabItems(session, "tabs", "im_ana")
  })
  observeEvent(input$switchtab_ana, {
    updateTabItems(session, "tabs", "seq_ana")
  })
  observeEvent(input$switchtab_depth_up, {
    updateTabItems(session, "tabs", "depth_ana")
  })


}


#runApp()
#shinyAppDir(".")
shinyApp(ui,server)
#runApp(host='127.0.0.1', port=2000)
