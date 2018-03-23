library(shiny)
library(ggplot2)
library(data.table)
library(DESeq2)

source("../common/ggplot_ext/eupath_default.R")
source("../common/ggplot_ext/eupath_functions.R")
source("../common/tooltip/abundance_tt.R")
source("../common/mbiome/mbiome-reader.R")
source("../common/mbiome/mbiome-stats.R")
source("../common/mbiome/mbiome-utils.R")
source("../common/config.R")

shinyServer(function(input, output, session) {
  
  mstudy_obj<-NULL
  
  
  df_abundance <- NULL
  df_sample <- NULL
  df_sample.formatted <- NULL
  
  # global objects to read in more than one function
  columns <- NULL
  hash_sample_names<- NULL
  hash_count_samples <- NULL
  ggplot_object <- NULL
  ggplot_build_object <- NULL
  
  SAMPLE <- NULL
  OTU<-NULL
  TAX<-NULL
  
  loaded_category <- FALSE
  
  hash_colors <- NULL
  
  abundance_otu <- NULL
  abundance_otu_relative <- NULL
  abundance_taxa <- NULL
  
  maximum_samples_without_resizing <- 65
  minimum_height_after_resizing <- 6
  
  NO_METADATA_SELECTED <- "No Metadata Selected"
  
  load_microbiome_data <- reactive({
    if(is.null(mstudy_obj)){
      mstudy_obj <<- import.biom(biom_file, metadata_details, aggregate_by = input$taxonLevel)
      
      filtered_categories <- mstudy_obj$get_filtered_categories()
      
      selected_category<-which.min(mstudy_obj$sample_table$get_count_category())
      
      updateSelectizeInput(session, "category",
                           choices = c(filtered_categories),
                           selected = filtered_categories[selected_category],
                           options = list(placeholder = 'Choose the metadata to calculate the differential abundance'),
                           server = TRUE)
      
      unique_values<-mstudy_obj$sample_table$get_unique_details(filtered_categories[selected_category])
      
      if(length(unique_values)>2){
              choices_level <- c(paste("not",unique_values[1]), unique_values)
              selected_level<-choices_level[1]
              updateSelectizeInput(session, "factor1",
                                   choices = c("Not Factor 2", unique_values),
                                   selected = unique_values[1],
                                   options = list(placeholder = 'Choose the first factor'),
                                   server = T)
              updateSelectizeInput(session, "factor2",
                                   choices = c("Not Factor 1", unique_values),
                                   selected = unique_values[2],
                                   options = list(placeholder = 'Choose the second factor'),
                                   server = T)

            }else{
              updateSelectizeInput(session, "factor1",
                                   choices = unique_values,
                                   selected = unique_values[1],
                                   options = list(placeholder = 'Choose the first factor'),
                                   server = T)
              updateSelectizeInput(session, "factor2",
                                   choices = unique_values,
                                   selected = unique_values[2],
                                   options = list(placeholder = 'Choose the second factor'),
                                   server = T)
            }
      
      
      
    }
    
    mstudy_obj
  })
  
  observeEvent(input$exchangeBtn, {
    isolate(factor1<-input$factor1)
    isolate(factor2<-input$factor2)
    if(identical(factor1, "Not Factor 2")){
      factor1<-"Not Factor 1"
    }
    if(identical(factor2, "Not Factor 1")){
      factor2<-"Not Factor 2"
    }
    updateSelectizeInput(session, "factor1", selected = factor2, server = F)
    updateSelectizeInput(session, "factor2", selected = factor1, server = F)
  })

  observeEvent(input$category, {
    category<-input$category
    if(!identical(category, "")){
      
      unique_values<-mstudy_obj$sample_table$get_unique_details(category)

      if(length(unique_values)>2){
        updateSelectizeInput(session, "factor1",
                             choices = c("Not Factor 2", unique_values),
                             selected = unique_values[1],
                             options = list(placeholder = 'Choose the first factor'),
                             server = F)
        updateSelectizeInput(session, "factor2",
                             choices = c("Not Factor 1", unique_values),
                             selected = unique_values[2],
                             options = list(placeholder = 'Choose the second factor'),
                             server = F)
      }else{
        updateSelectizeInput(session, "factor1",
                             choices = unique_values,
                             selected = unique_values[1],
                             options = list(placeholder = 'Choose the first factor'),
                             server = F)
        updateSelectizeInput(session, "factor2",
                             choices = unique_values,
                             selected = unique_values[2],
                             options = list(placeholder = 'Choose the second factor'),
                             server = F)
      }
    }
  })
  
  mainChart <- function(){}
  output$mainChart <- renderUI({
    mstudy_obj <- load_microbiome_data()
    
    isolate(category<-input$category)
    
    factor1<-input$factor1
    factor2<-input$factor2
    taxon_level<-input$taxonLevel
    if(!identical(category, "") & !identical(factor1, "") & !identical(factor2, "")){
      if(identical(factor1, factor2) | (identical(factor1, "Not Factor 2") & identical(factor2, "Not Factor 1") )){
        output$datatableOutput<-renderDataTable(NULL)
        shinyjs::disable("btnDownloadPNG")
        shinyjs::disable("btnDownloadSVG")
        shinyjs::disable("btnDownloadEPS")
        shinyjs::disable("btnDownloadCSV")
        return(
          h5(class="alert alert-danger", "Please choose different factors to calculate the differential abundance.")  
        )
      }else{
        deseq_result <- runDeseq()
        
        if(!is.null(deseq_result)){
          # the x axis will be from -max_fold_change to +max_fold_change
          max_fold_change<-max(abs(deseq_result$log2FoldChange))
          # category_column <- hash_count_samples[[category]]
          
          limits_plot<-c(levels(deseq_result[[taxon_level]]),"","")
          
          chart<-ggplot(deseq_result, aes_string(x="log2FoldChange", y=taxon_level, color="Phylum")) +
            annotate("text", x = 0, y = nrow(deseq_result)+1, size=5, colour = "red", parse=T,
                     label = sprintf("atop(bold(\"%s - %s vs %s\"))", category, factor2, factor1) )+
            xlim(-max_fold_change, max_fold_change)+
            # ylim(0, nrow(deseq_result)+1.5)+
            geom_segment(aes_string(x = 0, y = taxon_level, xend = "log2FoldChange", yend = taxon_level), color = "grey50") +
            geom_point(aes(size = log10(0.5+1/pvalue)))+
            scale_size(range = c(3, 9), guide = 'none')+
            # theme_eupath(legend.position = "right", legend.direction="vertical")+
            theme_eupath(legend.position = "bottom")+
            scale_y_discrete(position = "right", limits=limits_plot, breaks = levels(deseq_result[[taxon_level]]) )+
            guides(colour = guide_legend(override.aes = list(size=8)))
          
          ggplot_object<<-chart
          
          ggplot_build_object<<-ggplot_build(chart)
          shinyjs::enable("btnDownloadPNG")
          shinyjs::enable("btnDownloadSVG")
          shinyjs::enable("btnDownloadEPS")
          shinyjs::enable("btnDownloadCSV")
          
          output$plotWrapper<-renderPlot({
            chart
          })
          
          result_to_show<-plotOutput("plotWrapper", hover = hoverOpts("plot_hover"))
          
        }else{ 
          # output$datatableOutput<-renderDataTable(NULL)
          # output$mainChart<-renderUI(NULL)
          # ggplot_object<<-NULL
          # ggplot_build_object<<-NULL
          shinyjs::disable("btnDownloadPNG")
          shinyjs::disable("btnDownloadSVG")
          shinyjs::disable("btnDownloadEPS")
          shinyjs::disable("btnDownloadCSV")
          result_to_show<-h5(class="alert alert-warning", "Sorry, but there is no OTU with differential abundance using your search parameters.")
        }
        result_to_show
      }
    }else{
      if(!identical(category, "") & (identical(factor1,"")|identical(factor2,"") )){
        output$datatableOutput<-renderDataTable(NULL)
        ggplot_object<<-NULL
        ggplot_build_object<<-NULL
      }
    }
  })
  
  runDeseq <- reactive({
    isolate(
      category<-input$category
    )

    factor1<-input$factor1
    factor2<-input$factor2
    taxon_level <- input$taxonLevel

    chart<-NULL
    shinyjs::hide("chartContent")
    shinyjs::show("chartLoading")
    if(!identical(factor1, factor2) &  !(identical(factor1, "Not Factor 2") & identical(factor2, "Not Factor 1") ) ){
      # category_column <- hash_sample_names[[hash_count_samples[[category]]]]
      
      # df_sample_filter <- SAMPLE[,category_column]
      df_sample_filter <- mstudy_obj$get_metadata_as_column(category)
      # print(head(OTU))
      # print(head(TAX))
      # print(head(df_sample_filter))
      if(identical(factor2, "Not Factor 1")){
        # df_sample_filter[df_sample_filter[,category]!=factor1,category]<-"Not Factor 1"
        df_sample_filter[get(category)!=factor1,category]<-"Not Factor 1"
      }else if(identical(factor1, "Not Factor 2")){
        df_sample_filter[get(category)!=factor2,category]<-"Not Factor 2"
      }else{
        df_sample_filter<-df_sample_filter[get(category)==factor1 | get(category)==factor2]
        # df_sample_filter<-df_sample_filter[df_sample_filter[,category_column]==factor1 | df_sample_filter[,category_column]==factor2]
      }
      df_sample_filter<-as.data.frame(df_sample_filter)
      df_sample_filter[[category]] <- factor(df_sample_filter[[category]], levels=c(factor1,factor2))
      rownames(df_sample_filter)<-df_sample_filter$SampleName
      
      dcasted_otu <- mstudy_obj$otu_table$get_sample_as_column_by_otu(taxon_level)
      selected_levels<-get_columns_taxonomy(taxon_level)
      dt_taxonomy<-dcasted_otu[,selected_levels,with=F]
      dt_otu<-as.data.frame(dcasted_otu[,df_sample_filter$SampleName,with=F])
      
      diagdds <- DESeqDataSetFromMatrix(countData=dt_otu, colData=df_sample_filter, as.formula(paste0("~", category)))
      

      # diagdds = phyloseq2Deseq2(new_physeq_obj, as.formula(paste0("~", category_column)))
      # calculate geometric means prior to estimate size factors
      gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }

      geoMeans = apply(counts(diagdds), 1, gm_mean)
      diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
      diagdds = DESeq(diagdds, fitType="local", quiet=T)
      #investigating test result
      res = results(diagdds)
      rownames(res)<-1:nrow(res)
      # res = res[order(res$padj, na.last=NA), ]
      alpha = 0.01
      sigtab = res[!is.na(res$padj) & (res$padj < alpha), ]
      # print("helloooo")
      if(nrow(sigtab)>0){
        sigtab = cbind(as(sigtab, "data.frame"), as(dt_taxonomy[as.integer(rownames(sigtab)), ], "matrix"))
        # sigtab = cbind(as(sigtab, "data.frame"), as(TAX[rownames(sigtab), ], "matrix"))

        sigtab<-sigtab[order(sigtab$log2FoldChange),]
        sigtab[,taxon_level]<-factor(sigtab[,taxon_level], levels=sigtab[,taxon_level])

        cols_to_show<-sigtab[,c("baseMean", "log2FoldChange", "pvalue", taxon_level)]

        output$datatableOutput<-renderDataTable(cols_to_show,
                                                options = list(
                                                  order = list(list(1, 'desc'))
                                                )
        )
        # log data transformation by https://stats.stackexchange.com/questions/83914/how-to-log-transform-data-with-a-large-number-of-zeros
        chart<-sigtab
      }else{
        output$datatableOutput<-renderDataTable(NULL)
      }
    }
    shinyjs::hide("chartLoading")
    shinyjs::show("chartContent")
    chart
  })
  


  hovers <- function(){}

  output$hover_info <- renderUI({

    hover <- input$plot_hover

    if (is.null(hover$x) || is.null(hover$y) || is.null(ggplot_object))
      return(NULL)

    left_pct <-
      (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <-
      (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)

    left_px <-
      hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <-
      hover$range$top + top_pct * (hover$range$bottom - hover$range$top)

    style <-
      paste0(
        "position:absolute; background-color: rgba(245, 245, 245, 0.85); z-index:1000;",
        "left:",
        left_px + 2,
        "px; top:",
        top_px + 2,
        "px;"
      )
    near_points <- nearPoints(ggplot_object$data, hover)
    if(nrow(near_points)>0){
      isolate(taxa_level<-input$taxonLevel)
      isolate(category<-input$category)
      isolate(factor1<-input$factor1)
      isolate(factor2<-input$factor2)

      text_hover<-""
      if(nrow(near_points) == 1){
        
        # filtered_taxa <- abundance_taxa[abundance_taxa[,taxa_level] == near_points[1, taxa_level],]
        # filtered_otu <- abundance_otu_relative[rownames(filtered_taxa),]
      
        dcasted_otu <- mstudy_obj$otu_table$get_sample_as_column_by_otu(taxon_level)
        dcasted_otu <- dcasted_otu[as.integer(rownames(near_points)),]
        selected_levels<-get_columns_taxonomy(taxon_level)
        
        dt_otu<-as.data.frame(dcasted_otu[,mstudy_obj$get_sample_names(),with=F])
        
        otu_for_plot <- as.data.frame(t(dt_otu))
        otu_for_plot$SampleName <- rownames(otu_for_plot)
        
        otu_for_plot_filtered<-otu_for_plot[otu_for_plot[,1]>0,]

        colnames(otu_for_plot_filtered)<-c("Abundance", "SampleName")

        # category_column<-hash_sample_names[[hash_count_samples[[category]]]]
        # df_sample_selected <- df_sample.formatted[,c("SampleName", category_column)]
        df_sample_selected<-mstudy_obj$get_metadata_as_column(category)
        
        if(identical(factor2, "Not Factor 1")){
          # df_sample_filter[df_sample_filter[,category]!=factor1,category]<-"Not Factor 1"
          df_sample_selected[get(category)!=factor1,category]<-"Not Factor 1"
        }else if(identical(factor1, "Not Factor 2")){
          df_sample_selected[get(category)!=factor2,category]<-"Not Factor 2"
        }else{
          df_sample_selected<-df_sample_filter[get(category)==factor1 | get(category)==factor2]
          # df_sample_filter<-df_sample_filter[df_sample_filter[,category_column]==factor1 | df_sample_filter[,category_column]==factor2]
        }

        data_merged <- merge(df_sample_selected, otu_for_plot_filtered, by = "SampleName", all.y=F)
        
        data_merged[[category]] <- factor(data_merged[[category]], levels=c(factor1,factor2))

        chart<-ggplot(data_merged, aes_string(x=category, y="Abundance"))+geom_boxplot()+
          theme(
            axis.title = element_text(family = "Palatino", color="black", face="bold", size=16),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            text = element_text(family = "Palatino", size=13, face="bold", color="black"),
            panel.border = element_rect(colour="black", size=1, fill=NA),
            strip.text.x = element_text(family = "Palatino", size=13, face="bold", color="black"),
            strip.background = element_rect(fill="#F3F2F2")
          )+
          labs(x=stringi::stri_trans_totitle(hash_count_samples[[category]]), y="Relative Abundance")

        output$plotHover<-renderPlot(chart)

        text_hover <- paste0(text_hover, sprintf("<b>%s: </b>", taxa_level), near_points[1, taxa_level],
                       "<br><b>log2FoldChange: </b>", sprintf("%.3f",near_points[1,"log2FoldChange"]),
                       "<br><b>pvalue: </b>", sprintf("%0.6g",near_points[1,"pvalue"]) )
      }else{
        # for(i in 1:nrow(near_points)){
        # }
      }
      wellPanel(style = style,
        HTML(text_hover),
        plotOutput("plotHover", width = 250, height = 200)
      )
    }else{
      return(NULL)
    }

  })


  downloads <- function(){}

  output$btnDownloadPNG <- downloadHandler(
    filename = "plot.png",
    content = function(file) {
      png(file, width=1200,height=800,units="px")
      print(ggplot_object)
      dev.off()
    }
  )

  output$btnDownloadEPS <- downloadHandler(
    filename = "plot.eps",
    content = function(file) {
      setEPS()
      postscript(file, width=16,height=10.67, family = "Helvetica")
      print(ggplot_object)
      dev.off()
    }
  )

  output$btnDownloadSVG <- downloadHandler(
    filename = "plot.svg",
    content = function(file) {
      svg(file, width=16,height=10.67)
      print(ggplot_object)
      dev.off()
    }
  )

  output$btnDownloadCSV <- downloadHandler(
    filename = "data.csv",
    content = function(file) {
      write.csv(ggplot_object$data, file)
    }
  )
  
  shinyjs::hide(id = "loading-content", anim = TRUE, animType = "fade")
  shinyjs::show("app-content")
})