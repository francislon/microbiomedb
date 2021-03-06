# check biom version

# Declaring the packages
library(shiny)
library(ggplot2)
library(data.table)
source("functions.R")
source("../common/ggplot_ext/facet_even.R")
source("../common/ggplot_ext/eupath_default.R")
source("../common/tooltip/abundance_tt.R")
source("../common/mbiome/mbiome-reader.R")
source("../common/config.R")


shinyServer(function(input, output, session) {

  # Declaring global variables
  # df_abundance, df_sample and df_sample.formatted are declared global to avoid
  # multiple file reading in the reactive section
  mstudy_obj <- NULL
  hash_colors <- NULL

  global_otus<-NULL 
  selected_levels<-NULL
  # variables to define some plot parameters
  NUMBER_TAXA <- 10
  WIDTH <- 930
  MAX_SAMPLES_NO_RESIZE <- 40
  MIN_HEIGHT_AFTER_RESIZE <- 9.5

  NO_METADATA_SELECTED <- "Click to select sample details"

  ggplot_object <- NULL
  ggplot_build_object <- NULL

  ggplot_by_top_otu_object <- NULL
  ggplot_build_by_top_otu_object <- NULL

  ggplot_by_otu_object <- NULL
  ggplot_build_by_otu_object <- NULL

  eupath_pallete<-c("#999999", "#6a3d9a", "#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")

  hash_colors<-NULL

  load_microbiome_data <- reactive({
    if(is.null(mstudy_obj)){

      mstudy_obj <<- import.biom(biom_file, metadata_details)

      updateSelectizeInput(session, "category",
                           choices = c(NO_METADATA_SELECTED,
                                       mstudy_obj$get_filtered_categories()),
                           selected = NO_METADATA_SELECTED,
                           options = list(placeholder = 'Choose metadata to split the chart'),
                           server = TRUE)
    }
    mstudy_obj
  })

  observeEvent(input$taxonLevel,{
    mstudy<-load_microbiome_data()

    taxon_level <- input$taxonLevel
    if(!identical(taxon_level, "")){
      output$overviewDatatable <- renderDataTable(NULL)
      global_otus <<- mstudy$get_otus_by_level(taxon_level)
      mstudy$calculate_relative_abundance()
      selected_levels <<- get_columns_taxonomy(taxon_level)
      hash_colors <<- eupath_pallete
      top_ten<-mstudy$get_top_n_by_mean(taxon_level, NUMBER_TAXA)
      ordered<-c(mstudy$otu_table$get_ordered_otu(NUMBER_TAXA))
      rev_ordered<-c("Other", rev(ordered))
      names(hash_colors) <<- rev_ordered

      updateSelectizeInput(session, "filterOTU",
                           choices = global_otus,
                           selected = global_otus[1],
                           options = list(placeholder = 'Choose a OTU'),
                           server = TRUE)
    }
  })

  overviewChart <- function(){}

  output$overviewChart <- renderUI({
    mstudy <- load_microbiome_data()
    category<-input$category
    taxon_level<-input$taxonLevel
    result_to_show<-NULL
    is_numeric <- F
    if(!identical(input$category, "")){
      shinyjs::hide("divContent")
      shinyjs::show("chartLoading")

      # set_parameters()

      quantity_samples <- mstudy$get_sample_count()
      #TODO test for input without 10 different taxa
      top_ten<-mstudy$get_top_n_by_mean(taxon_level, NUMBER_TAXA)

      ordered<-c(mstudy$otu_table$get_ordered_otu(NUMBER_TAXA))

      rev_ordered<-c("Other", rev(ordered))

      wrapped_labels<-wrap_column(rev_ordered)

      top_ten[[taxon_level]]<-factor(top_ten[[taxon_level]], levels=rev_ordered)

      if(identical(category, NO_METADATA_SELECTED)){
        chart<-ggplot(top_ten, aes_string(x="SampleName", y="Abundance", fill=taxon_level))+
          geom_bar(stat="identity", position="stack", color="black")+
          theme_eupath_default(
            legend.title.align=0.4,
            legend.title = element_text(colour="black", size=rel(1), face="bold"),
            axis.text.y=element_blank(),
            axis.ticks.y = element_blank()
          )+
          scale_fill_manual(values=eupath_pallete, name=taxon_level,
                            labels = c(wrapped_labels),
                            guide = guide_legend(reverse=TRUE, keywidth = 1.7, keyheight = 1.7)
                          )+
          labs(x="Samples", y="Phylogenetic Relative Abundance")+
          coord_flip(expand=F)
      }else{
        dt_metadata<-mstudy$get_metadata_as_column(category)
        # dt_metadata<-merge(dt_metadata, top_ten, by="SampleName")

        if( is.numeric(dt_metadata[[category]])){
          is_numeric<-T
          samples<-mstudy$get_sample_names()
          combinations<-data.table(expand.grid(samples, rev_ordered))

          combinations[,Abundance:=0]
          colnames(combinations)<-c("SampleName", taxon_level, "Abundance")

          merged<-merge(combinations, top_ten, by=c("SampleName", taxon_level), all.x=T)

          selected_levels <- get_columns_taxonomy(taxon_level)

          merged[,c("Abundance.x", selected_levels[1:(length(selected_levels)-1)]):=NULL]
          colnames(merged)<-c("SampleName", taxon_level, "Abundance")

          merged[is.na(Abundance),Abundance:=0]

          merged1<-merge(dt_metadata, merged, by="SampleName")
          merged1<-merged1[complete.cases(merged1),]

          merged1[,SampleName:=NULL]

          setorderv(merged1, category)

          uniq_days<-unique(merged1[[category]])
          qt_days<-length(uniq_days)
          if(qt_days > 25){
            min<-min(merged1[[category]])
            max<-max(merged1[[category]])
            steps <- floor(max/20)
            print(min)
            print(max)
            print(steps)
            groups<-seq(min,max,steps)
            lg<-length(groups)
            groups[lg]<-max

            merged1[,(category):=cut(merged1[[category]], breaks = 25, labels = F)]
            merged1[,(category):=as.integer(groups[merged1[[category]]])]
            # category<-"age in days"
            # merged1[,(category) := NULL]
            merged1[,Abundance:=mean(Abundance), by=c(category, taxon_level)]

            merged1<-unique(merged1)
          }else{
            # colnames(merged1)<-c("group", taxon_level, "Abundance")
            merged1[,Abundance:=mean(Abundance), by=c(category, taxon_level)]
            merged1<-unique(merged1)
          }

          merged1[[taxon_level]]<-factor(merged1[[taxon_level]], levels = rev_ordered)
          chart<-ggplot(merged1, aes_string(x=sprintf("`%s`",category),
                                       y="Abundance", fill=taxon_level)) +
            geom_area(alpha=0.9, size=1)+
            theme_eupath_default(
              legend.title.align=0.4,
              legend.title = element_text(colour="black", size=rel(1), face="bold")
            )+
            scale_fill_manual(values=eupath_pallete, name=taxon_level,
                              labels = c(wrapped_labels),
                              guide = guide_legend(reverse=TRUE, keywidth = 1.7, keyheight = 1.7)
            )+
            labs(x=category, y="Phylogenetic Relative Abundance")+
            coord_cartesian(expand=F)
        }else{
          dt_metadata<-merge(dt_metadata, top_ten, by="SampleName")
          dt_metadata<-subset(dt_metadata, !is.na(get(input$category)))
          chart<-ggplot(dt_metadata, aes_string(x="SampleName", y="Abundance", fill=taxon_level))+
            geom_bar(stat="identity", position="stack", color="black")+
            facet_even(as.formula(sprintf("~ `%s`", category)), ncol=1, scales='free_y')+
            theme_eupath_default(
              legend.title.align=0.4,
              legend.title = element_text(colour="black", size=rel(1), face="bold"),
              axis.text.y=element_blank(),
              axis.ticks.y = element_blank()
            )+
            scale_fill_manual(values=eupath_pallete, name=taxon_level,
                              labels = c(wrapped_labels),
                              guide = guide_legend(reverse=TRUE, keywidth = 1.7, keyheight = 1.7)
            )+
            labs(x="Samples", y="Phylogenetic Relative Abundance")+
            coord_flip(expand=F)
        }
      }

      ggplot_object<<-chart
      ggplot_build_object<<-ggplot_build(chart)

      output$plotWrapper<-renderPlot({
        chart
      })

      if(quantity_samples<MAX_SAMPLES_NO_RESIZE | is_numeric){
        result_to_show<-plotOutput("plotWrapper",
               hover = hoverOpts("overview_hover", delay = 100, delayType = "debounce"),
               click = clickOpts("overview_click"),
               dblclick = dblclickOpts("plot_dblclick"),
               width = paste0(WIDTH,"px"), height = "500px"
             )
      }else{
        result_to_show<-plotOutput("plotWrapper",
           hover = hoverOpts("overview_hover", delay = 100, delayType = "debounce"),
           click = clickOpts("overview_click"), width = paste0(WIDTH,"px"),
           dblclick = dblclickOpts("plot_dblclick"),
           height = quantity_samples*MIN_HEIGHT_AFTER_RESIZE
         )
      }
      shinyjs::hide("chartLoading", anim = TRUE, animType = "slide")
      shinyjs::show("divContent")
    }
    result_to_show
  })


  topAbundance <- function(){}

  output$chartByTopOTU <- renderUI({

    mstudy <- load_microbiome_data()
    category<-input$category
    taxon_level<-input$taxonLevel

    result_to_show<-NULL

    if(!identical(category, "")){

      shinyjs::hide("topTabContent")
      shinyjs::show("topTabLoading")

      quantity_samples <- mstudy$get_sample_count()

      top_ten<-mstudy$get_top_n_by_mean(taxonomy_level = taxon_level, n = NUMBER_TAXA,
                                        add_other = F)

      ordered<-mstudy$otu_table$get_ordered_otu(NUMBER_TAXA)

      top_ten[[taxon_level]]<-factor(top_ten[[taxon_level]], levels=ordered)

      ADDITIONAL<-1

      if(identical(category, NO_METADATA_SELECTED)){
        chart<-ggplot(top_ten, aes_string(x=taxon_level, y="Abundance"))+
          geom_boxplot()+
          theme_eupath_default(
            legend.title = element_text(colour="black", size=rel(1), face="bold"),
            axis.text.x = element_text(angle = 45, hjust = 1)
          )+
          labs(x=taxon_level, y="Relative Abundance")+
          scale_y_continuous(limits = c(0, 1))+
          scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))
      }else{ # end if identical(category, NO_METADATA_SELECTED)
        dt_metadata<-mstudy$get_metadata_as_column(category)
        dt_metadata<-merge(dt_metadata, top_ten, by="SampleName")

        if(is.numeric(dt_metadata[[category]])){
          ADDITIONAL <- 1.5
          dt_metadata<-dt_metadata[!is.na(get(category)),]
          qt_days <- length(unique(dt_metadata[[category]]))

          if(qt_days > 25){
            min<-min(dt_metadata[[category]])
            max<-max(dt_metadata[[category]])
            steps <- floor(max/20)
            groups<-seq(min,max,steps)
            lg<-length(groups)
            groups[lg]<-max

            dt_metadata[,(category):=cut(dt_metadata[[category]], breaks = 25, labels = F)]
            dt_metadata[,(category):=as.integer(groups[dt_metadata[[category]]])]
          # #   # category<-"age in days"
          # #   # merged1[,(category) := NULL]
          #   dt_metadata[,Abundance:=mean(Abundance), by=c(category, taxon_level)]
          # #
          #   dt_metadata<-unique(merged1)
          # # }else{
          # #   # colnames(merged1)<-c("group", taxon_level, "Abundance")
          # #   merged1[,Abundance:=mean(Abundance), by=c(category, taxon_level)]
          # #   merged1<-unique(merged1)
          }
          dt_metadata<-subset(dt_metadata, !is.na(get(category)))
          chart<-ggplot(dt_metadata, aes(x=get(category), y=Abundance, group=get(category)))+
            geom_boxplot()+
            facet_even(as.formula(sprintf("~ `%s`", taxon_level)), ncol=1, scales='free_y')+
            theme_eupath_default(
              legend.title = element_text(colour="black", face="bold"),
              axis.text.x = element_text(angle = 45, hjust = 1)
            )+
            labs(x=category, y="Relative Abundance")

        }else{
          dt_metadata<-subset(dt_metadata, !is.na(get(category)))
          chart<-ggplot(dt_metadata, aes(x=get(taxon_level), y=Abundance, fill=get(category)))+
            geom_boxplot()+
            theme_eupath_default(
              legend.title = element_text(colour="black", face="bold"),
              axis.text.x = element_text(angle = 45, hjust = 1)
            )+
            guides(fill = guide_legend(keywidth = 1.7, keyheight = 1.7))+
            labs(x=taxon_level, y="Relative Abundance",fill=category)+
            scale_y_continuous(limits = c(0, 1))+
            scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))
        }

        # calculating stats
        samples_with_details <- mstudy$sample_table$get_samples_details(category)
        unique_details <- mstudy$sample_table$get_unique_details(category)
        data_frame_table <- data.frame()

        if(length(unique_details)==2){
          # the merged solution is to avoid comparing empty vectors
          # now we only have all samples being compared and if there is no taxa for
          # that particular sample we place 0
          for(i in 1:length(ordered)){
            taxa_data <- subset(dt_metadata, get(taxon_level)==ordered[i], select=c("SampleName","Abundance"))
            merged<-merge(samples_with_details, taxa_data, by="SampleName", all = TRUE)
            set(merged,which(is.na(merged[[4L]])),4L,0)
            result<-wilcox.test(merged$Abundance ~ merged$Value, conf.level = 0.95)
            df<-data.frame(a=ordered[i], "W"=result$statistic, "P-value"=result$p.value)
            data_frame_table<-rbind(data_frame_table, df)
          }
          colnames(data_frame_table)<-c(taxon_level, "W", "P-Value")

          data_frame_table[,3]<-format(data_frame_table[,3], scientific = F)

          sketch <- tags$table(
            tags$thead(
              tags$tr(
                tags$th(style="text-align:center;", rowspan = 2, taxon_level),
                tags$th(style="text-align:center;",colspan = 2, 'Wilcoxon rank sum test')
              ),
              tags$tr(
                tags$th(style="text-align:center;","W"),
                tags$th(style="text-align:center;","P-Value")
              )
            )
          )

          data_frame_table[,taxon_level] <- paste0("<a class='link_table'onclick='goToOtuTab(\"",data_frame_table[,taxon_level],"\")'>",data_frame_table[,taxon_level],"</a>")
          output$by_top_otu_datatable <- DT::renderDataTable(data_frame_table, escape = F, selection = 'none',
                                                             container = sketch, rownames = FALSE)


        }else{
          for(i in 1:length(ordered)){
            taxa_data <- subset(dt_metadata, get(taxon_level)==ordered[i], select=c("SampleName","Abundance"))
            merged<-merge(samples_with_details, taxa_data, by="SampleName", all = TRUE)
            # TODO replace other ways to change 0 to the follow line
            set(merged,which(is.na(merged[[4L]])),4L,0)
            result<-kruskal.test(merged$Value ~ merged$Abundance)
            df<-data.frame("a"=ordered[i], "chi-squared"=unname(result$statistic), "df"=result$parameter, "P-value"=result$p.value)
            data_frame_table<-rbind(data_frame_table, df)
          }
          colnames(data_frame_table)<-c(taxon_level, "chi-squared", "df", "P-Value")
          data_frame_table[,4]<-format(data_frame_table[,4], scientific = F)

          sketch <- tags$table(
            tags$thead(
              tags$tr(
                tags$th(style="text-align:center;", rowspan = 2, taxon_level),
                tags$th(style="text-align:center;",colspan = 3, 'Kruskal-Wallis rank sum test')
              ),
              tags$tr(
                tags$th(style="text-align:center;","chi-squared"),
                tags$th(style="text-align:center;","df"),
                tags$th(style="text-align:center;","P-Value")
              )
            )
          )

          data_frame_table[,taxon_level] <- paste0("<a class='link_table'onclick='goToOtuTab(\"",data_frame_table[,input$taxonLevel],"\")'>",data_frame_table[,taxon_level],"</a>")
          output$by_top_otu_datatable <- DT::renderDataTable(data_frame_table, escape = F, selection = 'none',
                                                             container = sketch, rownames = FALSE)

        }
      } # end else: if(identical(category, NO_METADATA_SELECTED))


      ggplot_by_top_otu_object <<- chart
      ggplot_build_by_top_otu_object <<- ggplot_build(chart)

      output$topOtuPlotWrapper<-renderPlot({
        chart
      })

      result_to_show<-plotOutput("topOtuPlotWrapper", width = "100%", height = paste0(500*ADDITIONAL,"px"),
                                 dblclick = dblclickOpts("plot_dbclick_top_otu"),
                                 hover = hoverOpts("plot_hoverByTopOTU", delay = 60, delayType = "throttle"))


      shinyjs::hide("topTabLoading", anim = TRUE, animType = "slide")
      shinyjs::show("topTabContent")
    }


    result_to_show
  })

  abundanceByOTU <- function(){}
  output$chartByOTU <- renderUI({

    mstudy <- load_microbiome_data()
    category<-input$category
    otu_picked <- input$filterOTU
    isolate(taxon_level<-input$taxonLevel)
    result_to_show<-NULL

    if(is.null(hash_colors)){
      random_color<-sample(eupath_pallete)
      localPallete <- c(random_color[1], "#363636")
    }else{
      if(is.na(match(otu_picked, names(hash_colors)))){
        localPallete <- c("#363636", hash_colors[["Other"]])
      }else{
        localPallete <- c("#363636", hash_colors[[otu_picked]])
      }
    }

    if(!identical(category, "") & !identical(otu_picked, "") & !identical(taxon_level,"")){
      shinyjs::hide("singleOtuContent")
      shinyjs::show("singleOtuLoading")

      quantity_samples <- mstudy$get_sample_count()
      other<-mstudy$otu_table$OTHER_TEXT
      if(identical(category, NO_METADATA_SELECTED)){
        shinyjs::hide("result_tests")
        data_to_show <- mstudy$get_single_otu(taxon_level, otu_picked, T, T)

        data_to_show[[taxon_level]]<-factor(data_to_show[[taxon_level]],
                                            levels=c(other, otu_picked))

        chart<-ggplot(data_to_show, aes_string(x="SampleName", y="Abundance", fill=taxon_level))+
          # geom_bar(position = position_stack(reverse = TRUE), stat = "identity", color="black")+
          geom_bar(stat = "identity", color="black")+
          theme_eupath_default(
            legend.title.align=0.4,
            legend.title = element_text(colour="black", size=rel(1), face="bold")
          )+
          scale_fill_manual(values=localPallete, labels=c(other, otu_picked), name=taxon_level,
                            guide = guide_legend(reverse=T, keywidth = 1.7, keyheight = 1.7))+
          labs(x="Samples", y=paste(otu_picked,"Relative Abundance"))+
          coord_flip(expand=F)

        output$singleOtuBarplotWrapper<-renderPlot({
          chart
        })

        if(quantity_samples<MAX_SAMPLES_NO_RESIZE){
          result_to_show<-plotOutput("singleOtuBarplotWrapper", width = paste0(WIDTH,"px"),
                                     height = "500px",
                                     hover = hoverOpts("singleOtuTooltip", delay = 100, delayType = "debounce")
                                     )
        }else{
          result_to_show<-plotOutput("singleOtuBarplotWrapper", width = paste0(WIDTH,"px"),
                                     height = quantity_samples*MIN_HEIGHT_AFTER_RESIZE,
                                     hover = hoverOpts("singleOtuTooltip", delay = 100, delayType = "debounce")
                                     )
        }
      }else{ # end if(identical(category, NO_METADATA_SELECTED))
        otu_data <- mstudy$get_single_otu(taxon_level, otu_picked, F)
        metadata_as_column <- mstudy$get_metadata_as_column(category)

        col_renamed <- make.names(category)
        colnames(metadata_as_column)<-c("SampleName", col_renamed)

        # to plot we don't show samples with 0 abundance
        merged_to_plot<-merge(otu_data, metadata_as_column, by="SampleName")
        # to calculate the statistics we work with all samples
        merged_to_stats<-merge(metadata_as_column, otu_data, by="SampleName", all=T)

        chart<-ggplot(merged_to_plot, aes_string(x=col_renamed, y="Abundance"))+geom_boxplot()+
          theme_eupath_default()+
          labs(x=stringi::stri_trans_totitle(category), y=paste(otu_picked, "Relative Abundance"))

        output$singleOtuPlotWrapper<-renderPlot({
          chart
        })

        result_to_show<-plotOutput("singleOtuPlotWrapper", width = paste0(WIDTH,"px"),
                                   height = "500px",
                                   hover = hoverOpts("singleOtuTooltip", delay = 60, delayType = "throttle")
        )

        merged_to_stats[is.na(merged_to_stats)] <- 0
        unique_details <- mstudy$sample_table$get_unique_details(category)

        quantity <- length(unique_details)
        abundance_col <- ncol(merged_to_stats)
        category_col <- 2
        if(quantity==2){
          result<-wilcox.test(merged_to_stats[[abundance_col]]~merged_to_stats[[category_col]], conf.int = T, conf.level = 0.95)
          html_formatted<-HTML(sprintf("<ul class=\"shell-body\"> <li>Wilcoxon rank sum test: W = %f, p-value = %.8f</li></ul>", result$statistic, result$p.value))
        }else{
          result<-kruskal.test(merged_to_stats[[category_col]]~merged_to_stats[[abundance_col]])
          html_formatted<-HTML(sprintf("<ul class=\"shell-body\"> <li>Kruskal-Wallis rank sum test: chi-squared = %f, df = %f, p-value = %.8f</li></ul>",
                                       result$statistic, result$parameter, result$p.value))
        }
        shinyjs::show("result_tests")
        output$result_tests <-renderUI({html_formatted})
      } # end else (if(identical(category, NO_METADATA_SELECTED)))

      ggplot_by_otu_object <<- chart
      ggplot_build_by_otu_object <<- ggplot_build(chart)

      shinyjs::hide("singleOtuLoading", anim = TRUE, animType = "slide")
      shinyjs::show("singleOtuContent")
    }else{
      shinyjs::hide("result_tests")
      output$by_otu_datatable <- renderDataTable(NULL)
    }
    result_to_show
  })

  hoversFunctions <- function(){}

  output$overviewTooltip <- renderUI({
    hover <- input$overview_hover
    if (is.null(hover) || is.null(hover$x) || is.null(hover$y) || round(hover$x) <0 || round(hover$y)<0 ) {
      return(NULL)
    }

    if(!identical(input$category, NO_METADATA_SELECTED)){
      if(is.numeric(ggplot_object$data[[input$category]])){
        return(NULL)
      }
    }
    barplot_points(ggplot_object, hover, WIDTH, T, 6, 6, 6)
  })

  output$hoverByOTU <- renderUI({
    hover <- input$singleOtuTooltip
    if (is.null(hover) || is.null(hover$x) || is.null(hover$y) || round(hover$x) <=0 || round(hover$y)<0 ) {
      return(NULL)
    }
    if(identical(input$category, NO_METADATA_SELECTED)){
      barplot_points(ggplot_by_otu_object, hover, WIDTH, T, 20, 20, 6)
    }else{
      smp_details<-ggplot_build_by_otu_object$layout$panel_ranges[[1]]$x.labels # x axis
      if(round(hover$x)>length(smp_details)){
        return(NULL)
      }
      line_data<-ggplot_build_by_otu_object$data[[1]][round(hover$x),]
      boxplot_hover <- get_single_boxplot_hover(hover, line_data$xmax, line_data$x, line_data$ymin, line_data$lower, line_data$middle, line_data$upper, line_data$ymax, WIDTH) # removed text to treat better
      # in the future
      HTML(boxplot_hover)
    }
  })



  output$hoverByTopOTU <- renderUI({
    hover <- input$plot_hoverByTopOTU

    topotus<-ggplot_build_by_top_otu_object$layout$panel_ranges[[1]]$x.labels # x axis
    topotus<-unlist(topotus)

    if (is.null(hover$x) || round(hover$y) < 0 || round(hover$x) < 1 || round(hover$x) > length(topotus))
      return(NULL)

    if(identical(input$category, NO_METADATA_SELECTED)){
      hover_otu <- topotus[round(hover$x)]
      line_data<-ggplot_build_by_top_otu_object$data[[1]][round(hover$x),]
      boxplot_hover <- get_single_boxplot_hover(hover, line_data$xmax, line_data$x, line_data$ymin, line_data$lower, line_data$middle, line_data$upper, line_data$ymax, WIDTH) # removed text to treat better
                                                                  # in the future
      HTML(boxplot_hover)
    }else{
      hover_otu <- topotus[round(hover$x)]
      group_colors <- ggplot_build_by_top_otu_object$plot$scales$scales[[3]]$palette.cache

      group_names <- ggplot_build_by_top_otu_object$plot$scales$scales[[3]]$range$range

      names(group_names)<-group_colors

      interval_value <- ggplot_build_by_top_otu_object$data[[1]][1,"xmax"]-ggplot_build_by_top_otu_object$data[[1]][1,"xmin"]
      line_data<-subset(ggplot_build_by_top_otu_object$data[[1]], hover$x>=xmin & hover$x<=xmax)

      if(nrow(line_data)==1){
        category_hover <- group_names[[line_data$fill]]
        text_hover <- paste0(hover_otu, " [ ", category_hover, " ]")
        boxplot_hover <- get_single_boxplot_hover(hover, line_data$xmax, line_data$x, line_data$ymin, line_data$lower, line_data$middle, line_data$upper, line_data$ymax, WIDTH) # removed text to treat better
                                                                            # in the future
        HTML(boxplot_hover)
      }
    }
  })


  clicks <- function(){}

  observeEvent(input$overview_click, {
    click <- input$overview_click
    if (is.null(click$y))
      return(NULL)

    if(!identical(input$category, NO_METADATA_SELECTED)){
      if(is.numeric(ggplot_object$data[[input$category]])){
        return(NULL)
      }
    }
    mstudy <- load_microbiome_data()
    sample_names <- mstudy$get_sample_names()
    
    sample<-get_sample_from_action(ggplot_object, click)
    
    # # print(length(sample_names)) # 191 samples
    # if(identical(input$category, NO_METADATA_SELECTED)){
    #   # sample <- sample_names[round(click$y)]
    #   sample <- ggplot_build_object$layout$panel_ranges[[1]]$y.labels[round(click$y)]
    # }else{
    #   panel_layout<-ggplot_build_object$layout$panel_layout # panel_layout returns a data.frame
    #   facet_column <- colnames(panel_layout)[4]
    #   # hover$panelvar1 returns the name of the facet where the mouse are hovering
    #   if(is.null(click$panelvar1)){
    #     # panel_hover<-panel_layout[is.na(panel_layout[[facet_column]]),]
    #     panel_hover <- subset(panel_layout, is.na(get(facet_column))) 
    #     
    #   }else{
    #     panel_hover <- subset(panel_layout, get(facet_column)==click$panelvar1) 
    #   }
    #   panel_index <- panel_hover$PANEL
    #   sample<-ggplot_build_object$layout$panel_ranges[[panel_index]]$y.labels[round(click$y)]
    # }
    
    # sample
    raw_data<-mstudy$get_otu_by_sample(sample)[,c("SampleName", input$taxonLevel , "Abundance"),with=F]

    colnames(raw_data)<-c("Sample", colnames(raw_data)[2:length(colnames(raw_data))])

    raw_data$Abundance<-format(raw_data$Abundance, scientific = F)

    raw_data[, (input$taxonLevel) := lapply(.SD, function(x){
      sprintf("<a class='link_table' onclick='goToOtuTab(\"%s\")'>%s</a>", x, x)
      } ), .SDcols=input$taxonLevel]

    output$overviewDatatable <- renderDataTable(raw_data,
                                                  escape = F,
                                                  options = list(
                                                    order = list(list(2, 'desc'),list(1, 'asc'))
                                                  )
    )
  })

  observeEvent(input$plot_dblclick, {
    click <- input$plot_dblclick

    if (is.null(click$y))
      return(NULL)

    mstudy <- load_microbiome_data()
    sample_names <- mstudy$get_sample_names()

    sample<-""
    otu<-NULL

    
    sample<-get_sample_from_action(ggplot_object, click)
    
    if(identical(input$category, NO_METADATA_SELECTED)){
      # sample <- sample_names[round(click$y)]

      click_data<-subset(ggplot_object$data, SampleName==sample & Abundance>0)
      
      layer_data_click<-subset(layer_data(ggplot_object), x==round(click$y))
      
      unique_y<-unique(layer_data_click$y)

      abundances_filtered <- get_abundances_from_plot(unique_y)

      abundances_joined <- join_abundance(abundances_filtered, click_data)

      all_sum <- cumsum(abundances_joined$Abundance)

      index_abundance_click = get_abundance_index(all_sum, click$x)

      if(index_abundance_click == -1)
        otu<-NULL

      otu<-abundances_joined[index_abundance_click,get(input$taxonLevel)]
    }else{
      pnl_layout <- ggplot_build_object$layout$panel_layout
      renamed_category <- make.names(input$category)
      panel_index <- pnl_layout[ pnl_layout[[renamed_category]] == click$panelvar1 , ]$PANEL
      if(length(panel_index) > 0){
        lvls <- ggplot_build_object$layout$panel_ranges[[panel_index]]$y.labels
        # sample <- lvls[round(click$y)]
        if(!is.na(sample)){
          click_data<-subset(ggplot_object$data, SampleName==sample & Abundance>0)

          layer_data_click<-subset(layer_data(ggplot_object), x==round(click$y) & PANEL==panel_index)
          unique_y<-unique(layer_data_click$y)
          unique_y<-unique_y[unique_y>0]
          abundances_filtered <- get_abundances_from_plot(unique_y)

          abundances_joined <- join_abundance(abundances_filtered, click_data)
          all_sum <- cumsum(abundances_joined$Abundance)

          index_abundance_click = get_abundance_index(all_sum, click$x)

          if(index_abundance_click == -1)
            otu<-NULL

          otu<-abundances_joined[index_abundance_click,get(input$taxonLevel)]

          abundances_joined_fill <- abundances_joined

        }else{
          otu<-NULL
        }
      }else{
        otu<-NULL
      }
    }
    if(!is.null(otu)){
      updateSelectizeInput(session, "filterOTU",
                           choices = global_otus,
                           selected = otu,
                           options = list(placeholder = 'Choose a OTU'),
                           server = TRUE)
      updateTabsetPanel(session, "tabs", selected = "byOTU")
    }
  })

  observeEvent(input$plot_dbclick_top_otu, {
    click <- input$plot_dbclick_top_otu

    if (is.null(click$y))
      return(NULL)

    topotus<-ggplot_build_by_top_otu_object$layout$panel_ranges[[1]]$x.labels # x axis
    topotus<-unlist(topotus)

    otu <- topotus[round(click$x)]
    otu <- gsub("\n", " ", otu)

    if(!is.null(otu)){
      updateSelectizeInput(session, "filterOTU",
                           choices = global_otus,
                           selected = otu,
                           options = list(placeholder = 'Choose a OTU'),
                           server = TRUE)
      updateTabsetPanel(session, "tabs", selected = "byOTU")
    }

  })

  downloads<-function(){}
  output$btnDownloadPNG <- downloadHandler(
    filename = "plot.png",
    content = function(file) {
      png(file, width=1200,height=800,units="px")
      if(identical(input$tabs, "bySample")){
        print(ggplot_object)
      }else if(identical(input$tabs, "byOTU")){
        print(ggplot_by_otu_object)
      }else{
        print(ggplot_by_top_otu_object)
      }
      dev.off()
    }
  )

  output$btnDownloadEPS <- downloadHandler(
    filename = "plot.eps",
    content = function(file) {
      setEPS()
      postscript(file, width=16,height=10.67)
      if(identical(input$tabs, "bySample")){
        print(ggplot_object)
      }else if(identical(input$tabs, "byOTU")){
        print(ggplot_by_otu_object)
      }else{
        print(ggplot_by_top_otu_object)
      }
      dev.off()
    }
  )

  output$btnDownloadSVG <- downloadHandler(
    filename = "plot.svg",
    content = function(file) {
      svg(file, width=16,height=10.67)
      if(identical(input$tabs, "bySample")){
        print(ggplot_object)
      }else if(identical(input$tabs, "byOTU")){
        print(ggplot_by_otu_object)
      }else{
        print(ggplot_by_top_otu_object)
      }
      dev.off()
    }
  )

  output$btnDownloadCSV <- downloadHandler(
    filename = "data.csv",
    content = function(file) {
      if(identical(input$tabs, "bySample")){
        write.csv(ggplot_object$data, file)
      }else if(identical(input$tabs, "byOTU")){
        write.csv(ggplot_by_otu_object$data, file)
      }else{
        write.csv(ggplot_by_top_otu_object$data, file)
      }
    }
  )

  shinyjs::hide(id = "loading-content", anim = TRUE, animType = "fade")
  shinyjs::show("app-content")
})
