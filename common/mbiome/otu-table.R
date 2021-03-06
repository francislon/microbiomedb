OTUClass <- R6Class("OTUClass",
    private = list(
      summarized_otu = NULL,
      summarized_means = NULL,
      last_taxonomy = "",
      use_relative_abundance = T,
      ordered_n_otu = NULL,
      ordered_all_otu = NULL,
      sample_names = NULL
    ),

    public = list(
      otu_dt = NULL,
      OTHER_TEXT = "Other",
      initialize = function(otu_dt, aggregate_by = "Phylum", use_relative_abundance=T) {
        self$otu_dt <- otu_dt
        private$sample_names <- unique(self$otu_dt$Sample)
        self$reshape(aggregate_by, use_relative_abundance)

      },
      reshape = function(taxonomy_level, use_relative_abundance=NULL){

        condA <- is.null(use_relative_abundance)
        condB <- identical(taxonomy_level, private$last_taxonomy)
        if( condA & condB ){
          return(self)
        }else if(!condA){
          condC <- identical(private$use_relative_abundance, use_relative_abundance)
          if(condB & condC){
            return(self)
          }
        }

        if(!condA){
          private$use_relative_abundance <- use_relative_abundance
        }

        private$last_taxonomy <- taxonomy_level

        selected_levels<-get_columns_taxonomy(taxonomy_level)

        if(private$use_relative_abundance){
          private$summarized_otu<-self$otu_dt[, sum(RelativeAbundance), by=c("Sample", selected_levels)]
        }else{
          private$summarized_otu<-self$otu_dt[, sum(AbsoluteAbundance), by=c("Sample", selected_levels)]
        }

        # private$summarized_otu[,  .(selected_levels) := fix_taxonomy_data_table(.SD), .SDcols=selected_levels, by=1:nrow(private$summarized_otu)]
        # the solution above was slow then I replaced for the follow lines just for performance
        length_levels <- length(selected_levels)
        for(i in 2:length(selected_levels)){
          private$summarized_otu[ get(selected_levels[i])=="N/A",selected_levels[i:length_levels]] <-
            paste("unclassified", private$summarized_otu[get(selected_levels[i])=="N/A",get(selected_levels[i-1])])
        }
        # end of workaroud for perfomance
        colnames(private$summarized_otu)<-c("SampleName", selected_levels, "Abundance")

        private$summarized_means<-private$summarized_otu[,list(Abundance=mean(Abundance)), by=taxonomy_level]

        setorderv(private$summarized_means, c("Abundance", taxonomy_level), c(-1,1))

        private$ordered_all_otu <- private$summarized_means[[taxonomy_level]]

        self
      },

      get_sample_names = function(){
        private$sample_names
      },

      get_sample_count = function(){
        length(private$sample_names)
      },
      
      get_sample_as_column_by_otu = function(taxonomy_level = NULL){
        if(!is.null(taxonomy_level) & !identical(taxonomy_level, private$last_taxonomy)){
          self$reshape(taxonomy_level = taxonomy_level)
        }
        selected_levels<-get_columns_taxonomy(private$last_taxonomy)
        formula_dcast <- paste(selected_levels, collapse ='+')
        formula_dcast <- paste0(formula_dcast, "~", "SampleName")
        dcast(data = private$summarized_otu,formula = as.formula(formula_dcast),fun.aggregate = sum, fill=0,value.var = "Abundance")
      },
      
      get_otu_by_sample = function(sample_name, taxonomy_level = NULL){
        if(!is.null(taxonomy_level) & !identical(taxonomy_level, private$last_taxonomy)){
          self$reshape(taxonomy_level = taxonomy_level)
        }
        subset(private$summarized_otu, SampleName==sample_name)
      },

      get_single_otu = function(taxonomy_level=NULL, otu_name, add_other=F, keep_all=F) {
        if(!is.null(taxonomy_level) & !identical(taxonomy_level, private$last_taxonomy)){
          self$reshape(taxonomy_level)
        }
        
        otu <- subset(private$summarized_otu, get(private$last_taxonomy)==otu_name)
        if(keep_all){
          all_samples<-data.table(SampleName=private$sample_names)
          merged_otus<-merge(all_samples, otu, by="SampleName", all=T)
          set(merged_otus,which(is.na(merged_otus[,c("Abundance")])), c(private$last_taxonomy, "Abundance"), list(otu_name, 0))
          otu<-merged_otus
        }
        
        if(add_other){
          grouped_dt<-copy(otu)
          grouped_dt[,Abundance:=1-Abundance]
          grouped_dt[,get_columns_taxonomy(private$last_taxonomy):=self$OTHER_TEXT]
          rbindlist(list(otu,  grouped_dt))
        }else{
          otu
        }
      },

      get_ordered_otu = function(n = NULL){
        if(is.null(n)){
          private$ordered_all_otu
        }else{
          if(n != length(private$ordered_n_otu)){
            top_n <- head(private$summarized_means, n)
            private$ordered_n_otu <- top_n[[private$last_taxonomy]]
          }
          private$ordered_n_otu
        }
      },

      get_otus_by_level = function(taxonomy_level=NULL){
        if(!is.null(taxonomy_level) & !identical(taxonomy_level, private$last_taxonomy)){
          self$reshape(taxonomy_level)$get_ordered_otu()
        }else{
          self$get_ordered_otu()
        }
      },

      get_top_n_by_mean = function(taxonomy_level=NULL, n=10, add_other=T){
        if(!is.null(taxonomy_level) & !identical(taxonomy_level, private$last_taxonomy)){
          self$reshape(taxonomy_level = taxonomy_level)
        }
        
        top_n <- head(private$summarized_means, n)
        top_n$Abundance<-format(top_n$Abundance, scientific = F)
        
        private$ordered_n_otu <- top_n[[private$last_taxonomy]]
        
        filtered_ten<-subset(private$summarized_otu,
                             get(private$last_taxonomy)%chin%top_n[[private$last_taxonomy]] )
        
        if(add_other){
          grouped_dt<-filtered_ten[,list(Abundance=sum(Abundance)),by="SampleName"]
          grouped_dt[,Abundance:=1-Abundance]
          grouped_dt[,get_columns_taxonomy(private$last_taxonomy):="Other"]
          setcolorder(grouped_dt, colnames(filtered_ten))
          filtered_ten <- rbindlist(list(filtered_ten,  grouped_dt))
        }
        filtered_ten
      },

      get_otu_as_column = function(taxonomy_level = NULL, filter_samples = NULL){
        if(!is.null(taxonomy_level) & !identical(taxonomy_level, private$last_taxonomy)){
          tmp<-self$reshape(taxonomy_level)
        }

        if(is.null(filter_samples) || length(filter_samples)==0){
          dcast(data = private$summarized_otu, fun.aggregate = "sum",
                formula = as.formula(paste0("SampleName~",private$last_taxonomy)),
                value.var = "Abundance", fill=0)
        }else{
          print(filter_samples)
          dcast(data = subset(private$summarized_otu, SampleName != filter_samples),
                formula = as.formula(paste0("SampleName~",private$last_taxonomy)),
                fun.aggregate = "sum", value.var = "Abundance", fill=0)
        }

      },
      
      get_summarized_otu = function(taxonomy_level = NULL){
        if(!is.null(taxonomy_level) & !identical(taxonomy_level, private$last_taxonomy)){
          self$reshape(taxonomy_level = taxonomy_level)
        }
        private$summarized_otu
      },
      calculate_relative_abundance = function(){
        private$summarized_otu[,Abundance:=Abundance/sum(Abundance), by="SampleName"]
      }
    )
)
