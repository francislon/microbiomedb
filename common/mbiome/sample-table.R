SampleClass <- R6Class("SampleClass",
      private = list(
        dcast_sample_dt = NULL,
        quantitative_metadata = NULL,
        qualitative_metadata = NULL,
        sample_names = NULL,
        is_eupath = NULL
      ),
      public = list(
        sample_dt = NULL,
        initialize = function(eupath_dt = NULL, biom_dt=NULL, details_dt=NULL) {
          if(!is.null(eupath_dt)){
            #to do remove that after iodice fix the duplicated metadata
            self$sample_dt <- subset(eupath_dt, Property!="source" | (Property=="source" & Value!="Pet" ) )

            private$sample_names <- unique(self$sample_dt$SampleName)
            private$quantitative_metadata <-
              unique(subset(self$sample_dt, Type=="number")[["Property"]])

            private$qualitative_metadata <-
              unique(subset(self$sample_dt, Type!="number")[["Property"]])
            private$is_eupath <- T
          }else if(!is.null(biom_dt)){
            private$sample_names<-biom_dt$SampleName
            private$dcast_sample_dt = biom_dt[,c("SampleName", details_dt$Property), with=F]
            if(is.null(details_dt)){
              private$qualitative_metadata <- details_dt$Property
            }else{
              private$quantitative_metadata <- details_dt[details_dt[,Type=="number"],Property]
              private$dcast_sample_dt[, (private$quantitative_metadata) := lapply(.SD, as.numeric), .SDcols = private$quantitative_metadata]

              private$qualitative_metadata <- details_dt[details_dt[,Type!="number"],Property]
            }
            private$is_eupath <- F
          }
        },
        get_samples_details = function(metadata){
          if(private$is_eupath){
            subset(self$sample_dt, Property==metadata, select = c("SampleName", "Property", "Value"))
          }else{
            result<-private$dcast_sample_dt[,c("SampleName", metadata), with=F]
            result$Property<-metadata
            colnames(result)<-c("SampleName", "Value", "Property")
            result
          }
        },
        get_unique_details = function(metadata){
          sample_details<-self$get_samples_details(metadata)
          unique(sample_details$Value)
        },
        get_sample_count = function(){
          length(private$sample_names)
        },
        get_qualitative_metadata = function(filter = NULL){
          private$qualitative_metadata
        },
        get_quantitative_metadata = function(filter = NULL){
          private$quantitative_metadata
        },
        # TODO: remover filter_samples do objeto e filtrar o dataframe antes de passar para o construtor
        get_metadata_as_column = function(metadata = NULL, filter_samples = NULL) {
          if(is.null(private$dcast_sample_dt)){
            if(is.null(filter_samples)){
              private$dcast_sample_dt <-dcast(data = self$sample_dt,formula = SampleName~Property, value.var = "Value")
            }else{
              private$dcast_sample_dt <-dcast(
                data = subset(self$sample_dt, SampleName != filter_samples),formula = SampleName~Property, value.var = "Value")
            }
            qm<-self$get_quantitative_metadata()
            private$dcast_sample_dt[, (qm):=lapply(.SD, as.numeric), .SDcols=qm]
          }
          if(is.null(metadata)){
            private$dcast_sample_dt
          }else{
            private$dcast_sample_dt[,c("SampleName", metadata), with=F]
          }

        },

        get_sample_names = function(){
          private$sample_names
        },

        get_filtered_categories = function(){
          filtered_categories<-NULL
          if(private$is_eupath){
            setkey(self$sample_dt, Property)
            columns <- unique(self$sample_dt[,Property])
            for(i in 1:length(columns)){
              sample_part <- subset(self$sample_dt, Property==columns[i])
              if(identical(sample_part[1,Type],"number")){
                unique_factors <- as.factor(as.numeric(sample_part$Value))
              }else{
                unique_factors <- as.factor(sample_part$Value)
              }
              if(length(levels(unique_factors)) > 1){
                new_columns <- paste0(columns[i], " (",length(levels(unique_factors)), ")")
                filtered_categories[[new_columns]] <- columns[i]
              }
            }
            setkey(self$sample_dt, NULL)
          }else{
            columns <- colnames(private$dcast_sample_dt)[2:length(private$dcast_sample_dt)]
            for(i in 1:length(columns)){
              if(columns[i] %chin% private$quantitative_metadata){
                unique_factors <- as.factor(as.numeric(private$dcast_sample_dt[[columns[i]]]))
              }else{
                unique_factors <- as.factor(private$dcast_sample_dt[[columns[i]]])
              }
              if(length(levels(unique_factors)) > 1){
                new_columns <- paste0(columns[i], " (",length(levels(unique_factors)), ")")
                filtered_categories[[new_columns]] <- columns[i]
              }
            }
          }
          filtered_categories
        }
      )
)
