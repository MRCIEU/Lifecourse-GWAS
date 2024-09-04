#Functions for analysing phenotypes

#set up packages



read_phenotype_data <- function(phecode, input_dir, agebins) {
  # Check if phenotype is present
  filename <- file.path(input_dir, paste0(phecode, ".txt"))
  if(!file.exists(filename)) {
    message(phecode, " not present")
    return(NULL)
  }
  
  # Read in phenotype
  phen <- data.table::fread(filename)
  
  # Check columns are there as expected
  column_names <- c("FID", "IID", "age", "value", "sex")
  if(!all(column_names %in% names(phen))) {
    print(head(phen))
    stop("expect 'FID', 'IID', 'age' and 'value' columns to be present")
  }

  # Keep only required columns
  #phen <- phen %>% dplyr::select(all_of(column_names))

  # Remove duplicates by FID,IID,age
  phen <- subset(phen, !duplicated(paste(FID, IID, age)))
  phen$V1 <- NULL
  
  # Cut into age bins
  # phen$agebin <- cut(phen$age+1, breaks=c(0:19, seq(20, 120, by=5)))
  # levels(phen$agebin) <- levels(phen$agebin) %>% gsub("\\(", "", .) %>% gsub("]", "", .) %>% gsub(",", "-", .)
  phen <- make_agebin(phen, agebins)
  return(phen)
}

make_agebin <- function(phen, agebins) {
  lsbins <- subset(agebins, ! grepl("^[0-9]", category))
  agebins <- subset(agebins, grepl("^[0-9]", category))
  a <- lapply(1:nrow(agebins), \(i) {
    subset(phen, age >= agebins$min[i] & age < agebins$max[i]) %>%
      mutate(agebin = agebins$category[i])
  }) %>% bind_rows()
  b <- lapply(1:nrow(lsbins), \(i) {
    subset(phen, age >= lsbins$min[i] & age < lsbins$max[i]) %>%
      mutate(lsbin = lsbins$category[i])
  }) %>% bind_rows()
  inner_join(a, select(b, FID, IID, age, lsbin))
}

plot_phen_all <- function(phen, phecode) {
  ggplot(phen, aes(x=agebin, y=value)) +
    geom_boxplot() +
    labs(title=phecode)
}

plot_phen_traj <- function(phen, phecode) {
  ggplot(phen, aes(x=agebin, y=value, group=IID)) +
    geom_line(aes(group=IID), alpha=0.2) +
    labs(title=phecode)
}


detect_outliers <- function(phen, phecode) {
  min <- df %>% 
    filter(pheno_id == phecode) %>% {.$min}
  max <- df %>%
    filter(pheno_id == phecode) %>% {.$max}
  
  outliers <- phen %>% 
    filter(!between(phen$value, min, max)) %>% 
    select('value')
  return(outliers)
}

rank_transform <- function(x) {
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}


remove_outliers <- function(phen, phecode) {
  min <- df %>% 
    filter(pheno_id == phecode) %>% {.$min}
  max <- df %>%
    filter(pheno_id == phecode) %>% {.$max}
  
  outlier_rm <- phen %>% 
    filter(between(phen$value, min, max))
  message("Detected ", nrow(phen) - nrow(outlier_rm), " outliers based on phenotype definition file")
  return(outlier_rm)
}


summarise_phen <- function(data) {
  data %>% 
    summarise(
      n=sum(!is.na(value)), 
      m=mean(value, na.rm=T),
      sd=sd(value, na.rm=T),
      q5=quantile(value, 0.05, na.rm=T),
      q10=quantile(value, 0.10, na.rm=T),
      q15=quantile(value, 0.15, na.rm=T),
      q20=quantile(value, 0.20, na.rm=T),
      q25=quantile(value, 0.25, na.rm=T),
      q30=quantile(value, 0.30, na.rm=T),
      q35=quantile(value, 0.35, na.rm=T),
      q40=quantile(value, 0.40, na.rm=T),
      q45=quantile(value, 0.45, na.rm=T),
      q50=quantile(value, 0.50, na.rm=T),
      q55=quantile(value, 0.55, na.rm=T),
      q60=quantile(value, 0.60, na.rm=T),
      q65=quantile(value, 0.65, na.rm=T),
      q70=quantile(value, 0.70, na.rm=T),
      q75=quantile(value, 0.75, na.rm=T),
      q80=quantile(value, 0.80, na.rm=T),
      q85=quantile(value, 0.85, na.rm=T),
      q90=quantile(value, 0.90, na.rm=T),
      q95=quantile(value, 0.95, na.rm=T), 
      male_sex1=sum(as.numeric(pheno_out$sex==1)),
      female_sex2=sum(as.numeric(pheno_out$sex==2)),
      m_age=mean(age, na.rm=T),
      sd_age=sd(age, na.rm=T)
    )
}



if (!require(polspline)) install.packages("polspline", repos="https://cloud.r-project.org")
if (!require(quantreg)) install.packages("quantreg", repos="https://cloud.r-project.org")
if (!require(rms)) install.packages("rms", repos="https://cloud.r-project.org")
#------------------------------------
# Establish the function "phenoplot":
#------------------------------------
phenoplot <- function(Yvbl, Xvbl, Quantiles=c(0.25,0.5,0.75), knots=NA, Nknots=0, jitter_X=0, jitter_Y=0, plotname=""){
  #-------------------------------------------------------
  # Check the function arguments and prepare the datasets:
  #-------------------------------------------------------
  # Put the input variables into a dataframe:
  PP_data <- data.frame(Y=Yvbl,X=Xvbl)
  # Check that one and only one of knots and Nknots is specified:
  if (is.na(knots[1]) & Nknots==0) stop("One of knots and Nknots must be specified")
  if (!is.na(knots[1]) & Nknots!=0) stop("Only one of knots and Nknots must be specified")
  # If knots is specified, calculate the knot values and check that they are appropriate:
  if (!is.na(knots[1])){
    PP_knotvals <- knots
    if ((min(PP_knotvals)>max(PP_data$X)) | (max(PP_knotvals)<min(PP_data$X))) stop("Knots do not overlap the X data range")
    # The model will fail if there are more than two knots above the maximum X or below the minimum X. Trim the knots accordingly:
    keepme <- rep(TRUE,length(PP_knotvals))
    for (kn in 1:length(PP_knotvals)){
      # Drop the knot if there are >=2 knots at lower values which are greater than the maximum X:
      if (sum(PP_knotvals>max(PP_data$X) & PP_knotvals<PP_knotvals[kn])>=2) keepme[kn] <- FALSE
      # Drop the knot if there are >=2 knots at higher values which are less than the minimum X:
      if (sum(PP_knotvals<min(PP_data$X) & PP_knotvals>PP_knotvals[kn])>=2) keepme[kn] <- FALSE
    }
    if (sum(keepme)<length(keepme)) warning("WARNING: Some knots removed outside observed range of X")
    PP_knotvals <- PP_knotvals[keepme]
    rm(keepme)
    if (length(PP_knotvals)<3 | length(PP_knotvals)>7) stop("knots must contain 3-7 values")  
  }
  # If Nknots is specified, calculate the knot values and check that they are appropriate:
  if (Nknots>0){
    if (Nknots<3 | Nknots>7) stop("Nknots must be zero, or between 3 and 7")
    if (Nknots==3) PP_knotvals <- quantile(PP_data$X,probs=c(0.1,0.5,0.9),type=6)
    if (Nknots==4) PP_knotvals <- quantile(PP_data$X,probs=c(0.05,0.35,0.65,0.95),type=6)
    if (Nknots==5) PP_knotvals <- quantile(PP_data$X,probs=c(0.05,0.275,0.5,0.725,0.95),type=6)
    if (Nknots==6) PP_knotvals <- quantile(PP_data$X,probs=c(0.05,0.23,0.41,0.59,0.77,0.95),type=6)
    if (Nknots==7) PP_knotvals <- quantile(PP_data$X,probs=c(0.025,0.183,0.342,0.5,0.658,0.817,0.975),type=6)
  }
  # Check whether the splines are going to work:
  test_cspl <- try(rms::rcs(PP_data$X,PP_knotvals),silent=TRUE)
  if (class(test_cspl)=="try-error"){
    warning("Cubic splines fail with the specified knots. This may be because there is insufficient variation in age.")
    warning("Fitted values displayed on graph will be simple percentiles of the phenotype, not functions of age")
  }
  if (class(test_cspl)!="try-error"){
    message("Using the following knot placements:")
    print(PP_knotvals)
  }
  #---------------------------------
  # Make fitted values for plotting:
  #---------------------------------
  # Make a data frame containing 1000 equally spaced values of X for the plotting of fitted values:
  # (This is so that plotted fitted values don't reflect individual X data)
  PP_plotdata <- data.frame(X=seq(from=min(PP_data$X),to=max(PP_data$X),length.out=1000))
  for(qn in 1:length(Quantiles)){
    q <- Quantiles[qn]
    if (class(test_cspl)=="try-error"){
      PP_plotdata$YQfit <- quantile(PP_data$Y,probs=q,type=6)
    }
    if (class(test_cspl)!="try-error"){
      # Create the spline values each time, since they must have the same names whether derived from full data or plot data:
      PP_cspl_X <- rms::rcs(PP_data$X,PP_knotvals)
      # Model the quantile using quantile regression:
      qmodel <- quantreg::rq(Y~PP_cspl_X, data=PP_data, tau=q, model=TRUE)
      # Calculate fitted values in plotdata:
      PP_cspl_X <- rms::rcs(PP_plotdata$X,PP_knotvals)
      PP_plotdata$YQfit <- predict(qmodel,PP_cspl_X)
    }  
    names(PP_plotdata)[names(PP_plotdata)=="YQfit"]<-paste0("q",qn)  
  }
  #---------------
  # Make the plot:
  #---------------
  # Set it up to save the graph as a pdf, or show it in R:
  if (plotname!=""){
    print("Saving plot in the working directory")
    pdf(paste0(plotname,".pdf"), width=8, height=5,bg="white", colormodel="cmyk")  
  }
  if (plotname=="") print("Showing plot in R")
  # Make a jittered version of X (+/- a uniformly distributed value between 0 and jitter_X):
  PP_data$Xj <- PP_data$X+runif(n=dim(PP_data)[1],min=(-1)*jitter_X,max=jitter_X)
  # Make a jittered version of Y (+/- a uniformly distributed value between 0 and jitter_Y*IQR(Y)):
  if (class(test_cspl)=="try-error"){
    YPC25 <- rep(quantile(PP_data$Y,probs=0.25,type=6),length(PP_data$Y))
    YPC75 <- rep(quantile(PP_data$Y,probs=0.75,type=6),length(PP_data$Y))
  }
  if (class(test_cspl)!="try-error"){
    PP_cspl_X <- rms::rcs(PP_data$X,PP_knotvals)
    qmodel <- quantreg::rq(Y~PP_cspl_X, data=PP_data, tau=0.25, model=TRUE)
    YPC25 <- predict(qmodel,PP_cspl_X)
    qmodel <- quantreg::rq(Y~PP_cspl_X, data=PP_data, tau=0.75, model=TRUE)
    YPC75 <- predict(qmodel,PP_cspl_X)
    rm(qmodel)
  }
  randoms <- runif(n=dim(PP_data)[1],min=-1,max=1)
  PP_data$Yj <- PP_data$Y+randoms*jitter_Y*(YPC75-YPC25)
  rm(list=c("randoms", "YPC25", "YPC75"))
  # Make a scatter plot of the jittered original data:
  xlab <- "Age (upward ticks are knots, vertical lines are lifestage boundaries)"
  ylab <- "Phenotype value"
  plot(PP_data$Xj,PP_data$Yj,col="grey",pch=1,cex=0.25,lwd=0.5,xlab=xlab,ylab=ylab,yaxt="n",xaxt="n")
  title(paste0("(fitted quantiles at: ",paste(Quantiles, collapse = ", "),")"),font.main=1,cex.main=0.75,line=0.25,adj=0)
  axis(side=1,tck=-0.01)
  axis(side=2,tck=-0.01)
  # Plot the "life stage" boundaries:
  abline(v=c(3,8,12,18,40,65),lwd=0.25,lty=2,col="grey")
  # Plot the knots used in the fitted values, as upward tick marks:
  if (class(test_cspl)!="try-error"){
    axis(side=1,at=PP_knotvals, tck=0.02, lwd.ticks=2, col.ticks="red", labels=FALSE)
  }
  # Add the fitted quantiles:
  for(qn in 1:length(Quantiles)){
    q <- Quantiles[qn]
    PlotQ <- PP_plotdata[,paste0("q",qn)]
    points(PP_plotdata$X,PlotQ,col="red",lwd=0.25,type="l", lty="solid")
    rm(list=c("q","PlotQ"))
  }
  # Finalise the pdf, if appropriate:
  if (plotname!="") dev.off()
}
