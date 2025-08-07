#Functions for analysing phenotypes

#set up packages
nthreads <- as.numeric(Sys.getenv("env_threads"))

# Added safe_hist function
safe_hist <- function(x, main="Histogram") {
  if(length(x) > 0 && sum(!is.na(x)) > 0 && length(unique(round(x))) > 1) {
    hist(x, breaks=length(unique(round(x))), main=main)
  } else {
    message("Skipping histogram: insufficient data for ", main)
  }
}


read_phenotype_data <- function(phecode, input_dir, agebins, covlist=NULL) {
  # Check if phenotype is present
  filename <- file.path(input_dir, paste0(phecode, ".txt"))
  if(!file.exists(filename)) {
    message(phecode, " not present")
    return(NULL)
  }
  
  # Read in phenotype
  phen <- data.table::fread(filename, nThread = nthreads)
  
  # Check columns are there as expected
  column_names <- c("FID", "IID", "age", "value")
  if(!all(column_names %in% names(phen))) {
    print(head(phen))
    stop("expect 'FID', 'IID', 'age', 'value' columns to be present")
  }

  # Remove duplicates by FID,IID,age
  phen <- subset(phen, !duplicated(paste(FID, IID, age)))
  
  phen$age <- as.numeric(phen$age)
  phen$value <- as.numeric(phen$value)
  phen$FID <- as.character(phen$FID)
  phen$IID <- as.character(phen$IID)
  

  # Check expected age-varying covariates are there
  if(!is.null(covlist)) {
    missing_covs <- covlist[!covlist %in% names(phen)]
    if(length(missing_covs) > 0) {
      stop(paste(missing_covs, collapse=","), " - these age specificshould be present for this phenotype")
    }
  }

  # Keep only required columns
  phen <- subset(phen, select=c("FID", "IID", "age", "value", covlist))

  # Cut into age bins
  phen <- make_agebin(phen, agebins)

  return(phen)
}


read_gen_covs <- function(file, npcs) {
  dat <- fread(file, nThread = nthreads)
  if(!all(c("FID", "IID", paste0("PC", 1:npcs)) %in% names(dat))) {
    stop("expected FID, IID, PC1, ..., PC", npcs, " in ", file)
  }
  dat <- subset(dat, select=c("FID", "IID", paste0("PC", 1:npcs)))
  dat <- subset(dat, !duplicated(paste(FID, IID)))
  dat$FID <- as.character(dat$FID)
  dat$IID <- as.character(dat$IID)
  return(dat)
}


read_covariate_data <- function(fn, covariate_list=c("sex", "yob")) {
  dat <- data.table::fread(fn, nThread = nthreads) %>% as_tibble()
  column_names <- c("FID", "IID", covariate_list)
  if(!all(column_names %in% names(dat))) {
    print(head(dat))
    stop("expected FID, IID, ", paste(covariate_list, collapse=", "), " in static_covariates.txt file")
  }
  dat <- subset(dat, select=c("FID", "IID", covariate_list))
  dat <- subset(dat, !duplicated(paste(FID, IID)))
  if("yob" %in% covariate_list) {
  dat$deob <- cut(dat$yob/10, breaks=seq(180,203, by=1))
  }
  if("sex" %in% covariate_list) {
    s <- unique(dat$sex[!is.na(dat$sex)])
    if(!setequal(c(1,2), s)) {
      stop("sex values must be 1=male and 2=female. Found: ", s)
    }
  }

  dat$FID <- as.character(dat$FID)
  dat$IID <- as.character(dat$IID)
  return(dat)
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
      male_sex1=sum(as.numeric(data$sex==1)),
      female_sex2=sum(as.numeric(data$sex==2)),
      m_age=mean(age, na.rm=T),
      sd_age=sd(age, na.rm=T),
      m_age2=mean(age^2,na.rm=T),
      sd_age2=sd(age^2, na.rm=T),
      m_age3=mean(age^3,na.rm=T),
      sd_age3=sd(age^3, na.rm=T)
    )
}


summarise_phen_msd <- function(data) {
  data %>% 
    summarise(
      m_t=mean(value, na.rm=T),
      sd_t=sd(value, na.rm=T),
    )
}





if (!require(polspline)) install.packages("polspline", repos="https://cloud.r-project.org")
if (!require(quantreg)) install.packages("quantreg", repos="https://cloud.r-project.org")
if (!require(rms)) install.packages("rms", repos="https://cloud.r-project.org")
#------------------------------------
# Establish the function "phenoplot":
#------------------------------------
phenoplot <- function(Yvbl, Xvbl, Quantiles=c(0.25,0.5,0.75), knots=NA, Nknots=0, jitter_X=0, jitter_Y=0, title="phenotype plot", plotname=""){
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
  plot(PP_data$Xj,PP_data$Yj,col="grey",pch=1,cex=0.25,lwd=0.5,xlab=xlab,ylab=ylab,yaxt="n",xaxt="n",main=title)
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


###############################
##Function to convert the phenotype data into what we need for the GWAS
#includes:
#-importing the phenotype data
#-categorising age groups
#-removing outliers
#-data transformation
#-dividing the data into the approprate age groups
#-checking sample size for each group
#-saving the file for the GWAS for each group that passes the threshold


##For most phenotypes outlier removal and transformation is done based on the phenotypes list csv file
##Pulse pressure and bmiz are generated from other data provided and the outlier removal values and transformations are set in this code. 


organise_phenotype <- function(phecode, phenotypes, df, gen_covs, cov_list, covdat, agebins, pl=TRUE) {
  
  ##read in the data and any age varying covariates (these should be listed in the phenotype list and be in the phenotypes file)
  
  if(phecode != "pp" & phecode != "bmiz"){
    type <- filter(df, pheno_id == phecode)$var_type
    str(filter(df, pheno_id == phecode))
    cs <- list()
    vary_covs <- unlist(subset(df, pheno_id == phecode)$covs %>% strsplit(., ":"))
    phen <- read_phenotype_data(phecode, Sys.getenv("phenotype_input_dir"), agebins, vary_covs)
  }
  
  ##read in data for dpb and sbp and generate pulse pressure
  #exclude outliers based on dbp, sbp or pp
  
  if(phecode == "pp") {
    type <- "cont"
    cs <- list()
    phecode1 <- "sbp"
    phen1 <- read_phenotype_data(phecode1, Sys.getenv("phenotype_input_dir"), agebins, "bp_med")
    if(is.null(phen1)) {
      return(NULL)
    }
    
    phen1$value <- phen1$value + (15 * phen1$bp_med)
    phen1 <- remove_outliers(phen1, phecode1)
    phen1 <- rename(phen1, sbp = value)
    
    phecode2 <- "dbp"
    phen2 <- read_phenotype_data(phecode2, Sys.getenv("phenotype_input_dir"), agebins, "bp_med")
    if(is.null(phen2)) {
      return(NULL)
    }
    
    phen2$value <- phen2$value + (10 * phen2$bp_med)
    phen2 <- remove_outliers(phen2, phecode2)
    phen2 <- rename(phen2, dbp = value)
   
    #join together the dbp and sbp data - these should only match if all of the other variables match (i.e. ID value, age and bp_med need to match for the variables to be merged)
    phen <- inner_join(phen1, phen2)
    phen$value <- phen$sbp - phen$dbp
    phen <- phen %>%
      select(!c(sbp, dbp))
    
    outliers <-  phen %>% 
      filter(!between(phen$value, 10, 180)) %>% 
      select('value')
    
    phen <- phen %>%
      filter(between(phen$value, 10, 180))
    
  }
  
  ##read in bmi data and generate bmiz for children - only keeps data for <= 18, adults are excluded for bmiz
  #outliers for bmiz are excluded based on having a bmiz <-5 or >5
  
  if(phecode=="bmiz"){
    type <- "cont"
    cs <- list()
    phen <- read_phenotype_data("bmi", Sys.getenv("phenotype_input_dir"), agebins)
    if(is.null(phen)) {
      return(NULL)
    }
    
    phen <- filter(phen, age <= 18)
    if(nrow(phen) == 0) {
      return(NULL)
    }
    phen$Month <- round(phen$age*12)
    
    filename <- here("resources", "phenotypes", "bmi-z-who-2007.csv")
    z_dat <- data.table::fread(filename, nThread = nthreads)

    phen <- inner_join(phen, covdat, by = join_by(FID, IID))
    
    phen <- left_join(phen,z_dat, by=join_by(Month))
    phen$bmiz_girls <- (((phen$value/phen$MG)^phen$LG) - 1)/(phen$SG*phen$LG)
    phen$bmiz_boys <- (((phen$value/phen$MB)^phen$LB) - 1)/(phen$SB*phen$LB)
    phen <- phen %>% 
      mutate(bmiz = case_when(sex==1 ~ bmiz_boys, sex==2 ~ bmiz_girls))
    
    phen <- phen %>% 
      select(FID, IID, age, bmiz, agebin, lsbin) %>%
      rename(value = bmiz)
    
    
    outliers <- phen %>% 
      filter(!between(phen$value, -5, 5)) %>% 
      select('value')
    
    phen <- phen %>% 
      filter(between(phen$value, -5, 5))
    
  }
  
  ##return null and move on to the next trait if this one is missing
  
  if(is.null(phen)) {
    return(NULL)
  }
  
  # Only one genetic observation per person
  setequal(unique(gen_covs$IID), gen_covs$IID)
  
  # Merge genetic covariates with phenotypic data
  dat <- phen %>%
    left_join(gen_covs, by = c("FID", "IID")) %>%
    na.exclude
  
  dat <- left_join(dat, covdat, by = c("FID", "IID")) %>%
    na.exclude

  # Age and ancestry distribution
  if(pl) {
    safe_hist(dat$age, paste0(phecode, " age distribution"))
  }
  cs$age_summary <- summary(dat$age)
  cs$age_quantile <- dat$age %>% quantile(probs=seq(0, 1, 0.01))
 
  if('yob' %in% cov_list) {
    cs$deob_summary <- table(dat$deob)
    phecode
    cs$deob_summary

    if(pl) {
      safe_hist(dat$yob, paste0(phecode, " year of birth distribution"))
    }
  }
  
  cs$sex_table <- table(dat$sex)
  phecode
  cs$sex_table

  if('yob' %in% cov_list) {
  cs$sex_yob <- dat %>% 
    group_by(sex) %>% 
    summarise(
        yob_mean = mean(yob, na.rm=T), 
        yob_sd = sd(yob, na.rm=T)
    )
  phecode
  cs$sex_yob
  }


  #remove outliers and plot the remaining data
  #A jitter is applied to the plot so no actual observations are plotted

  #remove outliers based on values in the phenotypes description for every trait except pp and bmiz
  #outliers for pp and bmiz are removed at the point the data is imported. 
  
  if(phecode != "pp" & phecode != "bmiz"){
  outliers <- detect_outliers(dat, phecode)
  analysis_data <- remove_outliers(dat, phecode)
  } else {
    analysis_data <- dat
  }
  
  
  ##This will break the code if all values are removed as outliers - this breaks rather than moving on to the next trait as it should only happen if there is an error
  #e.g. incorrectly defined outliers, phenotype values entered in the wrong units etc. 
  
  if(length(analysis_data$value) == 0){
    message(paste(phecode,"All observations removed as outliers, skipping."))
    return(NULL)
  }

  #generate the plot 
  age_range = max(analysis_data$age) - min(analysis_data$age)
  y_sd = sd(analysis_data$value)

  phenoplot(analysis_data$value, analysis_data$age, Quantiles=c(0.05,0.25,0.5,0.75,0.95), knots=NA, Nknots=5, jitter_X=0.5, jitter_Y=0.1, plotname=ifelse(pl, "", phecode), title=phecode)
  
  
  # medication adjustment for ldl, sbp and dbp
  
  if(phecode=="ldl"){
    analysis_data$value <- analysis_data$value + (0.40 * analysis_data$value * analysis_data$cholesterol_med)
  }
  
  if(phecode=="sbp"){
    analysis_data$value <- analysis_data$value + (15 * analysis_data$bp_med)
  }
  
  if(phecode=="dbp"){
    analysis_data$value <- analysis_data$value + (10 * analysis_data$bp_med)
  }
  
  
  ##state the transformation that will be applied
  #no transformation is applied to bmiz as its already a standardised variable
  
  if(phecode != "pp" & phecode != "bmiz"){
  if (filter(df, pheno_id == phecode)$transformation=="log") {
    print("log transformed")
  } else if(filter(df, pheno_id == phecode)$transformation=="rank") {
    print("rank transformed")
  } else if(filter(df, pheno_id == phecode)$transformation=="none") {
    print("no transformation")
  }
  }
  
  if(phecode != "pp" & phecode != "bmiz"){
  if (filter(df, pheno_id == phecode)$standardisation=="yes") {
    print("standardised")
  } else if(filter(df, pheno_id == phecode)$transformation=="no") {
    print("no standardisation")
  } 
  }
  if(phecode == "pp") {
    print("standardised")
  } 
  

  ####By age group (1 year <20, 5 years >= 20)
  ############################################
  
  # Generate counts of sample size in each age group by sex and combined
  
  cats_sexspec <- analysis_data %>% 
    group_by(agebin, sex) %>%
    filter(!duplicated(paste(FID, IID))) %>%
    summarize(n())

  colnames(cats_sexspec) <- c("agebin", "sex", "n")
  
  cats <- analysis_data %>% 
    group_by(agebin) %>%
    filter(!duplicated(paste(FID, IID))) %>% 
    summarize(n())
   
  colnames(cats) <- c("agebin", "n")

  cats <- bind_rows(cats_sexspec, cats)
  cs$categories <- cats 

  #filter so know which age group have a large enough sample and make the sex combined sex=3
  cats <- filter(cats, n >= as.numeric(Sys.getenv("env_minumum_strat_n"))) %>%
    replace_na(list(sex=3))
  
  summary_output <- list()
  
  ###this section operates on each age/sex (including sex combined) group individually
  if(length(cats$n>=1)){
    for(k in 1:length(cats$n)){
    
      age_group <- cats$agebin[k]
      sex_group <- cats$sex[k]
  
      #select the individual observations from the phenotypes file that are the desired age/sex
      #check and remove one of any duplicated individuals in that group
      if(sex_group < 3) {
        pheno_out <- filter(analysis_data, (agebin == cats$agebin[k])) %>% 
          filter(sex == cats$sex[k]) %>%
          filter(!duplicated(paste(FID, IID))) %>% 
          select(FID, IID, value, age, sex)
      } else if(sex_group == 3) {
        pheno_out <- filter(analysis_data, (agebin == cats$agebin[k])) %>% 
          filter(!duplicated(paste(FID, IID))) %>% 
          select(FID, IID, value, age, sex)
      }

      
      #summarise the pretransformed variable 
      sumstats <- summarise_phen(pheno_out)
  
      ## apply transformation and standardisation if required
      #applied to all traits except pp and bmiz based on phenotypes description file
      if(phecode != "pp" & phecode != "bmiz"){
      if (filter(df, pheno_id == phecode)$transformation=="log") {
        pheno_out$value <- log(pheno_out$value)
      } else if(filter(df, pheno_id == phecode)$transformation=="rank") {
        pheno_out$value <- rank(pheno_out$value, ties.method = "average")
      } else if(filter(df, pheno_id == phecode)$transformation=="none") {
        pheno_out$value <- pheno_out$value
      } 
      }

      if(phecode != "pp" & phecode != "bmiz"){
      if (filter(df, pheno_id == phecode)$standardisation=="yes") {
        pheno_out$value <- (pheno_out$value - mean(pheno_out$value, na.rm = T))/sd(pheno_out$value,  na.rm = T)
      } else if(filter(df, pheno_id == phecode)$standardisation=="no") {
        pheno_out$value <- pheno_out$value
      } 
      }
      
      #pp is standardised bit not transformed
      if(phecode == "pp") {
        pheno_out$value <- pheno_out$value/sd(pheno_out$value)
      } 
    
      #define the sex groups for saving 
      if(sex_group == 1){
        sex_out = "m"
      } else if (sex_group == 2){
        sex_out = "f"
      } else if (sex_group == 3){
        sex_out = "both"
      }
    
      #save the (transformed) phenotype values for the GWAS
      write.table(subset(pheno_out, select=c("FID", "IID", "value")), file=file.path(Sys.getenv("phenotype_processed_dir"), paste0(phecode, "_", age_group, "_", sex_out, ".phen")), row=FALSE, col=FALSE, qu=FALSE)

      #define the covariates that are saved as adjustments for the GWAS
      if(phecode != "pp" & phecode != "bmiz"){
      cov_ids <- c(unlist(subset(df, pheno_id == phecode)$covs %>% strsplit(., ":")), cov_list)
      }
      if(phecode == "pp"){
        cov_ids <- c("bp_med", cov_list)
      }
      if(phecode == "bmiz"){
        cov_ids <- cov_list
      }
      
      #select these observations from the original data and remove duplicates
      covs <- analysis_data %>% 
        filter(agebin == age_group) %>%
        select(all_of(c(names(gen_covs), cov_ids))) %>% 
        filter(IID %in% pheno_out$IID) %>% 
        filter(!duplicated(paste(FID, IID)))
    
      write.table(covs, file=file.path(Sys.getenv("phenotype_processed_dir"), paste0(phecode, "_", age_group, "_", sex_out, ".covs")), row=FALSE, col=FALSE, qu=FALSE)
  
      sumstats_t <- summarise_phen_msd(pheno_out)
      summary_output[[k]] <- cbind(age_group, sex_group, length(outliers$value), sumstats, sumstats_t)
    }
  
    cs$sums <- bind_rows(summary_output)
    print(cs$sums)
  }
   
  
  ####By lifestage category
  ############################################
  #This section repeats the section above but for each lifestage category rather than each age group
  
  
  #generate counts of sample size in each lifestage group by sex and combined

  catsls_ss <- analysis_data %>% 
    group_by(lsbin, sex) %>%
    filter(!duplicated(paste(FID, IID))) %>%
    summarize(n())
  colnames(catsls_ss) <- c("lsbin", "sex", "n") 
  
  catsls <- analysis_data %>% 
    group_by(lsbin) %>%
    filter(!duplicated(paste(FID, IID))) %>%
    summarize(n())
  colnames(catsls) <- c("lsbin", "n")

  catsls <- bind_rows(catsls_ss, catsls)
  cs$categories_ls <- catsls
  
  #keep the categories where the sample size is greater than the minimum
  catsls <- filter(catsls, n >= as.numeric(Sys.getenv("env_minumum_strat_n"))) %>%
   replace_na(list(sex=3))
  
  summary_output <- list()

  
  ##From here the code is applied to keep age bin and sex combination
  
  if(length(catsls$n>=1)){
    for(k in 1:length(catsls$n)){
      
      age_group <- catsls$lsbin[k]    
      sex_group <- catsls$sex[k] 

      #keep the data from the main dataset that is part of this group and ensure that no individuals are repeated. 
      if(sex_group < 3) {
        pheno_out <- filter(analysis_data, (lsbin == catsls$lsbin[k])) %>% 
        filter(sex == catsls$sex[k]) %>%
        filter(!duplicated(paste(FID, IID))) %>% 
        select(FID, IID, value, age,sex)
      } else if(sex_group == 3) {
        pheno_out <- filter(analysis_data, (lsbin == catsls$lsbin[k])) %>% 
        filter(!duplicated(paste(FID, IID))) %>% 
        select(FID, IID, value, age, sex)
      }
      
      
      #summarise the pretransformed variable 
      sumstats <- summarise_phen(pheno_out)
      
      
      #apply the required transformation to the data based on the phenotypes info file 
      #pp is standardised below, nothing is done to bmiz as it is already standardised
      if(phecode != "pp" & phecode != "bmiz"){
      if (filter(df, pheno_id == phecode)$transformation=="log") {
        pheno_out$value <- log(pheno_out$value)
      } else if(filter(df, pheno_id == phecode)$transformation=="rank") {
        pheno_out$value <- rank(pheno_out$value, ties.method = "average")
      } else if(filter(df, pheno_id == phecode)$transformation=="none") {
        pheno_out$value <- pheno_out$value
      }
      }
          
      ## apply standardisation if required for traits other than bmiz and pp
      if(phecode != "pp" & phecode != "bmiz"){
      if (filter(df, pheno_id == phecode)$standardisation=="yes") {
        pheno_out$value <- (pheno_out$value - mean(pheno_out$value, na.rm = T))/sd(pheno_out$value,  na.rm = T)
      } else if(filter(df, pheno_id == phecode)$standardisation=="no") {
        pheno_out$value <- pheno_out$value
      }
      }
      #for pp the data is standardised but not transformed
      if(phecode == "pp"){
        pheno_out$value <- pheno_out$value/sd(pheno_out$value)
      }

      #define the sex based on the code 
      if(sex_group == 1){
        sex_out = "m"
      } else if (sex_group == 2){
        sex_out = "f"
      } else if (sex_group == 3){
        sex_out = "both"
      }
      
      ##save the phenotype data for the GWAS
      write.table(subset(pheno_out, select=c("FID", "IID", "value")), file=file.path(Sys.getenv("phenotype_processed_dir"), paste0(phecode, "_", age_group, "_", sex_out, ".phen")), row=FALSE, col=FALSE, qu=FALSE)

      #define the covariates that are saved as adjustments for the GWAS
      if(phecode != "pp" & phecode != "bmiz"){
        cov_ids <- c(unlist(subset(df, pheno_id == phecode)$covs %>% strsplit(., ":")), cov_list)
      }
      if(phecode == "pp"){
        cov_ids <- c("bp_med", cov_list)
      }
      if(phecode == "bmiz"){
        cov_ids <- cov_list
      }
      
      
      covs <- analysis_data %>% 
        filter(lsbin == age_group) %>%
        select(all_of(c(names(gen_covs), cov_ids, "age"))) %>% 
        filter(IID %in% pheno_out$IID) %>% 
        filter(!duplicated(paste(FID, IID)))
      
      #save the covariates for the GWAS
      write.table(covs, file=file.path(Sys.getenv("phenotype_processed_dir"), paste0(phecode, "_", age_group, "_", sex_out, ".covs")), row=FALSE, col=FALSE, qu=FALSE)
    
      sumstats_t <- summarise_phen_msd(pheno_out)
      summary_output[[k]] <- cbind(age_group, sex_group, length(outliers$value), sumstats, sumstats_t)
    }

    cs$sums_ls <- bind_rows(summary_output)
  }

  saveRDS(cs, file=file.path(Sys.getenv("results_dir"), "02", paste0(phecode,"_summary.rds")))
}


