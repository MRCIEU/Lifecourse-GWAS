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
  column_names <- c("FID", "IID", "age", "value")
  if(!all(column_names %in% names(phen))) {
    print(head(phen))
    stop("expect 'FID', 'IID', 'age' and 'value' columns to be present")
  }

  # Keep only required columns
  phen <- phen %>% dplyr::select(all_of(column_names))

  # Remove duplicates by FID,IID,age
  phen <- subset(phen, !duplicated(paste(FID, IID, age)))
  
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


summarise_phen <- function(data, age_group, anc_group,male,female) {
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
      female_sex2=sum(as.numeric(pheno_out$sex==2))
    )
}

mbw <- function(MBW_Y, MBW_X, knots=7, lpc=0.005, upc=0.995, multiple=1, plotme="none", comparator=0, jitter_X=0, jitter_Y=0){
  print(paste0("Winsorising at ",multiple," times the interval from the median to percentiles ",lpc," and ",upc))
  # Put the input variables into a dataframe:
  MBW_data <- data.frame(Y=MBW_Y,X=MBW_X)
  # knots must be either 0, or between 3 and 7 inclusive:
  if (knots!=0&(knots<3|knots>7)) stop("Number of knots must be 0,3,4,5,6 or 7") 
  # If knots=0, model simple percentiles of Y with no effect of age:
  if (knots==0){
    Ypc <- quantile(MBW_data$Y,probs=c(lpc,0.25,0.5,0.75,upc),type=6)
    print(Ypc)
    MBW_data$fit_lpc <- Ypc[1]
    MBW_data$fit_25 <- Ypc[2]
    MBW_data$fit_50 <- Ypc[3]
    MBW_data$fit_75 <- Ypc[4]
    MBW_data$fit_upc <- Ypc[5]
  }
  # If knots is 3-7, model Y against cubic splines of X:
  if (knots>=3&knots<=7){
    # Calculate Harrell's recommended percentiles of the exposure for the chosen number of knots:
    # (percentiles are type 6 to match Stata)
    if (knots==3) MBW_knotvals <- quantile(MBW_data$X,probs=c(0.1,0.5,0.9),type=6)
    if (knots==4) MBW_knotvals <- quantile(MBW_data$X,probs=c(0.05,0.35,0.65,0.95),type=6)
    if (knots==5) MBW_knotvals <- quantile(MBW_data$X,probs=c(0.05,0.275,0.5,0.725,0.95),type=6)
    if (knots==6) MBW_knotvals <- quantile(MBW_data$X,probs=c(0.05,0.23,0.41,0.59,0.77,0.95),type=6)
    if (knots==7) MBW_knotvals <- quantile(MBW_data$X,probs=c(0.025,0.183,0.342,0.5,0.658,0.817,0.975),type=6)
    print(MBW_knotvals)
    # Calculate the restricted cubic spline values from these knots:
    MBW_cspl_X <- rcs(MBW_data$X,MBW_knotvals)
    # Model percentiles of MBW_Y as a cubic spline function of MBW_X and get fitted values:
    model_lpc <- rq(Y~MBW_cspl_X, data=MBW_data, tau=lpc, model=TRUE)
    MBW_data$fit_lpc <- predict(model_lpc)
    model_25 <- rq(Y~MBW_cspl_X, data=MBW_data, tau=0.25, model=TRUE)
    MBW_data$fit_25 <- predict(model_25)
    model_50 <- rq(Y~MBW_cspl_X, data=MBW_data, tau=0.5, model=TRUE)
    MBW_data$fit_50 <- predict(model_50)
    model_75 <- rq(Y~MBW_cspl_X, data=MBW_data, tau=0.75, model=TRUE)
    MBW_data$fit_75 <- predict(model_75)
    model_upc <- rq(Y~MBW_cspl_X, data=MBW_data, tau=upc, model=TRUE)
    MBW_data$fit_upc <- predict(model_upc)
  }
  # If indicated, also fit a comparator model:
  if (comparator==1){
    # Comparator 1 is a quadratic function of age assuming constant variance and winsorising based on a normal assumption:
    MBW_data$X2 <- MBW_data$X^2
    c1model <- lm(Y~X+X2,data=MBW_data)
    MBW_data$c1fit_mn <- predict(c1model)
    MBW_data$c1fit_lpc <- predict(c1model)+qnorm(lpc)*summary(c1model)$sigma
    MBW_data$c1fit_upc <- predict(c1model)+qnorm(upc)*summary(c1model)$sigma
  }
  if (comparator==2){
    # Comparator 2 is simple winsorising within age bands (yearly to 19, then 5 years from 20):
    MBW_data$Xcat[MBW_data$X<20] <- floor(MBW_data$X[MBW_data$X<20])
    MBW_data$Xcat[MBW_data$X>=20] <- 5*floor(MBW_data$X[MBW_data$X>=20]/5)
    MBW_Xcats <- as.numeric(rownames(table(MBW_data$Xcat)))
    for (xcat in MBW_Xcats){
      MBW_data$c2fit_lpc[MBW_data$Xcat==xcat] <- quantile(MBW_data$Y[MBW_data$Xcat==xcat], probs=lpc)
      MBW_data$c2fit_50[MBW_data$Xcat==xcat] <- quantile(MBW_data$Y[MBW_data$Xcat==xcat], probs=0.5)
      MBW_data$c2fit_upc[MBW_data$Xcat==xcat] <- quantile(MBW_data$Y[MBW_data$Xcat==xcat], probs=upc)
    }
  }
  # Calculate the winsorising bounds (equal to the modelled percentiles if multiple=1):
  MBW_data$lb <- MBW_data$fit_50-multiple*(MBW_data$fit_50-MBW_data$fit_lpc)
  MBW_data$ub <- MBW_data$fit_50+multiple*(MBW_data$fit_upc-MBW_data$fit_50)
  # Winsorise any points outside the bounds and save this as a modified variable (BEFORE re-ordering):
  MBW_data$Winsorised <- MBW_data$Y<MBW_data$lb | MBW_data$Y>MBW_data$ub
  MBW_output <- MBW_data$Y
  MBW_output[MBW_output<MBW_data$lb] <- MBW_data$lb[MBW_output<MBW_data$lb]
  MBW_output[MBW_output>MBW_data$ub] <- MBW_data$ub[MBW_output>MBW_data$ub]
  # Report summary results of the winsorisation:
  MBW_data$Winsorised_amount <- MBW_output-MBW_data$Y
  mbw_report <- data.frame(matrix(nrow=0,ncol=2))
  colnames(mbw_report) <- c("Value","Description")
  mbw_report <- rbind(mbw_report,data.frame(Value=length(MBW_data$Winsorised),Description="Total sample size"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[,"X"])[3]),Description="Median X among all data"))
  mbw_report <- rbind(mbw_report,data.frame(Value=sum(MBW_data$Y<MBW_data$lb),Description="Number of values winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=sum(MBW_data$Y>MBW_data$ub),Description="Number of values winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y<MBW_data$lb,"Winsorised_amount"])[1]),Description="Minimum change when winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y<MBW_data$lb,"Winsorised_amount"])[2]),Description="25 %ile of change when winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y<MBW_data$lb,"Winsorised_amount"])[3]),Description="50 %ile of change when winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y<MBW_data$lb,"Winsorised_amount"])[4]),Description="Mean change when winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y<MBW_data$lb,"Winsorised_amount"])[5]),Description="75 %ile of change when winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y<MBW_data$lb,"Winsorised_amount"])[6]),Description="Maximum change when winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y<MBW_data$lb,"X"])[3]),Description="Median X among those winsorised up"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y>MBW_data$ub,"Winsorised_amount"])[6]),Description="Minimum change when winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y>MBW_data$ub,"Winsorised_amount"])[5]),Description="25 %ile of change when winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y>MBW_data$ub,"Winsorised_amount"])[3]),Description="50 %ile of change when winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y>MBW_data$ub,"Winsorised_amount"])[4]),Description="Mean change when winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y>MBW_data$ub,"Winsorised_amount"])[2]),Description="75 %ile of change when winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y>MBW_data$ub,"Winsorised_amount"])[1]),Description="Maximum change when winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(summary(MBW_data[MBW_data$Y>MBW_data$ub,"X"])[3]),Description="Median X among those winsorised down"))
  mbw_report <- rbind(mbw_report,data.frame(Value=knots,Description="Number of knots (If 0, intercept-only model fitted)"))
  if(knots>0){
    for (k in 1:knots){
      mbw_report <- rbind(mbw_report,data.frame(Value=as.numeric(MBW_knotvals[k]),Description=paste0("Position of knot number ",k)))
    }
  }
  print(mbw_report)
  # If indicated, plot the data points and fitted values on a scatterplot, along with the knots:
  if (plotme!="show"&plotme!="save") print("No plot option chosen (valid options are show or save)")
  if (plotme=="show"|plotme=="save"){
    # If indicated, set it up to save the graph as a pdf (and the report as a csv):
    if (plotme=="save"){
      print("Saving plot and report in the working directory")
      pdf("mbw_plot.pdf", width=8, height=7,bg="white", colormodel="cmyk", paper="A4")
      write.table(mbw_report, file="mbw_report.csv", sep=",", row.names=FALSE)
    }
    if (plotme=="show") print("Showing plot in R")
    # Make a dataframe containing fitted values for 1000 points uniformly spaced across the range of X:
    # (This is so that plotted fitted values don't reflect individual X data)
    MBW_plotdata <- data.frame(X=seq(from=min(MBW_data$X),to=max(MBW_data$X),length.out=1000))
    if(knots==0){
      MBW_plotdata$fit_lpc <- Ypc[1]
      MBW_plotdata$fit_50 <- Ypc[2]
      MBW_plotdata$fit_upc <- Ypc[3]
    }
    if(knots>=3&knots<=7){
      MBW_cspl_X <- rcs(MBW_plotdata$X,MBW_knotvals)
      MBW_plotdata$fit_lpc <- predict(model_lpc,MBW_cspl_X)
      MBW_plotdata$fit_50 <- predict(model_50,MBW_cspl_X)
      MBW_plotdata$fit_upc <- predict(model_upc,MBW_cspl_X)
    }
    MBW_plotdata$lb <- MBW_plotdata$fit_50-multiple*(MBW_plotdata$fit_50-MBW_plotdata$fit_lpc)
    MBW_plotdata$ub <- MBW_plotdata$fit_50+multiple*(MBW_plotdata$fit_upc-MBW_plotdata$fit_50)
    if (comparator==1){
      MBW_plotdata$X2 <- MBW_plotdata$X^2
      MBW_plotdata$c1fit_mn <- predict(c1model,MBW_plotdata)
      MBW_plotdata$c1fit_lpc <- predict(c1model,MBW_plotdata)+qnorm(lpc)*summary(c1model)$sigma
      MBW_plotdata$c1fit_upc <- predict(c1model,MBW_plotdata)+qnorm(upc)*summary(c1model)$sigma
      rm(c1model)
    }   
    if (comparator==2){
      MBW_plotdata$Xcat[MBW_plotdata$X<20] <- floor(MBW_plotdata$X[MBW_plotdata$X<20])
      MBW_plotdata$Xcat[MBW_plotdata$X>=20] <- 5*floor(MBW_plotdata$X[MBW_plotdata$X>=20]/5)
      for (xcat in MBW_Xcats){
        MBW_plotdata$c2fit_lpc[MBW_plotdata$Xcat==xcat] <- quantile(MBW_data$Y[MBW_data$Xcat==xcat], probs=lpc)
        MBW_plotdata$c2fit_50[MBW_plotdata$Xcat==xcat] <- quantile(MBW_data$Y[MBW_data$Xcat==xcat], probs=0.5)
        MBW_plotdata$c2fit_upc[MBW_plotdata$Xcat==xcat] <- quantile(MBW_data$Y[MBW_data$Xcat==xcat], probs=upc)
      }
    }
    # Make a jittered version of X (+/- a uniformly distributed value between 0 and jitter_X):
    MBW_data$Xj <- MBW_data$X+runif(n=dim(MBW_data)[1],min=(-1)*jitter_X,max=jitter_X)
    # Make a jittered version of Y (+/- a uniformly distributed value between 0 and jitter_Y*IQR(Y)):
    MBW_data$IQR <- MBW_data$fit_75-MBW_data$fit_25
    MBW_data$randoms <- runif(n=dim(MBW_data)[1],min=(-1),max=1)
    MBW_data$Yj <- MBW_data$Y+MBW_data$randoms*jitter_Y*MBW_data$IQR
    # Do the plot:
    MBW_data <- MBW_data[order(MBW_data$X),]
    plot(MBW_data$Xj,MBW_data$Yj,col="grey",pch=1,xlab="exposure",ylab="outcome",cex=0.5)
    points(MBW_data$Xj[MBW_data$Winsorised],pch=1,MBW_data$Yj[MBW_data$Winsorised],col="red",cex=0.5)
    points(MBW_plotdata$X,MBW_plotdata$fit_50,col="red",lwd=1,type="l")
    points(MBW_plotdata$X,MBW_plotdata$fit_lpc,col="red",lwd=1,type="l",lty=2)
    points(MBW_plotdata$X,MBW_plotdata$fit_upc,col="red",lwd=1,type="l",lty=2)
    points(MBW_plotdata$X,MBW_plotdata$lb,col="red",lwd=1,type="l",lty=2)
    points(MBW_plotdata$X,MBW_plotdata$ub,col="red",lwd=1,type="l",lty=2)
    if (knots>=3&knots<=7) abline(v=MBW_knotvals,lty=2,col="grey")
    # If indicated, also plot the comparator method:
    if (comparator==1){
      points(MBW_plotdata$X,MBW_plotdata$c1fit_mn,col="blue",lwd=1,type="l")
      points(MBW_plotdata$X,MBW_plotdata$c1fit_lpc,col="blue",lwd=1,type="l",lty=2)
      points(MBW_plotdata$X,MBW_plotdata$c1fit_upc,col="blue",lwd=1,type="l",lty=2)
    }
    if (comparator==2){
      points(MBW_plotdata$X,MBW_plotdata$c2fit_50,col="blue",lwd=1,type="l")
      points(MBW_plotdata$X,MBW_plotdata$c2fit_lpc,col="blue",lwd=1,type="l",lty=2)
      points(MBW_plotdata$X,MBW_plotdata$c2fit_upc,col="blue",lwd=1,type="l",lty=2)
    }
    # If saving, finalise the pdf:
    if (plotme=="save") dev.off()
  }
  # Return the modified variable:
  return(MBW_output)
}