#Functions for analysing phenotypes

#set up packages



read_phenotype_data <- function(phecode, input_dir) {
  # Check if phenotype is present
  filename <- file.path(input_dir, paste0(phecode, ".txt"))
  if(!file.exists(filename)) {
    message(phecode, " not present")
    return(NULL)
  }
  
  # Read in phenotype
  phen <- data.table::fread(filename)
  
  # Check columns are there as expected
  column_names <- c("id", "age", "value")
  if(!all(column_names %in% names(phen))) {
    print(head(phen))
    stop("expect 'id', 'age' and 'value' columns to be present")
  }
  
  # Cut into age bins
  phen$agebin <- cut(phen$age+1, breaks=c(0:19, seq(20, 120, by=5)))
  return(phen)
}

plot_phen_all <- function(phen, phecode) {
  ggplot(phen, aes(x=agebin, y=value)) +
    geom_boxplot() +
    labs(title=phecode)
}

plot_phen_traj <- function(phen, phecode) {
  ggplot(phen, aes(x=agebin, y=value, group=id)) +
    geom_line(aes(group=id), alpha=0.2) +
    labs(title=phecode)
}


detect_outliers <- function(phen, phecode) {
  min <- df %>% 
    filter(pheno_id == phecode) %>%
    select('min')
  min <- as.vector(min)
  max <- df %>%
    filter(pheno_id == phecode) %>%
    select('max')
  max <- as.vector(max)
  
  outliers <- phen %>% 
    filter(!between(phen$value, min, max)) %>% 
    select('value')
  
}


remove_outliers <- function(phen, phecode) {
  min <- df %>% 
    filter(pheno_id == phecode) %>%
    select('min')
  min <- as.vector(min)
  max <- df %>%
    filter(pheno_id == phecode) %>%
    select('max')
  max <- as.vector(max)
  
  outlier_rm <- phen %>% 
    filter(between(phen$value, min, max)) 
  
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



mbw <- function(MBW_Y, MBW_X, knots=5, lpc=0.005, upc=0.995, multiple=1, plotme="show", comparator=0){
  print(paste0("Winsorising at ",multiple," times the interval from the median to percentiles ",lpc," and ",upc))
  # Put the input variables into a dataframe:
  MBW_data <- data.frame(Y=MBW_Y,X=MBW_X)
  
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
  MBW_data$fit_lpc <- fitted(rq(Y~MBW_cspl_X, data=MBW_data, tau=lpc, model=TRUE))
  MBW_data$fit_50 <- fitted(rq(Y~MBW_cspl_X, data=MBW_data, tau=0.5, model=TRUE))
  MBW_data$fit_upc <- fitted(rq(Y~MBW_cspl_X, data=MBW_data, tau=upc, model=TRUE))
  # If indicated, also fit a comparator model:
  if (comparator==1){
    # Comparator 1 is a quadratic function of age assuming constant variance and winsorising based on a normal assumption:
    MBW_data$X2 <- MBW_data$X^2
    c1model <- lm(Y~X+X2,data=MBW_data)
    MBW_data$c1fit_mn <- fitted(c1model)
    MBW_data$c1fit_lpc <- predict(c1model)+qnorm(lpc)*summary(c1model)$sigma
    MBW_data$c1fit_upc <- predict(c1model)+qnorm(upc)*summary(c1model)$sigma
    rm(c1model)
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
  # Calculate the winsorising bounds:
  MBW_data$lb <- MBW_data$fit_50-multiple*(MBW_data$fit_50-MBW_data$fit_lpc)
  MBW_data$ub <- MBW_data$fit_50+multiple*(MBW_data$fit_upc-MBW_data$fit_50)
  # Winsorise any points outside the bounds and save this as a modified variable (BEFORE re-ordering):
  MBW_data$Winsorised <- MBW_data$Y<MBW_data$lb | MBW_data$Y>MBW_data$ub
  MBW_output <- MBW_data$Y
  MBW_output[MBW_output<MBW_data$lb] <- MBW_data$lb[MBW_output<MBW_data$lb]
  MBW_output[MBW_output>MBW_data$ub] <- MBW_data$ub[MBW_output>MBW_data$ub]
  # Report the winsorisation:
  print(paste0(sum(MBW_data$Winsorised)," of ",length(MBW_data$Winsorised)," data points winsorised"))
  print(head(MBW_data[MBW_data$Winsorised,],n=20))
  # If indicated, plot the fitted values on a scatterplot, along with the knots:
  if (plotme!="show"&plotme!="save") print("No plot option chosen (valid options are show or save)")
  if (plotme=="show"|plotme=="save"){
    # If indicated, set it up to save as a pdf:
    if (plotme=="save"){
      print("Saving plot in the working directory")
      pdf("Winsorising_plot.pdf", width=8, height=7,bg="white", colormodel="cmyk", paper="A4")
    }
    # Do the plot:
    if (plotme=="show") print("Showing plot in R")
    MBW_data <- MBW_data[order(MBW_data$X),]
    plot(MBW_data$X,MBW_data$Y,col="grey",pch=1,xlab="exposure",ylab="outcome",cex=0.5)
    points(MBW_data$X[MBW_data$Winsorised],pch=1,MBW_data$Y[MBW_data$Winsorised],col="red",cex=0.5)
    points(MBW_data$X,MBW_data$fit_50,col="red",lwd=1,type="l")
    points(MBW_data$X,MBW_data$fit_lpc,col="red",lwd=1,type="l",lty=2)
    points(MBW_data$X,MBW_data$fit_upc,col="red",lwd=1,type="l",lty=2)
    points(MBW_data$X,MBW_data$lb,col="red",lwd=1,type="l",lty=2)
    points(MBW_data$X,MBW_data$ub,col="red",lwd=1,type="l",lty=2)
    abline(v=MBW_knotvals,lty=2,col="grey")
    # If indicated, also plot the comparator method:
    if (comparator==1){
      points(MBW_data$X,MBW_data$c1fit_mn,col="blue",lwd=1,type="l")
      points(MBW_data$X,MBW_data$c1fit_lpc,col="blue",lwd=1,type="l",lty=2)
      points(MBW_data$X,MBW_data$c1fit_upc,col="blue",lwd=1,type="l",lty=2)
    }
    if (comparator==2){
      points(MBW_data$X,MBW_data$c2fit_50,col="blue",lwd=1,type="l")
      points(MBW_data$X,MBW_data$c2fit_lpc,col="blue",lwd=1,type="l",lty=2)
      points(MBW_data$X,MBW_data$c2fit_upc,col="blue",lwd=1,type="l",lty=2)
    }
    # If saving, finalise the pdf:
    if (plotme=="save") dev.off()
  }
  # Return the modified variable and remove temporary variables:
  return(MBW_output)
  rm(list=c("MBW_data","MBW_knotvals","MBW_cspl_X","MBW_output"))
}