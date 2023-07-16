# Example to calculate mean values over a moving window
library(climate4R.UDG)
library(loadeR)
library(transformeR)
library(gtools)
library(openxlsx)
library(tm)
c.ctl = as.data.frame(matrix(NA, nrow = 1, ncol = 1))
c.std = as.data.frame(matrix(NA, nrow = 1, ncol = 1))
#c.units = as.data.frame(matrix(NA, nrow = 1, ncol = 1))
files = list.files(path = getwd(), pattern='pr', all.files=TRUE)
#load("CORDEX-WAS-44_CCCma-CanESM2_historical_r1i1p1_IITM-RegCM4-4_v5_pr_1986_2005.Rdata")
coords <- list()
coords$x = 74.18
coords$y = 36.17
seql = seq(1:16)
j = 1
for (j in 1:length(seql)) {
  us <- files[seql[j]]
  dset <- load(files[seql[j]])
  #eobs1 <- bindGrid(hist1,hist2, dimension = "time")
  CORDEX_hist_pr <- interpGrid(bg, coords)
  getGridUnits(CORDEX_hist_pr)
  # Convert PPT from kg m-2 s-1 to mm
  #CORDEX_hist_pr <- gridArithmetics(CORDEX_hist_pr, 86400, operator="*")
  #c1 <- getGridUnits(dd)
  ##########################
  # create the vector of doys, with month and day
  doy <- substr(seq(as.Date("1986-03-01"), as.Date("1986-12-31"), by="days"),6,10) # exemplary 366-days year
  # you can also use here only the period you are working with: 1980-03-01 to 1980-10-30
  # get month and day of the full series
  my.dates <- substr(CORDEX_hist_pr$Dates$start,6,10)
  ##########################
  # prepare the calculation of indices for the desired days
  mw.length <- 31 # choose the moving window length
  mw <- floor(mw.length/2)
  index.doy <- list() # initialize list for each doy
  i = 16
  for (i in (mw+1):(length(doy)-mw-1)){
    index.doy[[i]] <- which(my.dates %in% doy[(i-mw):(i+mw)]) # these are the same for all simulations for a fixed period
  }
  ##########################
  # the following needs to be calculated for each simulation
  ls <-lapply((mw+1):(length(doy)-mw-1), function(i){
    dd.sub <- subsetDimension(CORDEX_hist_pr, dimension = "time", indices=index.doy[[i]]) 
    dd.clim <- climatology(dd.sub)
    return(dd.clim)
  })
  clim.ctl <- bindGrid(ls, dimension = "time")
  ls.std <-lapply((mw+1):(length(doy)-mw-1), function(i){
    dd.sub <- subsetDimension(CORDEX_hist_pr, dimension = "time", indices=index.doy[[i]]) 
    dd.std <- climatology(dd.sub, clim.fun = list(FUN = "sd", na.rm = TRUE))
    return(dd.std)
  })
  clim.std <- bindGrid(ls.std, dimension = "time")
  #Output path
  #OutPath<- "F:/CORDEX_WAS44_Naltar_Mean/Historical_1986-2005/Precipitation/"
  #dir.create(OutPath)
  c.ctl = cbind(c.ctl,clim.ctl$Data)
  c.std = cbind(c.std,clim.std$Data)
  #c.units = cbind(c.units,c1)
  save(clim.ctl ,file =  paste0("/",us,"MA_CTL.Rdata"))
  save(clim.std ,file =  paste0("/",us,"MA_STD.Rdata")) 
  #save(clim.ctl ,file =  paste0(us,"MA_CTL.Rdata"))
  #save(clim.std ,file =  paste0(us,"MA_STD.Rdata"))
  write.csv(c.ctl,file =  paste0("/","hist_1991-2010_CTL.csv"))
  write.csv(c.std,file =  paste0("/","hist_1991-2010_STD.csv"))
}




##########################
# Do the same for the scenario simulation (e.g. clim.rcp). Consider that index.doy changes when the period is different.
# Compute additive or multiplicative deltas using:
#deltaADD <-gridArithmetics(clim.rcp, clim.ctl, operator="+" )
# deltaMUL <-gridArithmetics(clim.rcp, clim.ctl, operator="/" )

# ADDITIVE and Multipilicative Delta Change
clim.ctl_files_hist <- list.files(path = '', pattern='.MA_CTL', all.files=TRUE)
clim.std_files_hist <- list.files(path = '', pattern='.MA_STD', all.files=TRUE)
clim.ctl_files_rcp <- list.files(path = '/', pattern='.MA_CTL', all.files=TRUE)
clim.std_files_rcp <- list.files(path = '/', pattern='.MA_STD', all.files=TRUE)

c_ctl = as.data.frame(matrix(NA, nrow = 1, ncol = 1))
c_std = as.data.frame(matrix(NA, nrow = 1, ncol = 1))
i = 1
#j = 1
seql = seq(1:16)
dir1=getwd()
dir2="/"
for (i in 1: length(clim.ctl_files_hist)) {
  # read historical cilimotolgy files
  us_ctl <- clim.ctl_files_hist[seql[i]]
  origonal_gridded = load(paste0(dir1,'/',clim.ctl_files_hist[seql[i]]))
  historical_ctl = clim.ctl
  # read historical standard devition files
  us_std <- clim.std_files_hist[seql[i]]
  origonal_gridded = load(paste0(dir1,'/',clim.std_files_hist[seql[i]]))
  historical_std = clim.std
  # read future cilimotolgy files
  CF_gridded = load(paste0(dir2,'/',clim.ctl_files_rcp[seql[i]]))
  future__ctl_rcp = clim.ctl
  # read future standard deviation files
  CF_gridded = load(paste0(dir2,'/',clim.std_files_rcp[seql[i]]))
  future__std_rcp = clim.std
  
  #deltaADD <-gridArithmetics(future_rcp, historical, operator="-" )
  deltaMUL_ctl <-gridArithmetics(future__ctl_rcp, historical_ctl, operator="/" )
  deltaMUL_std <-gridArithmetics(future__std_rcp, historical_std, operator="/" )
  # make dataframe for climotology and standard deavition
  c_ctl = cbind(c_ctl,deltaMUL_ctl$Data)
  c_std = cbind(c_std,deltaMUL_std$Data)
  
  write.csv(c_ctl,file =  paste0("/","SC_CTL_FF(45)_pr_2080_2099.csv"))
  write.csv(c_std,file =  paste0("/","SC_STD_FF(45)_pr_2080_2099.csv"))
  #deltaMUL <-gridArithmetics(future_rcp, historical, operator="/" )
  save(deltaMUL_ctl ,file =  paste0("/",us_ctl,"SC(CTL).Rdata"))
  save(deltaMUL_std ,file =  paste0("/",us_std,"SC(STD).Rdata"))
}







  ls <-lapply((mw+1):(length(doy)-mw-1), function(i){
  dd.sub <- subsetDimension(CORDEX_hist_pr, dimension = "time", indices=index.doy[[i]]) 
  dd.std <- climatology(dd.sub, clim.fun = list(FUN = "std", na.rm = TRUE))
  return(dd.std)
})
