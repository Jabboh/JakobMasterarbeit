rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\JakobMasterarbeit//Data")
#install.packages("ape")
#loading packages
library(readxl) 
library(sf)
library(ape) #for Moran's I

#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel("MonthlyData.xlsx")

#Reading in the spatial coordinates of the different trap locations
coords <- read_excel("Traps_coordenadas_geograficas.xls")

#Trap and Area mean the same thing, so we change the name in df from "area" to "trap"
names(df)[names(df)=="Area"] <- "trap"
#Canada is spelt differently, so I change the spelling in df to "Cañada"
df[,"trap"] <- lapply(df[,"trap"], gsub, pattern = "Ca?da", replacement = "Cañada", fixed = TRUE)

#adding lon-lat column to the data frame df
df_new <- merge(df, coords[, c("trap", "Norte", "Oeste")], by="trap", all.x= T)
#I do not have the coordinates for trap "M29" --> one NA in the coordinates >> We need to remove that row
df_new <- df_new[!is.na(df_new$Norte),]

#make an sf-object
projcrs <- 23029
dat <- st_as_sf(x = df_new,                         
               coords = c("Oeste", "Norte"),
               crs = projcrs)

#Make a distance matrix between the observations (this is needed for calculating Moran's I)
trap_d <- as.matrix(st_distance(dat))
trap_d_inv <- 1/trap_d #The zero entries get Inf
#I will set the inf values to 0, because I am not interested whether observation at the same location 
#are autocorrelated, but rather whether there is spatial autocorrelation. The within location correlation
#will be adressed in another script.

#We need to get rid of the units so that we can perform usual matrix operations
attr(trap_d_inv,"units") <- NULL
class(trap_d_inv) <- setdiff(class(trap_d_inv),"units")
trap_d_inv[trap_d_inv == Inf] <- 0

#calculating Moran's I
Moran.I(dat$Cxpippre, trap_d_inv)
#So, there is spatial autocorrelation und nuuuun????? I mean, I did weird stuff, throwing multiple observation 
#of the same location into the pot...
