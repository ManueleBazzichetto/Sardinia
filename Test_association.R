#Analyses done by Manuele Bazzichetto and Vojta Bartak (Czeck University of Life Sciences, Prague)
#Date: 8/12/2022

library(sf)
library(spatstat)
library(spatstat.random)
library(ggplot2)
library(ggpubr)
library(mapview)

#import all packages in a single call
invisible(lapply(c("sf", "spatstat", "spatstat.random", "ggplot2", "ggpubr", "mapview"), require, character.only = T))

#import data: map of vegetation series & human structures (change directory to wherever data are stored)
Sard_poly <- st_read("~/Documents/Sardegna/Ser_veg_diss/Serie_veg_diss.shp")
Veg_map <- st_read("~/Documents/Sardegna/serie_di_vegetazione/Serie_vegetazione_MM_Descr.shp")
Hum_str <- st_read("~/Documents/Sardegna/strutture megalitiche_points/Megalithic structures Sardinia.shp")

#transform CRS to UTM 32N
Sard_poly.utm <- st_transform(Sard_poly, crs = 32632)
Veg_map.utm <- st_transform(Veg_map, crs = 32632)
Hum_str.utm <- st_transform(Hum_str, crs = 32632)

#get rid of original data
rm(Sard_poly, Veg_map, Hum_str)

#have a look at the data
mapview(Veg_map.utm) + mapview(Hum_str.utm)

#check str
str(Hum_str.utm)

#subset Hum_str.utm to keep interesting columns
Hum_str.utm <- Hum_str.utm[c("content_id", "tipology")]

#subset Hum_str.utm to only keep nuraghe and villages
Nuraghe_villagi.utm <- Hum_str.utm[Hum_str.utm$tipology %in% c("nuraghe", "villaggi"), ]

#Create polygon to be used as window for ppp object (see below)
#transfom Sard_poly.utm from MULTIPOLYGON to POLYGON
Sardinia_polys.utm <- st_cast(Sard_poly.utm, "POLYGON")

#keep only NOTE column (which is made of NAs)
Sardinia_polys.utm <- Sardinia_polys.utm[c("NOTE")]

Sard_windows <- as.owin(Sardinia_polys.utm)

#----------------------------------------------------test clustering of Nuraghe + Villaggi

#get rid of problematic points
#NUR11084 (village) and NUR19459 (nuraghe)
#duplicated: NUR4855 duplicated of NUR4780 (the former is validated, the latter is not) -> delete the latter

#rectangular window (rw)
Nuraghe.ppp.rw <- as.ppp(Nuraghe_villagi.utm[!(Nuraghe_villagi.utm$content_id %in% c("NUR19459", "NUR11084", "NUR4780")), ])
sum(duplicated.ppp(Nuraghe.ppp.rw))

#divisor = "d"
plot(pcf.ppp(Nuraghe.ppp.rw, divisor = "d"), main = "Pair Correlation Function")
#divisor = "r"
plot(pcf.ppp(Nuraghe.ppp.rw), main = "Pair Correlation Function")

#Sardinia-shaped window
Nuraghe.ppp.sw <- ppp(x = st_coordinates(Nuraghe_villagi.utm[!(Nuraghe_villagi.utm$content_id %in% c("NUR19459", "NUR11084", "NUR4780")), ])[,1],
                   y = st_coordinates(Nuraghe_villagi.utm[!(Nuraghe_villagi.utm$content_id %in% c("NUR19459", "NUR11084", "NUR4780")), ])[,2], 
                   window = Sard_windows)

sum(duplicated.ppp(Nuraghe.ppp.sw)) #0
#which(duplicated.ppp(Nuraghe.ppp.sw))

#!! THIS TAKES AGES!!
#test clustering with Sardinian-shaped window -> this takes ages
Nuraghe_clust <- pcf.ppp(Nuraghe.ppp.sw, divisor = "d")
plot(Nuraghe_clust, "Pair Correlation Function")

#----------------------------------------------------test association
#shift point pattern keeping fixed the relative distance among points
#then compare observed vs expected number of points falling in vegetation series/categories

set.seed(20)
#use edge = "torus" but rectangular window -> this creates list ~ 20 GB (delete it afterwards)
#lower number for nsim should work anyway
Torus_list <- spatstat.random::rshift.ppp(Nuraghe.ppp.rw, nsim = 5000, edge = "torus")

#from ppp to sf obj st_as_sf()
#example: st_as_sf(Torus_list[[1]])

Min_N <- floor((nrow(Nuraghe_villagi.utm)-3)*.7)

Torus_list <- lapply(Torus_list, function(pnt_obj) {
  pnt_sf <- st_as_sf(pnt_obj)
  if(nrow(pnt_sf) == 1) return(NA)
  st_crs(pnt_sf) <- 32632
  colnames(pnt_sf)[c(1, 2)] <- c("geom_type", "pointID")
  #get rid of poly for window
  pnt_sf <- pnt_sf[pnt_sf$geom_type %in% c("point"), ]
  return(pnt_sf)
})

#check empty elements in Torus
sum(sapply(Torus_list, function(i) class(i)[1]) %in% c("sf")) #5000

#If needed -> get rid of non sf objs
#Torus_list <- Torus_list[sapply(Torus_list, function(i) class(i)[1]) %in% c("sf")] 

#----------------------------------------parametric way

Null_dist_chi.pr <- lapply(Torus_list, function(pnt_obj) {
  #join
  pnt_veg <- st_join(pnt_obj, Veg_map.utm[c("LEG_ITALIA")], join = st_intersects) #LEG_ITALIA is column with ID of vegetation series/categories
  smp_sz <- sum(!is.na(pnt_veg$LEG_ITALIA))
  if(smp_sz < Min_N) return(NA)
  #count pnts in series
  pnt_veg.cnt <- as.data.frame(table(pnt_veg[[c("LEG_ITALIA")]]))
  pnt_veg.cnt <- setNames(pnt_veg.cnt$Freq, pnt_veg.cnt$Var1)
  return(list(Exp = pnt_veg.cnt, N = smp_sz))
})

#Number of useful replicates 5000-2943 = 2057
length(Filter(is.logical, Null_dist_chi.pr)) #2943
#sum(sapply(Null_dist_chi.pr, function(i) is.na(i) && (length(i) == 1)))

#Get rid of NAs
Null_dist_chi.pr <- Filter(Negate(is.logical), Null_dist_chi.pr)

#average sample size 
mean(sapply(Null_dist_chi.pr, '[[', 2)) #4202.667
quantile(sapply(Null_dist_chi.pr, '[[', 2)) #3995, 4329 (1st, 3rd quartiles)

#average across categories
Exp_pr <- sapply(unique(Veg_map.utm$LEG_ITALIA), function(nm) {
  Veg_ser_id <- lapply(Null_dist_chi.pr, function(i) names(i[[1]]))
  Veg_ser_cnt <- Null_dist_chi.pr[which(sapply(Veg_ser_id, function(i) nm %in% i))]
  Veg_ser_cnt <- mean(sapply(Veg_ser_cnt, function(i) i[[1]][[nm]]))
  })

#Get observed number of nuraghe and villages falling in vegetation series/categories
Join_NuragheVill.Veg <- st_join(Nuraghe_villagi.utm[!(Nuraghe_villagi.utm$content_id %in% c("NUR19459", "NUR11084", "NUR4780")), ], Veg_map.utm, join = st_intersects)

ObservedNumNuragheVil <- as.data.frame(table(Join_NuragheVill.Veg$LEG_ITALIA))
ObservedNumNuragheVil <- setNames(object = ObservedNumNuragheVil$Freq, nm = ObservedNumNuragheVil$Var1)

#no points in "SA9"
unique(Veg_map.utm$LEG_ITALIA)[!unique(Veg_map.utm$LEG_ITALIA) %in% names(ObservedNumNuragheVil)]

ObservedNumNuragheVil <- c(ObservedNumNuragheVil, "SA9" = 0)

#compute Chi-sq statistics
Chi_stat <- sum(((ObservedNumNuragheVil - Exp_pr[names(ObservedNumNuragheVil)])^2)/(Exp_pr[names(ObservedNumNuragheVil)]))

qchisq(0.95, df = 29) #42.55697

pchisq(Chi_stat, df = 29, lower.tail = F) #0

#Create df with number of observed vs expected number of points falling in each vegetation category
Obs_vs_exp.df <- data.frame(Veg_class = rep(names(Exp_pr), 2),
                            Cnt = unname(c(ObservedNumNuragheVil[names(Exp_pr)], Exp_pr)))

Obs_vs_exp.df$OE <- rep(c("Obs", "Exp"), each = 30)

#Sort vegetation series/categories for decreasgin number of expected points (for data viz)
Order_veg_ser <- Obs_vs_exp.df[Obs_vs_exp.df$OE == "Exp", "Veg_class"][order(Obs_vs_exp.df[Obs_vs_exp.df$OE == "Exp", "Cnt"], decreasing = T)]

Obs_vs_exp.df$Veg_class <- factor(Obs_vs_exp.df$Veg_class, levels = Order_veg_ser)
Obs_vs_exp.df$OE <- factor(Obs_vs_exp.df$OE, levels = c("Exp", "Obs"))

#Create the bar-chart
ggplot(Obs_vs_exp.df, aes(y = Cnt, x = Veg_class, fill = OE)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Exp" = "#FFC107", "Obs" = "#1E88E5"), name = NULL) +
  ylab("Number of nuraghe and villages") + xlab(NULL) +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(size = 10, angle = 45, vjust = 1, hjust = 1))

#Save the plot
ggsave(filename = "~/Documents/Sardegna/Nuraghe_preference.jpeg", device = "jpeg", units = "cm",
       width = 18, height = 16, dpi = 300)

##REMEMBER TO GET RID OF Torus_list!!
rm(Torus_list)
