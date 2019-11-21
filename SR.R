
# upload the necissary data and import the dataset after 
# forking the data from gentry's github
packages<-c("mapproj", "cowplot", "dplyr", "geosphere", "ggplot2", "ggExtra", "maps", "maptools", "readxl", "rgdal", "rgeos", "sf", "sp", "spatialreg", "spdep", "tidyr", "viridis")
sapply(packages, require, character.only=T)

data <- read.csv('./Data/childpov18_southfull.csv', 
                 colClasses = c("character", "character", "character", 
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", "numeric",
                                "numeric", "numeric", "numeric", "numeric"))


head(data)

#############################################################
#change the names that started with dates that R replaxed 
# with x"" virsions. then subset the data to the state you
# want to use for the analysis (LA)
names(data)[names(data)=="X2016.child.poverty"] <- "child.pov.2016"

la_pov <- data %>% subset(State == "LA")

summary(la_pov)

##############################################################
# Ordinary Least Squares equation and lm model to determine 
# significant variables associated with child poverty
equation <- child.pov.2016 ~ rural + urban + lnmanufacturing + lnag + 
  lnretail + lnhealthss + lnconstruction + lnlesshs + 
  lnunemployment + lnsinglemom + lnblack + lnhispanic + 
  lnuninsured + lnincome_ratio + lnteenbirth + lnunmarried

options(scipen = 5)

ols <- lm(equation, data=la_pov)
summary(ols)
# lnconstruction and lnunemployment are both significant 
# with Inlesshs close. rural, lnhealthss, and lnteenbirth 
# may also be worth examining even though they are 
# non-significant.

###############################################################
# creat contiguity of neighbors within state and between the counties

#Obtain FIPS Codes by county 
fips <- county.fips

#Create county polygons
louisiana <- map(database = "county", regions = "louisiana", fill=T, plot=F)
IDs <- sub("^louisiana,","",louisiana$names)

#Add FIPS codes to the county polygons
fips.codes <- separate(data = fips, col = polyname, into = c("state", "county"), sep = ",")
la_fips <- subset(fips.codes, state=="louisiana", select=fips)
names <- fips.codes$county
la_IDs <- unique(la_fips$fips)

#Create spatial polygons
la_sp = map2SpatialPolygons(louisiana,la_fips$fips,CRS("+proj=longlat"))
names(la_sp@polygons) <- la_IDs

#Create neighbor weights using the queens case
neighb.data <- poly2nb(la_sp, queen=T)
names(neighb.data) <- names(la_sp@polygons)

#Create list of neighbors
cont.neighb <- nb2listw(neighb.data,style="W", zero.policy = TRUE)

##############################################################

lm.morantest(ols, cont.neighb)
# here we use Moran's Correlation to examine if their is 
# a spatial correlation between our residuals based on the 
# variables we selected. In this test we return a P-value 
# of 0.4027 meaning that there does not appear to be 
# significant spatial corellation between the residuals 
# for this subset of the data.

##############################################################

lm.LMtests(ols, cont.neighb, test="all")
# Here we run a LaGrange Multiplier test to give us a 
# simple estimate of LM tests for error dependence and 
# missing spatially lagged dependent variable (LMerr & 
# LMlag). In this case an LMlag appears to be the only 
# significant test based off of the P-value (0.01653).

# the R simply indicates robust or not, for the purposes 
#of this example we will focus on the non-robust versions. 
#this test also does not imply X or Y so you will need to 
#run both to determine which model works best.

##############################################################
# Spatial lagged X model which accounts for the 
# neighboring X values within the model.
SLX.model <- lmSLX(equation, data=la_pov, cont.neighb)
summary(SLX.model)

summary(impacts(SLX.model, cont.neighb), zstats = TRUE)[["pzmat"]]
# the SLX.model is the best fist based on the P-value 
# of 0.0000002991

################################################################
# Spatial lagged Y model which accounts for the 
# neighboring Y values within the model.
sp.lag.model <- spatialreg::lagsarlm(equation, data=la_pov, cont.neighb)
summary(sp.lag.model, Nagelkerke = TRUE)
#Where Rho is the spatially lagged y multiplier

summary(impacts(sp.lag.model, listw = cont.neighb, R=100), zstats = TRUE)[["pzmat"]]
#Because this is a simulation, R details the number of repetitions.

################################################################
# this model does not include lagged value and examines 
# if . . . . I actually dont really undersatnd this
sp.err.model <- spatialreg::errorsarlm(equation, data=la_pov, cont.neighb)
summary(sp.err.model, Nagelkerke = TRUE)
#Where Lambda is the error multiplier

spatialreg::Hausman.test(sp.err.model)

##############################################################
###############################
############################### figure out what to do at this point since he moved forward with the error model and you need to use the Xlm
###############################
sd.err <- spatialreg::errorsarlm(equation, la_pov, cont.neighb, etype = "emixed")
sdm <- spatialreg::lagsarlm(equation, la_pov, cont.neighb, type = "mixed")
summary(sd.err, Nagelkerke = TRUE)

#or
summary(spatialreg::impacts(sd.err, listw = cont.neighb, R = 100), zstats = TRUE)[["pzmat"]]

spatialreg::LR.sarlm(sd.err,sp.err.model)
#Where our df is equal to the number of lagged variables

##############################################################################

all.xy <- centroid(la_sp)
#tx_IDs <- unique(tx_fips$fips) this value was created in the contiguity section but would be needed here if only using distance functions. See "creating list of contiguity neighbors" for details.
rownames(all.xy) <- la_IDs
colnames(all.xy) <- cbind("x","y")

#Create neighbors
all.dist.k1 <- knn2nb(knearneigh(all.xy, k=1, longlat = TRUE))
all.dist.k3 <- knn2nb(knearneigh(all.xy, k=3, longlat = TRUE))
all.dist.k5 <- knn2nb(knearneigh(all.xy, k=5, longlat = TRUE))

#Determine max k distance value to neighbor
all.max.k1 <- max(unlist(nbdists(all.dist.k1, all.xy, longlat=TRUE)))
all.max.k3 <- max(unlist(nbdists(all.dist.k3, all.xy, longlat=TRUE)))
all.max.k5 <- max(unlist(nbdists(all.dist.k5, all.xy, longlat=TRUE)))

#Calculate neighbors based on distance
all.sp.dist.k1 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k1, longlat = TRUE)
all.sp.dist.k3 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k3, longlat = TRUE)
all.sp.dist.k5 <- dnearneigh(all.xy, d1=0, d2=1 * all.max.k5, longlat = TRUE)

#Create neighbor list
all.dist.neighb.k1 <- nb2listw(all.sp.dist.k1,style="W", zero.policy = TRUE)
all.dist.neighb.k3 <- nb2listw(all.sp.dist.k3,style="W", zero.policy = TRUE)
all.dist.neighb.k5 <- nb2listw(all.sp.dist.k5,style="W", zero.policy = TRUE)

##############################################################################

all.dist.lag.k1 <- lagsarlm(equation, data = la_pov, listw = all.dist.neighb.k1)
all.dist.lag.k3 <- lagsarlm(equation, data = la_pov, listw = all.dist.neighb.k3)
all.dist.lag.k5 <- lagsarlm(equation, data = la_pov, listw = all.dist.neighb.k5)

summary(all.dist.lag.k1, Nagelkerke = TRUE)

###################################################################

all.dist.err.k1 <- errorsarlm(equation, data = la_pov, listw = all.dist.neighb.k1)
all.dist.err.k3 <- errorsarlm(equation, data = la_pov, listw = all.dist.neighb.k3)
all.dist.err.k5 <- errorsarlm(equation, data = la_pov, listw = all.dist.neighb.k5)

summary(all.dist.err.k1, Nagelkerke = TRUE)

##############################################################

dist.err.data <- summary(all.dist.err.k1, correlation=TRUE, Nagelkerke = TRUE)

dist.err.output <- cbind.data.frame(la_pov$FIPS,
                                    dist.err.data$fitted.values, 
                                    dist.err.data$residual, 
                                    la_pov$child.pov.2016, 
                                    la_pov$lnconstruction, 
                                    la_pov$lnunemployment, 
                                    la_pov$lnlesshs, 
                                    la_pov$rural, 
                                    la_pov$lnteenbirth, 
                                    la_pov$lnhealthss,
                                    stringsAsFactors = FALSE)

#Renaming columns
colnames(dist.err.output) <- c("fips","fitted","resid","childpov",
                               "construction","unemployed","less_hs","rural","teen_birth","no_healthins")


#Create quantiles
quantiles_sm <- dist.err.output %>%
  pull(construction) %>%
  quantile(probs = seq(0, 1, length.out = 4), na.rm = TRUE)

quantiles_pov <- dist.err.output %>%
  pull(childpov) %>%
  quantile(probs = seq(0, 1, length.out = 4), na.rm = TRUE)

#Create ranks
sm_rank <- cut(dist.err.output$construction, 
               breaks= quantiles_sm, 
               labels=c("1", "2", "3"), 
               na.rm = TRUE, 
               include.lowest = TRUE)

pov_rank <- cut(dist.err.output$childpov, 
                breaks= quantiles_pov, 
                labels=c("1", "2", "3"), 
                na.rm = TRUE,
                include.lowest = TRUE)

#Join ranks and combined column to dataset
dist.err.output$cons_score <- as.numeric(sm_rank)
dist.err.output$pov_score <- as.numeric(pov_rank)
dist.err.output$cons_pov <- paste(as.numeric(dist.err.output$pov_score), 
                                 "-", 
                                 as.numeric(dist.err.output$cons_score))

###################################################################################

legend_colors <- tibble(
  x = c(3,2,1,3,2,1,3,2,1),
  y = c(3,3,3,2,2,2,1,1,1),
  z = c("#574249", "#627f8c", "#64acbe", "#985356", "#ad9ea5", "#b0d5df", "#c85a5a", "#e4acac", "#e8e8e8"))

xlabel <- "Poverty,Low \u2192 High"
xlabel <- gsub(",", "\n", xlabel)
ylabel <- "Construction Work Household,Low \u2192 High"
ylabel <- gsub(",", "\n", ylabel)

legend <- ggplot(legend_colors, aes(x,y)) + 
  geom_tile(aes(fill=z)) + 
  theme_minimal() + theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x = xlabel, y = ylabel) + 
  scale_fill_identity() +
  ggtitle("Legend") +
  theme(axis.title.y = element_text(face = "italic", hjust = 0.5, size = 8)) +
  theme(axis.title.x = element_text(face = "italic", hjust = 0.5, size = 8)) +
  theme(plot.title = element_text(face="bold", hjust = 0.5, size = 10))

####################################################################################

world <- map_data("world")
states <- map_data("state")
counties <- map_data("county")

counties$polyname <- paste(counties$region, counties$subregion, sep = ",")
counties <- counties %>% left_join(fips, by = c("polyname" = "polyname"))
counties$fips <- as.character(counties$fips)
counties <- counties %>% left_join(la_pov, by = c("fips" = "FIPS"))

southern_states <- subset(states, region %in% 
                            c("texas", "arkansas", "louisiana", "mississippi", 
                              "alabama", "georgia", "florida", "north carolina",
                              "south carolina", "tennessee", "oklahoma", 
                              "kentucky", "west virginia", "virginia", 
                              "maryland", "delaware", "district of columbia"))

southern_counties <- subset(counties, region %in% 
                              c("texas", "arkansas", "louisiana", "mississippi", 
                                "alabama", "georgia", "florida", "north carolina",
                                "south carolina", "tennessee", "oklahoma", 
                                "kentucky", "west virginia", "virginia", 
                                "maryland", "delaware", "district of columbia"))

louisiana_counties <- subset(southern_counties, region == "louisiana")

#########################################################################

#Attach the data via the FIPS column and fortify the polygon
la_poly <- louisiana_counties %>% 
  left_join(dist.err.output, by = c("fips" = "fips")) %>%
  fortify

#Add custom color scheme based on ranks
bivariate_color_scale <- tibble(
  "3 - 3" = "#574249", 
  "2 - 3" = "#627f8c",
  "1 - 3" = "#64acbe",
  "3 - 2" = "#985356",
  "2 - 2" = "#ad9ea5",
  "1 - 2" = "#b0d5df",
  "3 - 1" = "#c85a5a",
  "2 - 1" = "#e4acac",
  "1 - 1" = "#e8e8e8") %>%
  gather("group", "fill")

la_poly <- la_poly %>% 
  left_join(bivariate_color_scale, by = c("cons_pov" = "group"))

####################################################################

mom_pov_map <- ggplot() + 
  geom_polygon(data = world, aes(x=long,y=lat, group=group), fill = "gray95", color = "white") +
  geom_polygon(data = states, aes(x=long,y=lat, group=group), fill = "gray", color = "white") +
  geom_polygon(data = la_poly, aes(x=long, y=lat, group=group, fill = fill)) + 
  geom_polygon(data = southern_states, aes(x=long,y=lat, group=group), fill = NA, color = "white") +
  geom_polygon(data = louisiana_counties, aes(x=long,y=lat, group=group), fill = NA, color = "black", size = 0.05) +
  coord_map("conic", lat0 = 30, xlim=c(-95,-88), ylim=c(28,34)) +
  scale_fill_identity() +
  theme_grey() + theme(legend.position="bottom") + theme(legend.title.align=0.5) +
  theme(panel.background = element_rect(fill = 'deepskyblue'),
        panel.grid.major = element_line(colour = NA)) +
  labs(x = "Longitude", y = "Latitude", fill = "Child Poverty", 
       title = "Bivariate Map of Child Poverty and Construction Work Hourseholds") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
#mom_pov_map use to preview the map

###################################################################

final_map <- ggdraw() +
  draw_plot(mom_pov_map, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(legend, x = 0.65, y = 0.60, width = 0.2, height = 0.25) 

final_map






















