#-------------------------------------------------------------------------------
getwd() #check dir
source("CHIK_mapping/Functions/BurdenFunctions.R")
load("rasterDf.RData")
load("total_infection.RData")
age_vec <- read.csv('age_structure.csv')
#-------------------------------------------------------------------------------

#1. bring total population count data ------------------------------------------

pop_count <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/final/Worldpop_population_2020_5k.tif"
pop_count <- "D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping/Raster/landscan-global-2022.tif"

pop_count <- raster(pop_count)
plot(pop_count)

pop_count  <- projectRaster(pop_count, elev, method="bilinear")
plot(pop_count)

pop_spdf <- as(pop_count, "SpatialPixelsDataFrame")
pop_df <- as.data.frame(pop_count, xy = TRUE, na.rm = T)

#2. combine population with original FOI raster data----------------------------
foi_comb <- merge(raster_df, pop_df, by = c("x", "y"))
foi_comb$pop <- foi_comb$Worldpop_population_2020_5k * 5000

#3. combine country data -------------------------------------------------------
admin_boundaries <- st_read("World_shape/world-administrative-boundaries.shp")
admin_boundaries_sf <- st_as_sf(admin_boundaries, wkt = "geometry", crs = 4326)
foi_comb <- st_as_sf(foi_comb, coords = c("x", "y"), crs = 4326)
foi_comb <- st_join(foi_comb, admin_boundaries_sf, join = st_nearest_feature)

#4. burden calculation: total infection ----------------------------------------
burden <- age_strat_burden(age_vec, foi_comb)
burden_df <- burden$updated_df
burden_by_age <- as.data.frame(burden$infection_per_band)
burden_by_age$tot_infection <- rowSums(burden_by_age)

# change colnames
col_names <- colnames(age_vec)[4:ncol(age_vec)]
colnames(burden_by_age)[1:20] <- col_names

# combine with burden_df
total_infection_global <- cbind(burden_df, burden_by_age)

save("total_infection_global", "burden_df", "burden_by_age", file = "total_infection.RData")

#5. burden map 

total_infection_spatial <- as(total_infection_global, "Spatial")
TotInfspdf <- as(total_infection_spatial, "SpatialPixelsDataFrame")
TotInfdf <- as.data.frame(TotInfspdf)
crs(temp) <- crs(TotInfspdf)
total_infection_raster <- rasterize(TotInfspdf, temp, field="total_infection")
burdenmask <- total_infection_raster * mask_all

burden_spdf <- as(burdenmask, "SpatialPixelsDataFrame")
burden_df <- as.data.frame(burden_spdf)
names(burden_df)[1] <- "value"

min_val <- min(burden_df$value, na.rm = TRUE)
max_val <- max(burden_df$value, na.rm = TRUE)

p <- geom_tile(data=burden_df, aes(x=x, y=y, fill=value), alpha=0.8)

base_map(world, lake_mask=T) + p +
  scale_fill_viridis(option = "rocket", direction = -1, limits = c(0, max_val))


ggplot(data=burden_df, aes(x=x, y=y, fill=log(value + 1))) +
  geom_tile(alpha=0.8) +
  coord_fixed() +
  scale_fill_viridis(option = "rocket", direction = -1) 