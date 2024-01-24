library(viridis)
library(rnaturalearth)
library(sf)
library(ggplot2)

base_map <- function(world, lake_mask=T){
  # world map
  if (lake_mask) { 
    world <- st_read("Admin_shapefiles/world_no_lakes.shp", quiet=TRUE)
  } else { 
    world <- ne_countries(scale = 50, type= "countries", returnclass = "sf")
  }
  
  ggplot(data=world)+
    geom_sf(fill="transparent", color="grey70")+
    theme_bw()+
    coord_sf(xlim = c(-180, 180), ylim = c(-59, 75), expand = FALSE) +
    theme(panel.border = element_blank(),
          axis.title=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          #legend.title = element_text(size=15),
          #legend.text = element_text(size=12), 
          axis.ticks=element_blank(),
          panel.grid.major=element_blank())
}

mapRas <- function(pred_df) { 
  output <- template
  values(output) = as.numeric(pred_df[, 2])
  return(output)
}

plotRaster <- function(raster) { 
  raster_spdf <- as(raster, "SpatialPixelsDataFrame")
  raster_df <- as.data.frame(raster_spdf)
  names(raster_df)[1] <- "value"
  
  min_val <- min(raster_df$value, na.rm = TRUE)
  max_val <- max(raster_df$value, na.rm = TRUE)
  
  p <- geom_tile(data=raster_df, aes(x=x, y=y, fill=value), alpha=0.8)
  
  base_map(world, lake_mask=T) + p + 
    scale_fill_viridis(limits=c(min_val, max_val))
}
