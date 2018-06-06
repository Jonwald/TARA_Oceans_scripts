##########################
# load required libraries
##########################

install.packages("tmap")
install.packages("ggmap")
library(ggmap)
library(maptools)
library(maps)
library(reshape2)
display.brewer.all()

#######################################
# load data, (change paths as required)
#######################################

setwd("/media/jonny/TARA/ref_genes_only")
sample_meta <- read.csv("sample_meta.csv", sep="\t", header= T, stringsAsFactors = F)
all_cogs <- read.csv("TARA243.OG.profile.release", sep ="\t", header =T)
### all_taxa <- to be added

#####################################
# set up functions
#####################################

# select_cog: filters out single cog of interest from dataset and merge to metadata
# usage  select_cog("enter cog id here")

select_cog <- function(x){
  a <- subset(all_cogs, cog==x)
  b <- cbind(colnames(a[-1]), as.data.frame(t(as.matrix(a[-1]))))
  colnames(b) <- c("Sample", "abund")
  tmp <- full_join(sample_meta, b, by = "Sample")
  tmp$Environment <- factor(tmp$Environment, levels = c("surface water layer", "deep chlorophyll maximum layer", "marine epipelagic mixed layer", "mesopelagic zone"))
  tmp
  }

## map_cog: maps cog of interest from size fraction of interest using hexagonal bins,
## where multiple samples are take at a single a site the average abundance is used
## modify function line "fun=mean" to change this behaviour if desired

## Usage: map_cog("dataframe from output of select cog", "fraction: either Prokaryote, virus, or Girus")

## Note: environment "marine epipleagic mixed layer is removed due to low number of samples causing weird 
## behaviour during plotting, can be included by changing Environment != "marine epipelagic mixed layer" 
## to Environment != " ".

map_cog <- function(x,y){
  fract <- subset(x, Fraction_type == y)
  fract_sub <- subset(fract, Environment != "marine epipelagic mixed layer")
  mp <- NULL
  mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
  mp <- ggplot() + mapWorld
  mp +
    coord_cartesian() +
    stat_summary_hex(data=fract_sub, aes(x=Long, y=Lat, z=abund), fun=mean, alpha=0.9) +
    facet_wrap(~ Environment, nrow=2) +
    scale_fill_distiller(palette = "Spectral") +
    scale_x_continuous(limits = c(-180, +180)) +
    scale_y_continuous(limits = c(-90, +90)) +
    theme_gray() +
    ggtitle(paste("Relative Abundance Heatmap for",COG,"in TARA Oceans Metagenomes")) +
    theme(plot.title = element_text(hjust = 0.5))
}

## make boxplots: make boxplots for abundance of cog grouping by marine layer

make_boxplots_w_stats <- function(x) {
  ggboxplot(x, x="Environment", y="abund", fill="Environment",
            add = "jitter",
            add.params= list(alpha=0.7),
            outlier.shape = NA,
            legend="none",
            title =paste("Relative Abundance of",COG, "By Aquatic Layer")) +
    rotate_x_text(angle = 45) +
    stat_compare_means()
    }

make_depth_scatter <- function(x){
  ggplot(x, aes(x=Depth_m, y=abund)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_y_log10() +
    scale_x_log10() +
    theme(legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust =1),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 14),
          panel.border = element_rect(fill = NA, size=0.5),
          strip.background = element_rect(color = "black", size = 0.5)) +
  ggtitle(paste("Relative Abundance of",COG,"vs Sample Depth(m)")) +
    theme(plot.title = element_text(hjust = 0.5))
  }


make_depth_scatter_by_biome <- function(x){
  ggplot(x, aes(x=Depth_m, y=abund)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_y_log10() +
    scale_x_log10() +
    theme(legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust =1),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 14),
          panel.border = element_rect(fill = NA, size=0.5),
          strip.background = element_rect(color = "black", size = 0.5)) +
    facet_wrap(~ Biome, scales = "free_x") +
    ggtitle(paste("Relative Abundance of",COG,"vs Sample Depth(m), Split by Biome Type")) +
    theme(plot.title = element_text(hjust = 0.5))
}

make_depth_scatter_by_region <- function(x){
  ggplot(x, aes(x=Depth_m, y=abund)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_y_log10() +
    scale_x_log10() +
    theme(legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=14),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust =1),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 11),
          panel.border = element_rect(fill = NA, size=0.5),
          strip.background = element_rect(color = "black", size = 0.5)) +
    facet_wrap(~ Region, scales = "free_x") +
    ggtitle(paste("Relative Abundance of",COG,"vs Sample Depth(m), Split by Ocean Region")) +
    theme(plot.title = element_text(hjust = 0.5))
          }


#########################################
# analysis of particular genes goes here
#########################################

COG <- "NOG39321" # enter COG of interest here

pdf(paste(COG,".pdf"))
dat <- select_cog(COG)
map_cog(dat, "Prokaryote")
make_boxplots_w_stats(dat)
make_depth_scatter(dat)
make_depth_scatter_by_biome(dat)
make_depth_scatter_by_region(dat)
dev.off()






#############################
# additional features to add
#############################

# Add ability to parse tara taxonomy
# add ability to prodcue heta-tree from taxonomy








