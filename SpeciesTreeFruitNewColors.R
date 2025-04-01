# load required packages
library(ape)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(stringr)
library(gridExtra)
library(grid)
library(ggtext)

# Set complimentary color palettes for host-of-origin and clade affinities



hostcolors <- c(
  "mediumpurple4", #Avena
  "#E31A1C", # Bromus
  "green4", # Cenchrus
  "#20B2AA", # Cynodon
  "lightsalmon", # Digitaria,
  "goldenrod4", # Echinochloa
  "olivedrab3", # "Eleusine
  "#FB9A99", # Elionorus
  "firebrick1", # Eragrostis
  "khaki1", # Hakonechloa
  "darkorange1", # Leersia
  "yellow3", #Leptochloa
  "mediumpurple1", # Lolium
  "tan", # Luziola
  "azure4", # Melinis
  "#005F5F", #Oryza
  "deeppink1", # "Panicum"
  "lightgray", #"Paspalum"
  "steelblue1", # Setaria
  "orchid1", # Stenotaphrum
  "royalblue", #Triticum
  "maroon", # Urochloa
  "darkorange4" #zea
)

names(hostcolors) <- scan(text = "Avena Bromus Cenchrus Cynodon Digitaria Echinochloa Eleusine Elionorus Eragrostis Hakonechloa Leersia Leptochloa Lolium Luziola Melinis Oryza Panicum Paspalum Setaria Stenotaphrum Triticum Urochloa Zea", what = "")

teals <- c( "#005F5F", "#008080", "#20B2AA","#40E0D0", "#7DF9FF")

cladecolors <- c(
  "#20B2AA",  # Cynodon1
  "#7DF9FF", # Cynodon2
  "goldenrod4", # Echinochloa
  "darkolivegreen", # Eleusine1
  "olivedrab3", # Eleusine2
  "palegreen",# Eleusine3
  "firebrick2", # Eragrostis
  "darkorange1", # Leersia
  "mediumpurple1", # Lolium
  "mediumpurple4", # Lolium2,
  "purple4", # Lolium3
  "azure4", # Melinis
  "#005F5F", # Oryza
  "deeppink1", # Panicum
  "navyblue", #Panicum2
  "steelblue1", # Setaria
  "orchid1", # Stenotaphrum
  "royalblue", #Triticum
  "palevioletred", # Urochloa1
  "rosybrown1",# Urochloa2
  "maroon" #Urochloa3
)

names(cladecolors) <- scan(text = "Cynodon1 Cynodon2 Echinochloa Eleusine1 Eleusine2 Eleusine3 Eragrostis Leersia Lolium Lolium2 Lolium3 Melinis Oryza Panicum Panicum2 Setaria Stenotaphrum Triticum Urochloa1 Urochloa2 Urochloa3", what = "")

cladeOrder <- c("Oryza", "Leersia", "Panicum2", "Setaria", "Panicum", "Urochloa2", "Cynodon1", "Lolium3", "Stenotaphrum", "Urochloa3", "Urochloa1", "Cynodon2", "Melinis", "Eleusine3", "Echinochloa", "Lolium2", "Eragrostis", "Eleusine1", "Eleusine2", "Triticum", "Lolium")

# Read in tree data
Tree1 <- read.tree("~/Poryzae.support")

Tree1 <- drop.tip(Tree1, tip = "87-120")

# Read in metadata
tipdata <- read.table("~/PoryzaeMetadata.txt")
colnames(tipdata) <- c("strain", "hostID", "cladeID")
tipdata$strain = factor(tipdata$strain)
tipdata$hostID <- factor(tipdata$hostID)
tipdata$cladeID <- factor(tipdata$cladeID, levels = cladeOrder)

# create genus labels
genusLabels <- paste0("*", names(hostcolors), "*")

# generate metadata for fruit
fruitdata <- tipdata
fruitdata$clade <- fruitdata$cladeID

# List of tip labels to show in plot1
selected_tips <- paste(c("T6", "T21", "T29", "U247", "U248", "U249", "EiSM270", "EiSM276", "UbJA108", "UbJA110", "UbJA118", "UbJA119", "UbJA171", "UFVPY113", "UFVPY63", "UFVPY184", "UFVPY231", "UFVPY218", "UbJA121", "UbJA158", "UFVPY166", "UFVPY198", "UbJA174"), collapse = "|")


# Plot the tree
p1 <- ggtree(Tree1, layout = "rectangular")

p2 <- p1 %<+% tipdata + 
  geom_tiplab(aes(subset = str_detect(label, selected_tips)), size = 2.5, align = T, line.size = 0.25, offset = 0.01, show.legend = F) +
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) ==100), size = 1.5, show.legend=F) +
  geom_tippoint(aes(color = hostID), stroke = 0, alpha=0.6, size = 3, position = position_nudge(x = 0.005)) +
  scale_color_manual(values=hostcolors,
                     name ="Host Species",
                     labels = genusLabels,
                     guide=guide_legend(keywidth=0.8,
                                        keyheight=0.8,
                                        ncol=1,
                                        order=2,
                                        override.aes=list(size=4,alpha=0.6)))

p3 <- p2 +
  geom_fruit(data=fruitdata, geom=geom_tile, mapping = aes(y=strain, fill=clade), width = 0.03, offset = 0.25, alpha = 0.6) +
  scale_fill_manual(values=cladecolors,
                    name = "Clade ID",
                    labels = cladeOrder,
                    na.translate = F,
                    guide=guide_legend(keywidth=0.7,
                                       keyheight=0.4,
                                       ncol=1,
                                       order=1
                    )) +

  theme(plot.margin = unit(c(4, 0, 4, 0), "cm"), 
        legend.position = "right", 
        legend.box = "vertical", 
        legend.direction = "vertical",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0.5),
        legend.spacing.y = unit(18, "cm"),
 #       legend.margin = margin(r = 20, l = -10),
        title = element_text(size = rel(1.5)), 
        legend.text = element_markdown(size = rel(1)),
        legend.key.size = unit(0.8, "cm"))

p4 <- p3 #+ coord_fixed(ratio = 2)

# define species for plot2
species <- c("P.grisea", "P.oryzae", "P.penniseticola", "P.pennisetigena", "P.urashimae", "P.zingibericola", "Ps.cyperi", "Ps.javanica")

# Read tree
Tree2 <- read.tree("~/Pyricularia.support")

# Read in metadata
data <- read.table("~/PyriculariaMetadata.txt", header =T, sep = "", quote="",  comment.char = "", stringsAsFactors = FALSE)
speciesdata <- data %>% select(c("strain", "species", "color"))
speciesdata$species <- factor(speciesdata$species, levels = species)
                            
# define colors and legend labels
SPECIEScolors <- speciesdata[match(species,data$species), "color"]
SPECIESlabels <- levels(factor(speciesdata$species))
SPECIESlabels <- gsub("\\.", "\\. ", SPECIESlabels)

# Add asterisks to allow rendering of italics
SPECIESlabels <- paste0("*", SPECIESlabels, "*")

# select tips for highlighting in tree
new_selected_tips <- c("BP1", "UbJA125", "CD86", "UFVPY202", "UFVPY204", "UFVPY210", "UFVPY232", "UFVPY657", "UFVPY677", "RN1", "Cr9010", "CrA8401")

# reorder species data to match tree tips
speciesdata <- speciesdata[match(Tree2$tip.label, speciesdata$strain), ]

# specify sizes for selected and unselected tip labels
speciesdata$sizes <- ifelse(Tree2$tip.label %in% new_selected_tips, 5, 3)


# Plot the tree
p5 <- ggtree(Tree2, layout = "rectangular")
  
p6 <- p5 %<+% speciesdata +
  geom_tippoint(aes(color=species),
                alpha=0) +
  geom_nodepoint(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) ==100), size = 1.5, show.legend=F) +
  geom_tiplab(aes(color = species, size = sizes), align = T, linesize = 0.25, show.legend = F) +
  scale_size_continuous(range=c(3,5)) +
  scale_color_manual(name = "Species",
                     values=SPECIEScolors,
                     labels=SPECIESlabels,
                     guide=guide_legend(keywidth=0.8,
                                        keyheight=0.8,
                                        order=1,
                                        override.aes=list(size=4, alpha=1)) 
  )


#geom_nodelab(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 90), size = 4, hjust = -0.5, vjust = 2, show.legend=F) +
  #scale_color_manual(values=speciescolors) +  
p7 <- p6  + theme(plot.margin = unit(c(4, 0, 4, 0), "cm"), 
                 legend.position = c(0.85, 0.93),
                 legend.box = "vertical", 
                 legend.margin = margin(t = 0, r = 2, b = 0, l = 0),
                 legend.box.margin = margin(r = 2, l = -1),  
                 legend.direction = "vertical", 
                 title = element_text(size = rel(1.5)), 
                 legend.text = element_markdown(size = rel(1)),
                 legend.key.size = unit(0.8, "cm")) + ggplot2::xlim(0, 1.6)

p8 <- p7 #+ coord_fixed(ratio = 2)

# Arrange plots with labels
grid.arrange(p4, p8, ncol=2, widths = c(5,3))

pdf("Meyer_et_al_Fig1.pdf", 14, 18)
grid.arrange(p4, p8, ncol=2, widths = c(4,3))
dev.off()


