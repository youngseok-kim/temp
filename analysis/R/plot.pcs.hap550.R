# Create a nice plot showing the hap550 samples projected onto PCs 1 and 2.
library(hdf5r)
library(ggplot2)
library(cowplot)

# Load the PCA results.
out   <- H5File$new("../../data/hap550_pc.mat","r")
study <- out[["study"]][,]
pc    <- out[["pc"]][,]

# Plot the samples projected onto PCs 1 and 2.
dat        <- as.data.frame(cbind(factor(study),pc))
names(dat) <- c("study",paste0("PC",1:10))
ggplot(dat,aes(x = PC1,y = PC2,color = study,shape = study)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("dodgerblue","darkorange","forestgreen")) +
  scale_shape_manual(values = c(19,17,8))

