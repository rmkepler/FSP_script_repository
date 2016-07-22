

setwd("~/Documents/FSP_metagenomics/nematode_publication_analyses/Nematode_backbone/constraint")

library("ggtree", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
library(phylotools)
library("colorspace")

#### import original CIPRES "bipartitions" file ####
#### this is the tree that will have names changed to long form ####
#### import taxonomy table ####
#### make tree specific data frames and change tip labels to full length ####

bb_taxonomy<-read.csv("merged_constraint_taxonomy_update.csv", stringsAsFactors=FALSE, header = T, fill =T, na.strings=c(""))
bb_taxonomy[is.na(bb_taxonomy)] <- "other"

bb_tree<-read.tree(file = "RAxML_bipartitions.constraint_0608")
leaves<-bb_tree$tip.label
in_tree<-NULL
for (name in leaves) {
  in_tree<-rbind(in_tree, bb_taxonomy[grep(paste(name, "\\b", sep = ""), bb_taxonomy$bb_name),])
  in_tree_key<-in_tree[,c("bb_name", "key_name")]
}

annot_keys<-in_tree[,c("key_name", "family")]

fullname_tree<-sub.tip.label(bb_tree, in_tree_key)

fullname_root<-root(fullname_tree, node = 5293, resolve.root = FALSE)

#row.names(bb_taxonomy) <- NULL

#### output a tree with tippoint colored by family  ####

pdf(file = "FSP_only/constraint_tree_circle_tippoint_fig1.pdf", width = 10, height = 20)
p <- ggtree(fullname_root, layout = "circular", branch.length = "none", size = 0.1)

p <- p %<+% annot_keys + geom_tippoint(aes(color = family), size = 1, shape = 25) + 
  ggtitle ("Constraint") +
  theme(plot.title=element_text(size=40),
        legend.key = element_rect(colour = "white"), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.just = "left")
p
dev.off()

#### output a tree with branches colored by family, all together ####

#original grouping loop
fam_list<-list()
for (var in unique(in_tree$family)) {
  temp1<-in_tree[grep(paste(var, "\\b", sep = ""), in_tree$family),]
  fam_list[var]<-list(var = c(temp1$key_name))
}

#make a FSP specific data frame and get family groupings
FSP<-in_tree[grep(paste("FSP", "\\b", sep = ""), in_tree$source),]

fsp_list<-unique(FSP$family)
fsp_list<-sort(fsp_list)
fam_list1<-list()
for (var in fsp_list) {
  temp1<-FSP[grep(paste(var, "\\b", sep = ""), FSP$family),]
  fam_list1[var]<-list(var = c(temp1$key_name))
}

#Add reference group to list
ref_tips<-list()
test<-"reference"
for (var in test) {
  temp2<-FSP[grep(paste(var, "\\b", sep = ""), FSP$family),]
  ref_tips<-temp2$key_name
}
fam_list1<-append(fam_list1, list(reference = ref_tips))

#make the trees

group_tree<-groupOTU(fullname_root, fam_list1)

pdf(file = "FSP_only/constraint_FSP_only_families_v2.pdf", width = 8, height = 8)
bp <- ggtree(group_tree, layout = "circular", branch.length = "none", size = 0.25, aes(color=group)) +
  scale_color_manual(name = "Families", values=fam_colors) + 
  ggtitle ("Constraint: No Branch Lengths-FSP families") +
  theme(plot.title=element_text(size=20) ,
        legend.key = element_rect(colour = "white"), 
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.just = "left")
bp
dev.off()

####
pdf(file = "FSP_only/constraint_FSP_only_families_fig1_rect_v2.pdf", width = 8, height = 50)
bp <- ggtree(group_tree, size = 0.25, aes(color=group)) +
  scale_color_manual(name = "Families", values=fam_colors) + 
  geom_text2(aes(label=label, subset=isTip), fontface = 3, size=0.25, hjust = -0.04) +
  ggtitle ("Constraint- Branch Lengths") +
  theme(plot.title=element_text(size=20),
        #plot.margin = unit(c(1, 0.5, 0.5, 0.5), "in"),
        legend.position = "bottom")
bp
dev.off()

####

#### Batch output series of trees, colored with a single family ####

uni_fams<-unique(annot_keys$family)
sm_fams<- c("Cephalobidae", "Rhabditidae", "Dorylaimidae", "Tylenchidae", "Qudsianematidae", "Plectidae", "Aphelenchidae", "Other", "Dolichodoridae" )

tip_table<-in_tree[,c("key_name", "source")]


#DOUBLE CHECK TO MAKE SURE THIS IS THE LIST YOU WANT
for (var in fsp_list) {
  temp<-annot_keys[grep(paste(var, "\\b", sep = ""), annot_keys$family),]
  group<-assign(var, c(temp$key_name))
  test_fam_tree<-groupOTU(fullname_root, group)
  p<-ggtree(test_fam_tree, layout = "circular", branch.length = "none", size = 0.25, aes(color=group)) +
    #geom_text2(aes(subset=!isTip, label=node), size = 0.25, hjust=-.3) +
    #geom_text2(aes(label=label, subset=isTip), fontface = 3, size=0.25, hjust = -0.1) +
    scale_color_manual(values = c("black", "firebrick")) +
    ggtitle(var)
  #p<- p %<+% tip_table +
    #geom_tippoint(aes(shape = source), size = 0.25, alpha = 0.80) +
    #theme(panel.background = element_rect(colour = "pink")) +
    #theme(legend.key = element_rect(colour = "white"), 
          #legend.position = "left",
          #legend.direction = "vertical",
          #legend.box.just = "left") +
    #theme(plot.margin = unit(c(0,0,0,0), "in"))
  
  #pdf(file = paste("indiv_fam_trees/Nematoda_uni_fams_", var, ".pdf", sep ="") , width = 8, height = 50)
  pdf(file = paste("FSP_only/circle_FSP_uni_fams_", var, ".pdf", sep ="") , width = 8, height = 8)
  print(p)
  dev.off()
}

#### color tree by individual families, manually one at a time ####

uni_fams<-unique(annot_keys$family)
for (var in uni_fams) {
  temp<-annot_keys[grep(var,annot_keys$family),]
  assign(var, c(temp$key_name))
}
fam_tree<-groupOTU(fullname_root, Rhabditidae)
pdf(file = "FSP_only/Nematoda_rhabditidae.pdf" , width = 8, height = 50)
ggtree(fam_tree, size = 0.25, aes(color=group)) +
  geom_text2(aes(label=label, subset=isTip), fontface = 3, size=0.25, hjust = -0.04) +
  #pdf(file = "indiv_fam_trees/Nematoda_Rhabditidae_color.pdf", width = 10, height = 10)
  #ggtree(fam_tree, layout = "circular", aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  ggtitle("Rhabditidae")
dev.off()
plot(test_fam_tree)
#### original raxml tree from CIPRES, colored by bootstrap value, node labels ####

raxml <- read.raxml("RAxML_bipartitionsBranchLabels.backbone")

pdf(file = "oriname_backbone_tree.pdf", width = 15, height = 30)
raxml_p <- ggtree(raxml, aes(color = bootstrap)) + 
  geom_text2(aes(subset=!isTip, label=node), size = 1.5, hjust=-.3) + 
  geom_text2(aes(label=label, subset=isTip), fontface = 3, size=2.5, hjust = -0.04)+
  scale_color_continuous(name ='bootstrap', limits = c(70, 100), low = "red", high = "firebrick") +
  ggtitle ("Backbone bootstrap")
plot(raxml_p)  
dev.off()

#####

##### REMOVE #####

remove(p1)
rm(list = unique(annot_keys$family))
remove(in_tree, in_tree_key, temp1, var, name, tips, families, p, raxml, bb_tree, binom_root, binom_name_tree, bb_taxonomy)

####  Big list of colors for FSP only families ####

fam_colors<- c(
  "black", # 0: Reference
  "lightpink1", # 1: Alaimidae
  "grey40", # 2: Anatonchidae
  "grey", # 3: Anguinidae
  "gold", # 4: Aphelenchidae
  "yellow", # 5: Aphelenchoididae
  "dodgerblue4", #6: Cephalobidae
  "lightcyan2", #7: Chromadoridae
  "green", #8: Dolichodoridae
  "darkorchid4", #9:Dorylaimidae
  "lavenderblush3", #10: Heteroderidae
  "lightpink3", #11: Hoplolaimidae
  "lightpink2", #12: Ironidae
  "lightsalmon", #13: Longidoridae
  "lightsalmon3", #14: Monhysteridae
  "lightsteelblue3", #15: Mononchidae
  "orange", #16: Nygolaimidae
  "palegreen3", #17: other
  "lightskyblue1", #18: Pangrolaimidae
  "mediumorchid4", #19: Plectidae
  "peachpuff3", #20: Pratylenchidae
  "deeppink1", #21: Prismatolaimidae
  "springgreen4", #22: Qudsianematidae
  "chartreuse4", #23: Rhabditidae
  "bisque3", #24: Steinernematidae
  "burlywood3", #25: Strongyloididae
  "darkolivegreen2", #26: Telotylenchidae
  "lemonchiffon2", #27: Teratocephalidae
  "honeydew4", #28: Thornenematidae
  "seashell3", #29: Tripylidae
  "firebrick", #30: Tylenchidae
  "rosybrown3", #31: Tylenchorhynchidae
  "steelblue3" #32: Xiphinematidae
)