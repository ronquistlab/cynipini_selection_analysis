library(ape)
library(cowplot)
library(tibble)
library(tidyr)
library(tidytree)
library(treeio)
library(ggplot2)
library(ggtree)

draw_tree <- function(tree, hexp, rootNode)
{
    require(ggtree)
    require(ape)

    # Taxon names actually used in the trees
    taxonNames <- c(
    "A_curv","A_gros","A_qlan","Belon","Aulac","Isoco",
    "Aylax","Hedic","Qwaqw","S_gif","S_jap", "S_umb",
    "S_ito","Diast","Cerop","Proto","Calli","A_qram"
    )

    # Choose the corresponding display names we want to use
    displayNames <- c(
    "Andricus cur","Andricus gro","Andricus qln","Belonocnema kin","Aulacidea tav","Isocolus cen",
    "\"Aylax\" hyp","Hedickiana lev","Qwaqwaia sco","Synergus gif","Synergus jap", "Synergus umb",
    "Synergus ito","Diastrophus kin","Ceroptres mas","Protobalandricus spe","Neuroterus val","Andricus qrm"
    )

    # Order tip labels the way we want them
    tipOrder <- c(
    "A_curv","A_qram","A_gros","A_qlan","Belon","Calli",
    "Proto","Cerop","Diast","Aulac","Isoco","Aylax",
    "Hedic","Qwaqw","S_ito","S_gif","S_umb","S_jap"
    )
    tipOrder <- tipOrder[18:1]

    # Root the tree correctly
    tree <- ape::root(tree, node=rootNode, edgelabel=TRUE)
    tree$node.label[1]<-""
    tree <- ape::rotateConstr(tree, c(tipOrder))
    
    # Change display names and only show support < 100% 
    for ( i in 1:length(tree$tip.label) )
        tree$tip.label[i] <- displayNames[ match(tree$tip.label[i],taxonNames) ]
    for ( i in 1:length(tree$node.label) ) {
        if ( tree$node.label[i] == "100" )
            tree$node.label[i] <- ""
    }

    # Position bars and size text
    ggtree(tree, ladderize = FALSE) + 
        geom_nodelab(size=2.5,hjust=0, nudge_x=0.001) +
        geom_tree(size=0.8) +
        geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2.5) +
        geom_treescale(x = 0.0, y = 0.0, width = 0.02, linesize=0.5, fontsize=2.5) + 
        hexpand(hexp) 
}

t1 <- read.tree("../iqtree/A-50-18_C60/C60-A275-50-18-1.contree")
t2 <- read.tree("../iqtree/A-50-18_C60/C60-A275-50-18-2.contree")
t3 <- read.tree("../iqtree/A-50-18_C60/C60-A275-50-18-3.contree")
t4 <- read.tree("../iqtree/A-50-18_C60/C60-A275-50-18-4.contree")

p1 <- draw_tree(t1, 0.23, 27)
p2 <- draw_tree(t2, 0.23, 27)
p3 <- draw_tree(t3, 0.23, 27)
p4 <- draw_tree(t4, 0.23, 27)

labels <- c("A", "B", "C", "D")
cowplot::plot_grid(p1,p2,p3,p4,ncol=2,labels=labels, label_size=12, hjust=0.0) + theme(plot.margin=unit(c(3,3,3,3), "pt"))

ggsave("Fig_2.png", device="png")
