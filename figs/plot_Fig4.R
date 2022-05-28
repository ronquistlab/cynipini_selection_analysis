library(ape)
library(cowplot)
library(tibble)
library(tidyr)
library(tidytree)
library(treeio)
library(ggplot2)
library(ggtree)

draw_tree <- function(tree, hexp, hypothesis="1")
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
    tree <- ape::root(tree, node=27, edgelabel=TRUE)
    tree$node.label[1]<-""
    tree <- ape::rotateConstr(tree, c(tipOrder))
    
    # Get indices for clades and branches we want to color
    tb <- as_tibble(tree)
    cynipini_core1 <- MRCA(tb, "Calli", "A_curv")$node
    cynipini_core2 <- MRCA(tb, "Belon", "A_curv")$node
    aulacideini <- MRCA(tb, "Hedic", "Aulac")$node
    s_itoensis <- tb$node[match("S_ito",tb$label)]
    neuroterus <- tb$node[match("Calli",tb$label)]
    protobalandricus <- tb$node[match("Proto",tb$label)]
    diastrophini <- tb$node[match("Diast",tb$label)]
    qwaqwaiini <- tb$node[match("Qwaqw",tb$label)]

    # Add the lineage info to the tree (using the "group" attribute)
    if (hypothesis=="1") {
        int_clades <- c(Foreground=cynipini_core2, Background=aulacideini)
        tree <- groupClade(tree, int_clades)
        x <- as.character(attr(tree,"group"))
        for(i in 1:length(x))
            if (x[i]=="0") x[i]<-"Background"
        x[neuroterus] <- "Foreground"
        x[cynipini_core1] <- "0"
    } else {
        int_clades <- c(Foreground=cynipini_core1, Background=aulacideini)
        tree <- groupClade(tree, int_clades)
        x <- as.character(attr(tree,"group"))
        for(i in 1:length(x))
            if (x[i]=="0") x[i]<-"0"
        x[protobalandricus] <- "Background"
        x[diastrophini] <- "Background"
        x[s_itoensis] <- "Background"
    }
    attr(tree,"group") <- as.factor(x)

    # Set the colors
    cols <- c(Foreground="green", Background="blue2")
    
    # Change display names and only show support < 100% 
    for ( i in 1:length(tree$tip.label) )
        tree$tip.label[i] <- displayNames[ match(tree$tip.label[i],taxonNames) ]
    for ( i in 1:length(tree$node.label) ) {
        if ( tree$node.label[i] == "100" )
            tree$node.label[i] <- ""
    }

    # Position the legend
    leg.pos <- c(.21,.90)

    # Position bars and size text
    ggtree(tree, aes(color = group), ladderize = FALSE) + 
        geom_tree(size=0.8) +
        geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2.0, colour="black", offset=0.01) +
        guides(color = guide_legend(override.aes = list(size = 1.5, shape = 15))) +
        theme_tree(legend.position = leg.pos, legend.key.size = unit(0.06,'cm'), legend.title = element_text(size=9)) +
        scale_color_manual(values = c(cols, "black"), na.value = "black", name = "Evolutionary regime", 
            breaks = c("Foreground", "Background")) +
        hexpand(hexp) 
}

t1 <- read.tree("../iqtree/A-50-18_KOSI07-IG/A-50-18codGIb.contree")

p1 <- draw_tree(t1, 0.28)
p2 <- draw_tree(t1, 0.28, "2")

labels <- c("A", "B")
cowplot::plot_grid(p1,p2,ncol=2,labels=labels, label_size=12, hjust=0.0) + theme(plot.margin=unit(c(3,3,3,3), "pt"))

ggsave("Fig_4.png", device="png", height=3.5)

