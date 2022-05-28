library(ape)
library(cowplot)
library(tibble)
library(tidyr)
library(tidytree)
library(treeio)
library(ggplot2)
library(ggtree)

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


draw_tree <- function(tree, hexp, PSG=TRUE)
{
    require(ggtree)
    require(ape)

    # Get indices for clades and branches we want to color
    tb <- as_tibble(tree)
    cynipini_core <- MRCA(tb, "Calli", "A_curv")$node
    synergini <- MRCA(tb, "S_ito", "S_jap")$node
    ceroptresini <- tb$node[match("Cerop",tb$label)]
    s_itoensis <- tb$node[match("S_ito",tb$label)]

    # Add the lineage info to the tree (using the "group" attribute)
    int_clades <- c(Complex=cynipini_core, Inquiline=synergini)
    tree <- groupClade(tree, int_clades)
    x <- as.character(attr(tree,"group"))
    x[ceroptresini] = "Inquiline"
    for (i in 1:length(x)) if (x[i]=="0") x[i] <- "Simple"
    x[s_itoensis] <- "Simple"
    attr(tree,"group") <- as.factor(x)

    index <- c(1, 20, 2, 21, 23, 3, 4, 22, 24, 25, 26, 27, 28, 29, 30, 5, 6, 7, 8, 31, 9, 32, 33, 10, 34, 11, 12, 13, 14, 15, 16, 17, 18, 19)
    rev_index <- 1 + c(0, 2, 5, 6, 15, 16, 17, 18, 20, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 1, 3, 7, 4, 8, 9, 10, 11, 12, 13, 14, 19, 21, 22, 24)
    psg <- c(98, 40, 150, 28, 90, 165, 133, 114, 78, 44, 28, 88, 236, 10, 61, 107, 119, 142, 156, 8, 209, 306, 47, 136, 182, 88, 82, 118, 244, 301, 105, 230, 114, NA)
    nsg <- c(44, 7, 134, 0, 98, 183, 168, 84, 18, 9, 1, 15, 189, 0, 15, 74, 112, 111, 142, 0, 279, 393, 8, 70, 107, 31, 18, 53, 323, 371, 152, 290, 96, NA)
    
    psg <- psg[rev_index]
    nsg <- nsg[rev_index]

    if (PSG==TRUE)
        attr(tree,"N") <- as.factor(psg)
    else
        attr(tree,"N") <- as.factor(nsg)
    
    # Set the colors
    cols <- c(Complex="green", Simple="blue2", Inquiline="skyblue")

    # Change display names and only show support < 100% 
    for ( i in 1:length(tree$tip.label) )
        tree$tip.label[i] <- displayNames[ match(tree$tip.label[i],taxonNames) ]

    # Position the legend
    leg.pos <- c(.21,.90)

    # Position bars and size text
    ggtree(tree, aes(color=group), ladderize = FALSE) + 
        geom_tree(size=0.8) +
        geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2.0, colour="black", offset=0.01) +
        geom_label(aes(x=branch, label=N), size=2.5, color="black", label.size=0.2, label.padding=unit(0.05,'cm'), show.legend=FALSE, na.rm=TRUE) + 
        guides(color = guide_legend(override.aes = list(size = 2.5, shape = 15))) +
        theme_tree(legend.position = leg.pos, legend.key.size = unit(0.05,'cm'), legend.title = element_text(size=9)) +
        scale_color_manual(values = c(cols, "black"), na.value = "black", name = "Gall type", 
            breaks = c("Complex", "Simple", "Inquiline")) +
        hexpand(hexp) 
}

t1 <- read.tree("../iqtree/A-50-18_KOSI07-IG/A-50-18codGIb.contree")
t1 <- ape::root(t1, node=27, edgelabel=TRUE)
t1 <- ape::rotateConstr(t1, c(tipOrder))
t2 <- t1

p1 <- draw_tree(t1, 0.28)
p2 <- draw_tree(t2, 0.28, FALSE)

labels <- c("A", "B")
cowplot::plot_grid(p1,p2,ncol=2,labels=labels, label_size=12, hjust=0.0) + theme(plot.margin=unit(c(3,3,3,3), "pt"))

ggsave("Fig_3.png", device="png", height=3.5)

