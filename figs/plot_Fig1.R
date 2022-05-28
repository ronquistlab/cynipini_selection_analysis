library(ape)
library(cowplot)
library(tibble)
library(tidyr)
library(tidytree)
library(treeio)
library(ggplot2)
library(ggtree)

draw_tree <- function(tree, hexp)
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

    # Set the colors
    cols <- c(Complex="green", Simple="blue2", Inquiline="skyblue")
    
    # Change display names and only show support < 100% 
    for ( i in 1:length(tree$tip.label) )
        tree$tip.label[i] <- displayNames[ match(tree$tip.label[i],taxonNames) ]
    for ( i in 1:length(tree$node.label) ) {
        if ( tree$node.label[i] == "100" )
            tree$node.label[i] <- ""
    }

    # Position the legend
    leg.pos <- c(.1,.90)

    # Position bars and size text
    font_size <- 4.0    # default is 3.88
    cyn_offset <- 0.14
    cyn_offset_text <- 0.02
    ggtree(tree, aes(color = group), ladderize = FALSE) + 
        geom_nodelab(size=2,hjust=0) +
        geom_tree(size=1.5) +
        geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 3.0, colour="black", offset=0.01) +
        geom_treescale(x = 0.0, y = 0.0, width = 0.1, linesize=0.8, fontsize=3) + 
        guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
        theme_tree(legend.position = leg.pos, legend.key.size = unit(0.08,'cm'), legend.title = element_text(size=11)) +
        scale_color_manual(values = c(cols, "black"), na.value = "black", name = "Gall type", 
            breaks = c("Complex", "Simple", "Inquiline")) +
        geom_strip('Andricus cur', 'Protobalandricus spe', barsize=1, fontsize=font_size, color='black', label="Cynipini", offset=cyn_offset, offset.text=cyn_offset_text) +
        geom_strip('Ceroptres mas', 'Ceroptres mas', barsize=1, fontsize=font_size, color='black', label="Ceroptresini", offset=cyn_offset-0.22, offset.text=cyn_offset_text) +
        geom_strip('Diastrophus kin', 'Diastrophus kin', barsize=1, fontsize=font_size, color='black', label="Diastrophini", offset=cyn_offset-0.28, offset.text=cyn_offset_text) +
        geom_strip('Aulacidea tav', 'Hedickiana lev', barsize=1, fontsize=font_size, color='black', label="Aulacideini", offset=cyn_offset-0.25, offset.text=cyn_offset_text) +
        geom_strip('Synergus ito', 'Synergus jap', barsize=1, fontsize=font_size, color='black', label="Synergini", offset=cyn_offset-0.15, offset.text=cyn_offset_text) +
        geom_strip('Qwaqwaia sco', 'Qwaqwaia sco', barsize=1, fontsize=font_size, color='black', label="Qwaqwaiini", offset=cyn_offset-0.24, offset.text=cyn_offset_text) +
        hexpand(hexp) 
}

t1 <- read.tree("../iqtree/A-50-18_KOSI07-IG/A-50-18codGIb.contree")

p1 <- draw_tree(t1, 0.10)

p1

ggsave("Fig_1.png", device="png", height=7.0)

