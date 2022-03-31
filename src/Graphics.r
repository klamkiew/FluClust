suppressPackageStartupMessages({
#https://yulab-smu.top/treedata-book/chapter6.html 
require('ggpubr')
require('ggplot2')# || install.packages('ggplot2', dependencies = TRUE)
#require('BiocManager')# || install.packages('BiocManager', dependencies = TRUE)
require('ggtree')# || BiocManager::install('ggtree')
require('ggtreeExtra')# || BiocManager::install('ggtreeExtra')
require('treeio')
require('tidytree')
require('tidyverse')
require('dplyr')
require('ggstar')
#require('TDbook')
require('stringr')
require('ape')
require('plotly')
#require('viridis')
require('ggnewscale')
require('tanggle')
require('phangorn')
require('svglite')
})

name <- 'influenza_segment_6'
outfolder <- 'Result_Segment_6'
color <- 'N'

if (!substr(outfolder,nchar(outfolder),nchar(outfolder)) == '/'){
    outfolder <- paste0(outfolder, '/')
}

framecentroid <- read.csv(paste0(outfolder, 'centroids_', name, '.csv'), header = TRUE, sep = ',')
framemeta <- read.csv(paste0(outfolder, 'meta_', name, '.csv'), header = TRUE, sep = ',')
framecluster <- read.csv(paste0(outfolder, 'cluster_', name, '.csv'), header = TRUE, sep = ',')

framevectors <- read.csv(paste0(outfolder, 'vectors_', name, '.csv'), header = TRUE, sep = ',')
#, row.names = 1

cluster <- inner_join(framecluster, framemeta, by = 'accession')
vectors <- inner_join(framevectors, cluster, by = 'accession')
centroid <- inner_join(framecentroid, framemeta, by = 'accession')
framevectorsHD <- read.csv(paste0(outfolder, 'vectors_HD_', name, '.csv'), header = TRUE, sep = ',')
vectorsHD <- inner_join(framevectorsHD, cluster, by = 'accession')
pca <- read.csv(paste0(outfolder, 'pca_', name, '.csv'), header = TRUE, sep = ',')
info <- read.csv(paste0(outfolder, 'info_', name, '.csv'), header = TRUE, sep = ',')

#tree <- read.tree(paste0('centroids_', name, '.fasta.tree'))
tree <- read.tree(paste0(outfolder, 'RAxML_bestTree.', name, '.tree'))

listcolor <- c(
    "dodgerblue2", "#E31A1C",
    "green4",
    "#6A3D9A",
    "#FF7F00",
    "aquamarine4", "gold1",
    "skyblue2", "#FB9A99",
    "palegreen2",
    "#CAB2D6",
    "#FDBF6F",
    "gray70", "khaki2",
    "maroon", "orchid1", 
    "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
)

if(color != ''){

    count <- plyr::count(cluster[[color]])
    names(count)[names(count) == 'x'] <- color
    names(count) <- c(color, 'count')
    count[[color]] <- as.character(count[[color]])
    #count[[color]] <- count[[color]] %>% replace_na('mixed')

    subtype <- c(str_sort(levels(as.factor(count[[color]])), numeric = TRUE), 'mixed')
    count[[color]][is.na(count[[color]])] <- 'mixed'

    count[[color]] <- factor(count[[color]], levels=subtype)

    options(repr.plot.width=20, repr.plot.height=8)
    p8 <- ggplot(data=count, aes_string(x=color, y='count', label = 'count')) + 
        geom_bar(stat='identity') +
        theme_classic() +
        geom_col() +
        geom_text(nudge_y = 500, size = 5) + 
        xlab("Subtype") +
        ylab("Count") +
        theme(
            text = element_text(size=20),
            axis.text.x = element_text(size=15),
            axis.line=element_line(size=1),
            axis.ticks=element_line(size=1)
        )

    print(p8)
    ggsave(paste0(outfolder, 'overview_', name, '.pdf'), width = 50, height = 20, units = "cm", limitsize = FALSE)
    knitr::plot_crop(paste0(outfolder, 'overview_', name, '.pdf'))
}

if(color != ''){
    #frame <- centroid[ , c('accession', color)] %>% remove_rownames %>% column_to_rownames(var="accession") 
    subtype <- c(str_sort(levels(as.factor(centroid[[color]])), numeric = TRUE), 'mixed')

    centroid[[color]][is.na(centroid[[color]])] <- 'mixed'
    centroid[[color]] <- factor(centroid[[color]], levels=subtype)
    frame <- centroid[ , c('accession', color)] %>% remove_rownames %>% column_to_rownames(var="accession") 

    q <- length(subtype)-1
}

options(repr.plot.width=20, repr.plot.height=20)

p0 <- ggtree(tree, layout='circular') %<+% centroid

p1 <- open_tree(p0, 0) +
    geom_treescale(width = 0.25) + 
    geom_tippoint(aes(size=size), alpha=.6) +
    geom_tiplab2(aes(label=paste0(cluster, ':', label)), align=TRUE, offset = 0) +#, parse=T) +
    scale_size_continuous(
        name = 'size',
        breaks = c(2, 10, 100, 1000, 10000), 
        guide = guide_legend(
            title.theme = element_text(size = 20, face = "bold"),
            label.theme = element_text(size = 15)
        )
    )   

if(color != ''){
    p2 <- p1 +
         new_scale_fill()

    p3 <- gheatmap(
            p2, 
            frame, 
            offset=0.75, 
            width=0.05, 
            #font.size=3, 
            #colnames_angle=0, 
            colnames = TRUE, 
            #hjust=0
        ) +
        {if(color != '')scale_fill_manual(
            name = color,
            values = c(listcolor[0:q], 'black'),
            breaks = subtype, 
            labels = subtype,
            #na.value = "black",
            guide = guide_legend(
                override.aes = list(shape = 15, size = 10),
                title.theme = element_text(size = 20, face = "bold"),
                label.theme = element_text(size = 15)
            )
        )}
}else{
    p3 <- p1
}

print(p3)
ggsave(paste0(outfolder, 'centroid_', name, '.svg'), width = 50, height = 50, units = "cm", limitsize = FALSE)
ggsave(paste0(outfolder, 'centroid_', name, '.pdf'), width = 50, height = 50, units = "cm", limitsize = FALSE)
#knitr::plot_crop(paste0(outfolder, 'centroid_', name, '.pdf'))

msa <- read.dna(paste0(outfolder, 'centroids_', name, '.msa.fasta'), format = "fasta")
dm <- dist.ml(msa)
nnet <- neighborNet(dm)
#nnet <- read.nexus.networx(paste0('centroids_', name, '.nex'))

if(color != ''){
    k <- sapply(nnet$tip.label, (\(x) centroid[which(centroid$accession == x),][[color]]), simplify = TRUE, USE.NAMES = TRUE)
}

c <- sapply(nnet$tip.label, (\(x) centroid[which(centroid$accession == x),]$cluster), simplify = TRUE, USE.NAMES = TRUE)
s <- sapply(nnet$tip.label, (\(x) centroid[which(centroid$accession == x),]$size), simplify = TRUE, USE.NAMES = TRUE)

options(repr.plot.width=20, repr.plot.height=20)

p4 <- ggsplitnet(nnet) 

p5 <- p4 + 
    {if(color != '')geom_tiplab2(aes(label = paste0(c[label], ':', label), color = k[label]), hjust = -.5)} +
    {if(color == '')geom_tiplab2(aes(label = paste0(c[label], ':', label)), hjust = -.5)} +
    {if(color != '')geom_tippoint(aes(color = k[label], size = s[label]))} +
    {if(color == '')geom_tippoint(aes(size = s[label]))} +
    theme(legend.position = "right") +
    {if(color != '')scale_color_manual(
        name = color,
        values = c(listcolor[0:q], 'black'),
        breaks = subtype, 
        labels = subtype,
        #na.value = "black",
        guide = guide_legend(
            override.aes = list(shape = 15, size = 10),
            title.theme = element_text(size = 20, face = "bold"),
            label.theme = element_text(size = 15)
        )
    )} +
    scale_size_continuous(
        name = 'size',
        breaks = c(2, 10, 100, 1000, 10000), 
        guide = guide_legend(
            #override.aes = list(shape = 15, size = 10),
            title.theme = element_text(size = 20, face = "bold"),
            label.theme = element_text(size = 15)
        )
    ) +
    ggexpand(.1) + ggexpand(.1, direction=-1)

print(p5)
ggsave(paste0(outfolder, 'nexus_', name, '.svg'), width = 50, height = 50, units = "cm", limitsize = FALSE)
ggsave(paste0(outfolder, 'nexus_', name, '.pdf'), width = 50, height = 50, units = "cm", limitsize = FALSE)
#knitr::plot_crop(paste0(outfolder, 'nexus_', name, '.pdf'))
if(color != ''){
    subtype <- c(str_sort(levels(as.factor(vectorsHD[[color]])), numeric = TRUE), 'mixed')
    vectorsHD[[color]][is.na(vectorsHD[[color]])] <- 'mixed'
    #vectors[[color]] <- ifelse(vectors$accession %in% centroid$accession, 'centroid', vectors[[color]])
    vectorsHD[[color]] <- factor(vectorsHD[[color]], levels=subtype)
    
    q <- length(subtype)-1
}

#vectorsHD$centroid <- with(vectorsHD, ifelse(accession %in% centroid$accession, 'C', 'nC'))
#vectorsHD$centroid <- factor(vectorsHD$centroid)#z <- length(vectorsHD)-length(colnames(cluster))+1
z <- 4+1options(repr.plot.width=20, repr.plot.height=20)
o <- 0
for(pcax in colnames(vectorsHD)[2:z]){
    for(pcay in colnames(vectorsHD)[2:z]){
        assign(
            paste0("o", o), ggplot(vectorsHD, aes_string(x = pcax, y = pcay, color = color)) +
            geom_point() + 
            theme_classic() +
            scale_color_manual(
                name = color,
                values = c(listcolor[0:q], 'black'),
                breaks = subtype,
                labels = subtype,
                guide = guide_legend(
                    override.aes = list(shape = 15, size = 10),
                    title.theme = element_text(size = 20, face = "bold"),
                    label.theme = element_text(size = 15)
                )
            ) + 
            xlim(c(-1, 1)) + 
            ylim(c(-1, 1)) +
            theme(
                text = element_text(size=20),
                axis.text.x = element_text(size=15),
                axis.line=element_line(size=1),
                axis.ticks=element_line(size=1)
            )
        )
        print(paste0("o", o))
        ggsave(paste0(outfolder, pcax, '_vs_', pcay, '_', name, '.png'), width = 50, height = 50, units = "cm", limitsize = FALSE)
        ggsave(paste0(outfolder, pcax, '_vs_', pcay, '_', name, '.pdf'), width = 50, height = 50, units = "cm", limitsize = FALSE)
        o <- o + 1
    }
}
#options(repr.plot.width=20, repr.plot.height=20)
#ggarrange(o1, o2, o3, o4, ncol=2, nrow=2, common.legend = TRUE, legend="right")
#ggsave(paste0(outfolder, 'components_', name, '.pdf'), width = 50, height = 50, units = "cm", limitsize = FALSE)
#knitr::plot_crop(paste0(outfolder, 'nexus_', name, '.pdf'))

if(color != ''){
    subtype <- c(str_sort(levels(as.factor(vectors[[color]])), numeric = TRUE), 'mixed')
    vectors[[color]][is.na(vectors[[color]])] <- 'mixed'
    #vectors[[color]] <- ifelse(vectors$accession %in% centroid$accession, 'centroid', vectors[[color]])
    vectors[[color]] <- factor(vectors[[color]], levels=subtype)
    
    q <- length(subtype)-1
}

vectors$centroid <- with(vectors, ifelse(accession %in% centroid$accession, 'C', 'nC'))
vectors$centroid <- factor(vectors$centroid)

#as.matrix(vectors[vectors$centroid == 'C', -1])
vec1 <- subset(vectors, centroid == 'nC')
vec2 <- subset(vectors, centroid == 'C')

rownames(vec1) <- NULL
rownames(vec2) <- NULL

options(repr.plot.width=20, repr.plot.height=20)

fig <- if(color != ''){
    plot_ly(
        vec1,
        x = ~PCA1,
        y = ~PCA2,
        z = ~PCA3,
        #symbol = ~centroid,
        #symbols = c('circle','x'),
        color = vec1[[color]], 
        colors = c(listcolor[0:q], 'black'),
        type = 'scatter3d', 
        mode = 'markers',
        #size = 1,
        marker = list(symbol = 'circle', size = 1),
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster),
        width = 1200,
        height = 1200
    )
} else {
    plot_ly(
        vec1,
        x = ~PCA1,
        y = ~PCA2,
        z = ~PCA3,
        type = 'scatter3d', 
        mode = 'markers',
        #size = 1,
        marker = list(symbol = 'circle', size = 1),
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster),
        width = 1200,
        height = 1200,
        name = 'mixed',
        color = I('black')
    )
}

fig <- if(color != ''){
    fig %>% add_trace(
        data = vec2,
        x = ~PCA1,
        y = ~PCA2,
        z = ~PCA3,
        #color = I('black'),
        color = vec2[[color]], 
        colors = c(listcolor[0:q], 'black'),
        type = 'scatter3d', 
        mode = 'markers',
        marker = list(symbol = 'x', size = 1),
        #name = 'centroid',
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster)
    )
} else {
    fig %>% add_trace(
        data = vec2,
        x = ~PCA1,
        y = ~PCA2,
        z = ~PCA3,
        color = I('black'),
        marker = list(symbol = 'x', size = 1),
        name = 'mixed',
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster)
    )
}

fig <- fig %>% layout(
    scene = list(
        xaxis = list(
            title = 'PCA1',
            #gridcolor = 'rgb(255, 255, 255)',
            #zerolinewidth = 1,
            #ticklen = 5,
            #gridwidth = 2,
            range = c(-1.0, 1.0)
        ),
        yaxis = list(title = 'PCA2',
            range = c(-1.0, 1.0)
        ),
        zaxis = list(title = 'PCA3',
            range = c(-1.0, 1.0)
        ),
        #paper_bgcolor = 'rgb(243, 243, 243)',
        #plot_bgcolor = 'rgb(243, 243, 243)',
        aspectmode = 'cube'
    ),
    autosize = TRUE,
    margin = c(l=0, r=0, b=0, t=0), 
    #title = 'Life Expectancy v. Per Capita GDP, 2007',
    legend = list(
        orientation = 'v',
        yanchor = 'middle',
        y = 0.5,
        itemsizing='constant'
    )
)

embed_notebook(fig)

htmltools::save_html(fig, paste0(outfolder, 'sphere_', name, '.html'))

reduced <- info$components

options(repr.plot.width=20, repr.plot.height=8)

p6 <- ggplot(data=pca, aes(x=components, y=variance)) +
    geom_line(size = 1) +
    geom_point() +
    theme_classic() +
    geom_vline(xintercept=c(3, reduced), linetype='dashed', color=c('blue', 'red'), size=1) +
    theme(
        text = element_text(size=20),
        axis.text.x = element_text(size=15),
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1)
    )
    
p7 <- p6 + annotate(x=c(3, info$components),y=+Inf,label=c(pca[pca$components==3,]$variance, pca[pca$components==reduced,]$variance),vjust=1,geom="label", size= 6)

print(p7)
ggsave(paste0(outfolder, 'pca_', name, '.pdf'), width = 50, height = 20, units = "cm", limitsize = FALSE)
knitr::plot_crop(paste0(outfolder, 'pca_', name, '.pdf'))
#as.matrix(vectors[vectors$centroid == 'C', -1])
vec1 <- subset(vectors, centroid == 'C')
vec2 <- subset(vectors, centroid == 'nC')rownames(vec1) <- NULL
rownames(vec2) <- NULLfig <- plot_ly(
    data = vec1,
    x = ~x,
    y = ~y,
    z = ~z,
    color = vec1[[color]], 
    colors = c(listcolor[0:q], 'black'),
    type = 'scatter3d', 
    mode = 'markers',
    marker = list(symbole = 'circle', size = 1),
    text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster),
    width = 1200,
    height = 1200,
    legendgroup = "G1"
)%>% add_trace(
    data = vec2,
    x = ~x,
    y = ~y,
    z = ~z,
    color = vec2[[color]], 
    colors = c(listcolor[0:q], 'black'),
    type = 'scatter3d', 
    mode = 'markers',
    marker = list(symbol = 'x', size = 1),
    text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster),
    legendgroup = "G2"
)%>% layout(
    scene = list(
        xaxis = list(
            title = 'PCA1',
            range = c(-1.0, 1.0)
        ),
        yaxis = list(
            title = 'PCA2',
            range = c(-1.0, 1.0)
        ),
        zaxis = list(
            title = 'PCA3',
            range = c(-1.0, 1.0)
        ),
        aspectmode = 'cube'
    ),
    autosize = TRUE,
    margin = c(l=0, r=0, b=0, t=0), 
    legend = list(
        orientation = 'h',
        itemsizing='constant',
        legend_traceorder='grouped'
    )
)

embed_notebook(fig)subtype <- c(str_sort(levels(as.factor(centroid[[color]])), numeric = TRUE), 'mixed')

    centroid[[color]][is.na(centroid[[color]])] <- 'mixed'
    centroid[[color]] <- factor(centroid[[color]], levels=subtype)str(vectors)#if(exists('color')){
if(color != ''){
    cluster[[color]] <- factor(cluster[[color]], levels=str_sort(levels(as.factor(cluster[[color]])), numeric = TRUE))
}j <- 1
i <- 0
h <- 0.75
z <- c()

for(feature in vec){

    #frame <- data.frame('feature' = centroid[,c(feature)])
    #rownames(frame) <- centroid$accession

    frame <- centroid[ , c('accession', feature)] %>% remove_rownames %>% column_to_rownames(var="accession")
    
    subtype <- str_sort(levels(as.factor(centroid[[feature]])), numeric = TRUE)
    #subtype <- c(str_sort(levels(as.factor(frame[[feature]])), numeric = TRUE), 'mixed')
    #z <- c(z, length(subtype)-1)
    z <- c(z, length(subtype))
    frame[[feature]][is.na(frame[[feature]])] <- 'mixed'
    frame[[feature]] <- factor(frame[[feature]], levels=subtype)
    
    #k <- i+length(subtype)-1
    k <- i+length(subtype)
    
    assign(
        paste0("p", 1+j), get(paste0("p", j)) +
        new_scale_fill()
    )

    assign(
        paste0("p", 2+j), gheatmap(
            get(paste0("p", 1+j)), 
            frame, 
            offset=h, 
            width=0.05, 
            #font.size=3, 
            #colnames_angle=0, 
            colnames = TRUE, 
            #hjust=0
        ) +
        scale_fill_manual(
            name = feature,
            values = c(listcolor[i:k], 'black'),
            breaks = subtype, 
            labels = subtype,
            #na.value = "black",
            guide = guide_legend(override.aes = list(shape = 15, size = 5))
        )
    )
    
    h <- h + 0.25
    #comment out i if now color change for additional rings
    i <- k + 1
    j <- j + 2
}get(paste0("p", j))
ggsave(paste0('centroid_', name, '.pdf'), width = 40, height = 40, units = "cm", limitsize = FALSE)#fig2 <- plot_ly(
#    pca, 
#    x = ~components, 
#    y = ~variance, 
#    type = 'scatter', 
#    mode = 'lines',
#    width = 1200, 
#    height = 300
#)

#embed_notebook(fig2)
#htmltools::save_html(fig2, paste0('pca_', name, '.html'))p <- ggtree(tree, layout='circular')
p <- p %<+% centroidp1 <- p + geom_tree() +  
    #xlim(-1, 1) +
    #aes(color=continent) + 
    geom_tippoint(aes(size=size), alpha=.6) + 
    #geom_tiplab(aes(label=country), offset=.1) +
    #theme(plot.margin=margin(60,60,60,60))
    geom_tiplab(aes(label=paste0(cluster, ':', label)), parse=T, size = 5, align=TRUE) +
    #geom_tiplab(size = 5, align=TRUE, linesize=.25) + 
    #geom_text(aes(label=label)) + 
    theme(
        legend.title=element_text(size=10), 
        legend.text=element_text(size=15),
        #legend.position = 'right', 
        legend.spacing.x = unit(1, 'cm'),
        #legend.margin = margin(t = 0, l = 0, b = -2.5, r = 0, unit='cm'),
        #        legend.key.height=unit(3,"line"),
        #        legend.key.width=unit(3,"line")
        #plot.margin = margin(t = 0, l = -5, b = -2.5, r = -5, unit='cm')
    )p <- geom_tree()  
gheatmap(p, centr) #+
#scale_x_ggtree() + 
#scale_y_continuous(expand=c(0, 0.3))p2 <- p1 +
    geom_fruit(
        geom=geom_tile,
        mapping=aes_string(fill='N'),
        #mapping = aes(fill=group),
        width=0.2,
        offset=0.5
    ) +
    #{if(exists('feature'))geom_fruit(
    #    geom=geom_tile,
    #    mapping=aes_string(fill='H'),
        #mapping = aes(fill=group),
    #    width=0.2,
    #    offset=0.1
    #)} +
    #{if(exists('feature'))scale_fill_manual(
    #    name="Subtype",
    #    values=listcolor,
    #    guide=guide_legend(
    #        keywidth=1.5,
    #        keyheight=1.5,
    #        order=3,
    #        nrow=3
    #    ),
    #    na.translate=FALSE
    #)} +
    #geom_fruit(
    #    geom=geom_bar,
    #    mapping=aes(
    #        y=accession, 
    #        x=size, 
    #        fill=group
    #    ),
    #    pwidth=0.5, 
    #    orientation="y", 
    #    stat="identity",
    #    offset=0.1,
    #) + 
    
    theme(
        legend.title=element_text(size=20), 
        legend.text=element_text(size=15),
        legend.position = 'top', 
        legend.spacing.x = unit(1, 'cm'),
        #legend.margin = margin(t = 0, l = 0, b = -2.5, r = 0, unit='cm'),
        #        legend.key.height=unit(3,"line"),
        #        legend.key.width=unit(3,"line")
        #plot.margin = margin(t = 0, l = -5, b = -2.5, r = -5, unit='cm')
    )p3 <- p2 +
    geom_fruit(
        geom=geom_tile,
        mapping=aes_string(fill='H'),
        #mapping = aes(fill=group),
        width=0.2,
        offset=0.5
    ) +
    #{if(exists('feature'))geom_fruit(
    #    geom=geom_tile,
    #    mapping=aes_string(fill='H'),
        #mapping = aes(fill=group),
    #    width=0.2,
    #    offset=0.1
    #)} +
    #{if(exists('feature'))scale_fill_manual(
    #    name="Subtype",
    #    values=listcolor,
    #    guide=guide_legend(
    #        keywidth=1.5,
    #        keyheight=1.5,
    #        order=3,
    #        nrow=3
    #    ),
    #    na.translate=FALSE
    #)} +
    #geom_fruit(
    #    geom=geom_bar,
    #    mapping=aes(
    #        y=accession, 
    #        x=size, 
    #        fill=group
    #    ),
    #    pwidth=0.5, 
    #    orientation="y", 
    #    stat="identity",
    #    offset=0.1,
    #) + 
    
    theme(
        legend.title=element_text(size=20), 
        legend.text=element_text(size=15),
        legend.position = 'top', 
        legend.spacing.x = unit(1, 'cm'),
        #legend.margin = margin(t = 0, l = 0, b = -2.5, r = 0, unit='cm'),
        #        legend.key.height=unit(3,"line"),
        #        legend.key.width=unit(3,"line")
        #plot.margin = margin(t = 0, l = -5, b = -2.5, r = -5, unit='cm')
    )plot(p3)
ggsave(paste0('centroids_', name, '.png'), width = 40, height = 40, units = "cm", limitsize = FALSE)

framevectors <- read.csv(paste0('vectors_', name, '.csv'), header = TRUE, sep = ',')
framecluster <- read.csv(paste0('cluster_', name, '.csv'), header = TRUE, sep = ',')#fig <- htmltools::div(fig, align="center", height = 1200, width = 1200) geom_facet(panel = "Trait", data = df_bar_data, geom = geom_col, 
                aes(x = dummy_bar_value, color = location, 
                fill = location), orientation = 'y', width = .6) +p <- p + new_scale_fill() +
         geom_fruit(data=dat2, geom=geom_tile,
                  mapping=aes(y=ID, x=Sites, alpha=Abundance, fill=Sites),
                  color = "grey50", offset = 0.04,size = 0.02)+
         scale_alpha_continuous(range=c(0, 1),
                             guide=guide_legend(keywidth = 0.3, 
                                             keyheight = 0.3, order=5)) +
         geom_fruit(data=dat3, geom=geom_bar,
                    mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
                    pwidth=0.38, 
                    orientation="y", 
                    stat="identity",
         ) +
         scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                                "#800000", "#006400","#800080","#696969"),
                           guide=guide_legend(keywidth = 0.3, 
                                        keyheight = 0.3, order=4))+
         geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
         theme(legend.position=c(0.93, 0.5),
               legend.background=element_rect(fill=NA),
               legend.title=element_text(size=6.5),
               legend.text=element_text(size=4.5),
               legend.spacing.y = unit(0.02, "cm"),
             )
pp2 <- plotly::ggplotly(p, height = 1000, width=600)
embed_notebook(p)tip <- as.phylo(beast_tree)$tip.label
beast_tree <- groupOTU(beast_tree, tip[grep("Swine", tip)], 
                       group_name = "host")

p <- ggtree(beast_tree, aes(color=host), mrsd="2013-01-01", 
            yscale = "label", yscale_mapping = NAG) + 
  theme_classic() + theme(legend.position='none') +
  scale_color_manual(values=c("blue", "red"), 
                     labels=c("human", "swine")) +
  ylab("Number of predicted N-linked glycoslyation sites")viewClade(p, MRCA(p, c("Merc 230", "Merc 280C")))tree <- read.tree("tree_newick.nwk")
treeset.seed(2015-12-21)
tree2 = rtree(30)
p <- ggtree(tree2) + xlim(NA, 6)
p + geom_cladelabel(node=45, label="test label", align=T, color='red') +
    geom_cladelabel(node=34, label="another clade", align=T, color='blue')library("ape")
#library("Biostrings")
library("ggplot2")
library("ggtree")
nwk <- system.file("extdata", "sample.nwk", package="ggtree")
tree <- read.tree(nwk)
tree


stdout, stderr = raxml_cline()
## R### Packages#### Graphics%load_ext rpy2.ipython%%R -i outfolder -i color -i name
require('ggplot2')
require('ggtree')
require('ggtreeExtra')
require('treeio')
require('tidytree')
require('tidyverse')
require('dplyr')
require('ggstar')
require('stringr')
require('ape')
require('plotly')
require('ggnewscale')
require('tanggle')
require('phangorn')### Main Pipeline#### Graphics%%R
framecentroid <- read.csv(paste0(outfolder, 'centroids_', name, '.csv'), header = TRUE, sep = ',')
framemeta <- read.csv(paste0(outfolder, 'meta_', name, '.csv'), header = TRUE, sep = ',')
framecluster <- read.csv(paste0(outfolder, 'cluster_', name, '.csv'), header = TRUE, sep = ',')
framevectors <- read.csv(paste0(outfolder, 'vectors_', name, '.csv'), header = TRUE, sep = ',')%%R
cluster <- inner_join(framecluster, framemeta, by = 'accession')
vectors <- inner_join(framevectors, cluster, by = 'accession')
centroid <- inner_join(framecentroid, framemeta, by = 'accession')%%R
pca <- read.csv(paste0(outfolder, 'pca_', name, '.csv'), header = TRUE, sep = ',')
info <- read.csv(paste0(outfolder, 'info_', name, '.csv'), header = TRUE, sep = ',')%%R
#tree <- read.tree(paste0('centroids_', name, '.fasta.tree'))
tree <- read.tree(paste0(outfolder, 'RAxML_bestTree.', name, '.tree'))%%R
listcolor <- c(
    "dodgerblue2", "#E31A1C",
    "green4",
    "#6A3D9A",
    "#FF7F00",
    "aquamarine4", "gold1",
    "skyblue2", "#FB9A99",
    "palegreen2",
    "#CAB2D6",
    "#FDBF6F",
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown", 'blue'
)%%R
if(color != ''){
    #frame <- centroid[ , c('accession', color)] %>% remove_rownames %>% column_to_rownames(var="accession") 
    subtype <- c(str_sort(levels(as.factor(centroid[[color]])), numeric = TRUE), 'mixed')

    centroid[[color]][is.na(centroid[[color]])] <- 'mixed'
    centroid[[color]] <- factor(centroid[[color]], levels=subtype)
    frame <- centroid[ , c('accession', color)] %>% remove_rownames %>% column_to_rownames(var='accession') 

    q <- length(subtype)-1
}%%R -w 1200 -h 1200
#options(repr.plot.width=20, repr.plot.height=20)

p0 <- ggtree(tree, layout='circular') %<+% centroid

p1 <- open_tree(p0, 0) +
    geom_treescale(width = 0.25) + 
    geom_tippoint(aes(size=size), alpha=.6) +
    geom_tiplab2(aes(label=paste0(cluster, ':', label)), align=TRUE) +#, parse=T) +
    scale_size_continuous(
        name = 'size',
        breaks = c(2, 10, 100, 1000, 10000), 
    )   

if(color != ''){
    p2 <- p1 +
         new_scale_fill()

    p3 <- gheatmap(
            p2,
            frame,
            offset=0.5, 
            width=0.05, 
            #font.size=3, 
            #colnames_angle=0, 
            colnames = TRUE, 
            #hjust=0
        ) +
        {if(color != '')scale_fill_manual(
            name = color,
            values = c(listcolor[0:q], 'black'),
            breaks = subtype, 
            labels = subtype,
            #na.value = "black",
            guide = guide_legend(override.aes = list(shape = 15, size = 5))
        )}
}else{
    p3 <- p1
}

print(p3)
ggsave(paste0(outfolder, 'centroid_', name, '.pdf'), width = 40, height = 40, units = 'cm', limitsize = FALSE)%%R
msa <- read.dna(paste0(outfolder, 'centroids_', name, '.msa.fasta'), format = 'fasta')
dm <- dist.ml(msa)
nnet <- neighborNet(dm)
#nnet <- read.nexus.networx(paste0('centroids_', name, '.nex'))%%R
if(color != ''){
    k <- sapply(nnet$tip.label, (\(x) centroid[which(centroid$accession == x),][[color]]), simplify = TRUE, USE.NAMES = TRUE)
}

c <- sapply(nnet$tip.label, (\(x) centroid[which(centroid$accession == x),]$cluster), simplify = TRUE, USE.NAMES = TRUE)
s <- sapply(nnet$tip.label, (\(x) centroid[which(centroid$accession == x),]$size), simplify = TRUE, USE.NAMES = TRUE)%%R -w 1200 -h 1200
#options(repr.plot.width=20, repr.plot.height=20)

p4 <- ggsplitnet(nnet) 

p5 <- p4 + 
    {if(color != '')geom_tiplab2(aes(label = paste0(c[label], ':', label), color = k[label]), hjust = -.5)} +
    {if(color == '')geom_tiplab2(aes(label = paste0(c[label], ':', label)), hjust = -.5)} +
    {if(color != '')geom_tippoint(aes(color = k[label], size = s[label]))} +
    {if(color == '')geom_tippoint(aes(size = s[label]))} +
    theme(legend.position = 'right') +
    {if(color != '')scale_color_manual(
        name = color,
        values = c(listcolor[0:q], 'black'),
        breaks = subtype, 
        labels = subtype,
        #na.value = "black",
        guide = guide_legend(override.aes = list(shape = 15, size = 5))
    )} +
    scale_size_continuous(
        name = 'size',
        breaks = c(2, 10, 100, 1000, 10000), 
    ) +
    ggexpand(.1) + ggexpand(.1, direction=-1)

print(p5)
ggsave(paste0(outfolder, 'nexus_', name, '.pdf'), width = 40, height = 40, units = 'cm', limitsize = FALSE)%%R
if(color != ''){
    subtype <- c(str_sort(levels(as.factor(vectors[[color]])), numeric = TRUE), 'mixed')
    vectors[[color]][is.na(vectors[[color]])] <- 'mixed'
    #vectors[[color]] <- ifelse(vectors$accession %in% centroid$accession, 'centroid', vectors[[color]])
    vectors[[color]] <- factor(vectors[[color]], levels=subtype)
    
    q <- length(subtype)-1
}

vectors$centroid <- with(vectors, ifelse(accession %in% centroid$accession, 'C', 'nC'))
vectorscentroid <- factor(vectors$centroid)%%R
#as.matrix(vectors[vectors$centroid == 'C', -1])
vec1 <- subset(vectors, centroid == 'nC')
vec2 <- subset(vectors, centroid == 'C')%%R
rownames(vec1) <- NULL
rownames(vec2) <- NULL%%R
fig <- if(color != ''){
    plot_ly(
        vec1,
        x = ~x,
        y = ~y,
        z = ~z,
        #symbol = ~centroid,
        #symbols = c('circle','x'),
        color = vec1[[color]], 
        colors = c(listcolor[0:q], 'black'),
        type = 'scatter3d', 
        mode = 'markers',
        #size = 1,
        marker = list(symbol = 'circle', size = 1),
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster),
        width = 1200,
        height = 1200
    )
} else {
    plot_ly(
        vec1,
        x = ~x,
        y = ~y,
        z = ~z,
        type = 'scatter3d', 
        mode = 'markers',
        #size = 1,
        marker = list(symbol = 'circle', size = 1),
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster),
        width = 1200,
        height = 1200,
        name = 'mixed',
        color = I('black')
    )
}

fig <- if(color != ''){
    fig %>% add_trace(
        data = vec2,
        x = ~x,
        y = ~y,
        z = ~z,
        #color = I('black'),
        color = vec2[[color]], 
        colors = c(listcolor[0:q], 'black'),
        type = 'scatter3d', 
        mode = 'markers',
        marker = list(symbol = 'x', size = 1),
        #name = 'centroid',
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster)
    )
} else {
    fig %>% add_trace(
        data = vec2,
        x = ~x,
        y = ~y,
        z = ~z,
        color = I('black'),
        marker = list(symbol = 'x', size = 1),
        name = 'mixed',
        text = ~paste('Host:', host, '<br>Subtype:', subtype, '<br>Year:', year, '<br>Accession:', accession, '<br>Cluster:', cluster)
    )
}

fig <- fig %>% layout(
    scene = list(
        xaxis = list(
            title = 'PCA1',
            #gridcolor = 'rgb(255, 255, 255)',
            #zerolinewidth = 1,
            #ticklen = 5,
            #gridwidth = 2,
            range = c(-1.0, 1.0)
        ),
        yaxis = list(title = 'PCA2',
            range = c(-1.0, 1.0)
        ),
        zaxis = list(title = 'PCA3',
            range = c(-1.0, 1.0)
        ),
        #paper_bgcolor = 'rgb(243, 243, 243)',
        #plot_bgcolor = 'rgb(243, 243, 243)',
        aspectmode = 'cube'
    ),
    autosize = TRUE,
    margin = c(l=0, r=0, b=0, t=0), 
    #title = 'Life Expectancy v. Per Capita GDP, 2007',
    legend = list(
        orientation = 'v',
        yanchor = 'middle',
        y = 0.5,
        itemsizing='constant'
    )
)

embed_notebook(fig)

htmltools::save_html(fig, paste0(outfolder, 'sphere_', name, '.html'))%%R
reduced <- info$components%%R -w 1200 -h 400

p6 <- ggplot(data=pca, aes(x=components, y=variance)) +
    geom_line() +
    geom_point() +
    theme_light() + 
    geom_vline(xintercept=c(3, reduced), linetype='dashed', color=c('blue', 'red'))
    
p7 <- p6 + annotate(x=c(3, info$components),y=+Inf,label=c(pca[pca$components==3,]$variance, pca[pca$components==reduced,]$variance),vjust=2,geom="label")

ggplotly(p7)
ggsave(paste0(outfolder, 'pca_', name, '.png'), width = 40, height = 10, units = "cm", limitsize = FALSE)