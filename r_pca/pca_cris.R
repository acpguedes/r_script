  # Check if a LIB exist, if not, install it
  check_lib <- function(mypkg, install_pkg=1){
    if(is.element(mypkg, installed.packages()[,1])){
      suppressMessages(library(mypkg, character.only=TRUE))
    } else {
      if(install_pkg == 1){
        install.packages(mypkg, character.only=TRUE, repos="http://vps.fmvz.usp.br/CRAN/")
        suppressMessages(library(mypkg, character.only=TRUE))
      } else {
        stop(paste("Missing package:", mypkg, sep=" "), call.=FALSE)
      }
    }
  }
  
  
  # Load libs or install if it not exist
  check_lib("seqinr") # to deal with fasta format and fetch alignment info
  check_lib("tidyverse") # easy to work with dataframes also provides ggplot 
  check_lib("plotly")   # allow 3D plot 
  check_lib("heatmaply")
  check_lib("devtools") # toi load library in development as the one load below
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  
  # Load data
  # alignment
  AlnFile <- "/home/acpguedes/Downloads/prunMb.fa"
  GroupsFile <- "/home/acpguedes/Downloads/grupos_pca.txt"
  aln <- read.alignment(AlnFile, format = "fasta")
  # table with genomes and groups
  aln_groups <- read_tsv(GroupsFile, 
                         col_names = c("genome", "group"),
                         col_types = c(
                           genome=col_character(),
                           group=col_integer()
                           )) 
    
  # edit table to add color by group
  aln_groups %>%
    mutate(
      color = case_when(
        group == 1 ~ "red",
        group == 2 ~ "purple",
        group == 3 ~ "gray",
        group == 4 ~ "green",
        group == 5 ~ "pink",
        group == 6 ~ "black",
        group == 7 ~ "blue",
        TRUE ~ NA_character_
      ),
      lineage = case_when(
        group == 1 ~ "group 1",
        group == 2 ~ "group 2",
        group == 3 ~ "group 3",
        group == 4 ~ "group 4",
        group == 5 ~ "group 5",
        group == 6 ~ "group 6",
        group == 7 ~ "group 7",
        TRUE ~ NA_character_
      )
      ) -> aln_groups
  
  # extract a ordered vecto with colors
  aln_groups %>%
    pull(color) -> aln_groups_color
  
  # calculate the identity between each pair of genomes based on the alignment
  ali_dist <- dist.alignment(aln, matrix = "identity")
  #ali_dist2 <- dist.alignment(aln, matrix = "similarity")
  # calculate a distance matrix based on z-score (normalized) calculated up to the identity
  #ali_dist_matrix <- as.matrix(qnorm(ali_dist))
  
  # PCA : Spectral decomposition which examines the covariances 
                                   #/ correlations between variables
  # Used identity matrix
  b <- princomp(as.matrix(ali_dist))
  plot(b)
  # get a ordered vector based on pca order
  aln_groups_color <- aln_groups %>%
    arrange(match(genome, rownames(b$scores))) %>%
    pull(color)
  
  # Base plot (without ggplot, you may try it)
  plot(b$scores)
  points(b$scores[,1:2], col=aln_groups$color, pch=16, cex=0.9)
  
  # Convert pcs to dataframe and add colorgroups to customized plot
  b$scores %>% 
    as_tibble(rownames = "ID") %>% 
    select(ID, everything()) %>%  
    left_join( aln_groups, by = c("ID" = "genome") ) %>%
    select(Comp.1, Comp.2, Comp.3, ID, group, color, lineage) -> pca.cris.tmp
  
  # ggplot: easy to edit parameters
  pca.cris.tmp %>% 
    mutate(size=if_else(color == "green", 8, 1)) %>% 
    ggplot(aes(x=Comp.1, y=Comp.2, color=lineage)) +
    geom_point(color = pca.cris.tmp$color) +
#    scale_color_manual(c(
#      group = "red",
      group = "purle",
  d    group = "gray",
      group = "green",
      group = "pink",
      group = "black",
      group = "blue"
      )) +
    xlab("Comp.1 (43.9%)") +
    ylab("Comp.2 (27.4%)") +
    theme_classic()
  
          # 3D scaterplot to use the 3 first components
  plot_ly(pca.cris.tmp,
          x = ~Comp.1, y = ~Comp.2, z = ~Comp.3,
  #        group = ~group,
  #        marker = list(color = ~color),
          color = ~lineage, 
          colors = c("red", "purple", "grey", "grgreen", ey", "pink", "black", "blue") ,
          type = "scatter3d", mode = "markers", 
          #alpha = 0.5, 
          size = 1 ) 
  
  # SUGGESTION to cluster vizualization (if needed)
  
  # An another alternative to see the clusters is the supervised method which is k-means
  # k-means might be performed in any kind hierarquical clustering or tree (even ML and MP)
  
  # In this example below I clusterized it using wardD2 algorithim to build a 
  # distance based tree. Then, I used the identity values to build a heatmap.
  # Additionally I added the left bars to point where each leaf belongs 
  # throughout the manual and automatic (k-means with 4 groups) annotations.
  
  # The tree might be changed by the tree in the paper.
  
  my.plot <- function(matdist, alg){
    
    # perform cluster
    hc <- hclust(as.dist(matdist), method = alg)
    # get ordered vector of colors based in hclust
    aln_groups_color <- aln_groups %>%
        arrange(match(genome, rownames(hc$labels))) %>%
      pull(color)
    # apply k-means with 4 groups to the cladogram
    ct1 <- cutree(hc, k=4)
    # combine with previous groups  
    vec <- cbind( "4-means" = as.character(ct1),
                  "groups"  = aln_groups_color) 
    # transposed table
    vac <- t(vec)
    # plot heatmap
    heatmap.3(matdist, distfun = as.dist, 
              hclustfun = function(x) hclust(x, method = alg),
              # ColSideColors= vec, 
              RowSideColors= vac,
              #ColSideColorsSize=as.numeric(ncol(vec)) , 
              RowSideColorsSize=as.numeric(nrow(vac)), 
              density.info = "density")
    
  }
  
  my.plot(ali_dist, "average")
  
  
  alg <- "average"
  matdist <- ali_dist
  hc <- hclust(as.dist(matdist), method = alg)
  # get ordered vector of colors based in hclust
  aln_groups_color <- function(){
    aln_groups %>%
      arrange(match(genome, rownames(hc$labels))) %>%
      pull(color)
  }
  # apply k-means with 4 groups to the cladogram
  ct1 <- cutree(hc, k=4)
  # combine with previous groups  
  vec <- cbind( "4-means" = as.character(ct1),
                "groups"  = aln_groups_color) 
  # transposed table
  vac <- t(vec)
  # plot heatmap
  heatmap.3(matdist, distfun = as.dist, 
            hclustfun = function(x) hclust(x, method = alg),
            # ColSideColors= vec, 
            RowSideColors= vac,
            #ColSideColorsSize=as.numeric(ncol(vec)) , 
            RowSideColorsSize=as.numeric(nrow(vac)), 
            density.info = "density")
  
    
  
  
  matdist %>%
    as.matrix() %>%
    as_tibble(rownames = "ID") %>%
    left_join( aln_groups, by = c("ID" = "genome") ) %>%
    select(ID, group, color, lineage, everything()) -> x
  
  heatmaply(
    as.matrix(x[, -c(1:4)]), 
    seriate = 'mean', 
    row_side_colors = x[,3], 
    k_col = 5, 
    k_row = 4,
    dist_method = "euclidian",
    hclust_method = "complete"
  )
  
  # https://stackoverflow.com/questions/15153202/how-to-create-a-heatmap-with-a-fixed-external-hierarchical-cluster
  # https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html 
  # https://davetang.org/muse/2013/05/09/using-the-r-seqinr-package/
  
  
  
