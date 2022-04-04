library(bio3d)
library(tidyverse)
library(magrittr)
library(cluster)
library(factoextra)
library(patchwork)
library(ggtree)
library(ggtreeExtra)
library(ggplotify)
library(grid)
library(ggheatmap)
library(ComplexHeatmap)
# library(da)


#library(phangorn)
setwd("/home/acpguedes/projects/sig_trans/work/SBP_3/")
# pdbs <- read_tsv("/home/acpguedes/projects/sig_trans/work/SBP_3/pbps_pdb.rawtable")
setwd("/home/acpguedes/projects/sig_trans/work/SBP_3/pdb_ali/")
source("/home/acpguedes/projects/scripts/heatmap.3.R")

pdbs <- read_tsv("/home/acpguedes/projects/sig_trans/work/SBP_5/20190730/test2/sbps.rawtable.full")

pdbs %<>% 
    mutate(
        is.PDB = str_detect(pid, "[\\w\\d]{4}_[\\w\\d]")
    ) %>% 
    filter(is.PDB) 

pdbs %>% 
    group_by(id_0cov_0.8evalue_0.001) %>% 
    slice_min(1) %>% 
    pull(pid) -> ids

# ids <- pdbs %>% pull(pid) 
#ids <- c(ids, "5UB6_A")

raw.files <- get.pdb(ids)
files <-pdbsplit(raw.files, ids)

#sample <- sample(files, 5)
#dput(sample)

to_ali <- tidyr::crossing( fixed=files, mobile=files)
write_tsv(to_ali, "alignment_pair.tsv")

as_tibble(data.table::fread("/home/acpguedes/projects/sig_trans/work/SBP_5/20190730/test2/hhsuite_teste_SBP/sbps_node_annotation_leiden.tsv")) -> annot

# #,
#          col_names = c("node.name","node","rep","pid","acc","full",
#                        "e0.1","e0.001","e-5","e-6","e-10","e-15","e-20",
#                        "id_0cov_0.8evalue_0.001","id_0cov_0.8evalue_1e-05",
#                        "id_0cov_0.8evalue_1e-06","id_0cov_0.8evalue_1e-10",
#                        "id_0cov_0.8evalue_1e-15","id_0cov_0.8evalue_1e-20",
#                        "id_0cov_0.8evalue_1e-30","id_0cov_0.8evalue_1e-50",
#                        "id_0cov_0.8evalue_1e-100","id_0cov_0.8evalue_1e-150",
#                        "id_0cov_0.8evalue_1e-200","nr","node.size","architecture.x",
#                        "regions.x","evalue.x","architecture.y","regions.y","evalue.y",
#                        "Cluster","Subcluster","Protein","Ligand","k2","k3","k7","k19",
#                        "PDB","block_id","ST","cache_neighbor","STrel")) -> annot

annot %>%
    filter(acc %in% ids) %>%
    group_by(Cluster) %>%
    mutate(Clusternum = group_indices()) %>%
    # ungroup() %>%
    # group_by(Subcluster) %>%
    # mutate(Subclusternum = group_indices()) %>%
    ungroup() -> annot

annot %>% pull(acc) %>% unique() -> annot.acc

#' SYS: time tail -n +2 alignment_pair.tsv | parallel -j 40 --colsep '\t' TMalign {1} {2} > alignment_pair.TMalign.out 
#' SYS: (echo -e "fixed\tmobile\trmsd" ; perl -ne 'print $1 . "\t" if /^Name of Chain_1:\s+split_chain\/(\S+)\.pdb/; print $1 . "\t" if /^Name of Chain_2:\s+split_chain\/(\S+).pdb/; print $1 . "\n" if /RMSD=\s+(\S+),/' alignment_pair.TMalign.out) > alignment_pair.TMalign.tsv 

read_tsv("alignment_pair.TMalign.tsv") -> pbps_rmsd
#read_tsv("sbps_pdb.rmsd.tsv") -> pbps_rmsd

pbps_rmsd %>%
    filter(mobile %in% annot.acc & fixed %in% annot.acc) %>%
    inner_join(pbps_rmsd, by = c("fixed"="mobile", "mobile"="fixed") ) %>%
    mutate(
           rmsd = (rmsd.x + rmsd.y)/2
    #    fixed = str_extract(str_extract(fixed, '[A-Za-z0-9]{4}_\\w\\.'),'[^\\.]+'),
    #    mobile = str_extract(str_extract(mobile, '[A-Za-z0-9]{4}_\\w\\.'),'[^\\.]+')
           ) %>% 
    #filter(mobile %in% annot.acc & fixed %in% annot.acc) %>%View()
    select(fixed, mobile, rmsd) -> pbps_rmsd

pbps_rmsd.acc <- pbps_rmsd %>%
    filter(fixed %in% ids) %>% 
    filter(mobile %in% ids) %>% 
    pull(fixed)

pbps_rmsd %>%
    filter(fixed %in% ids) %>% 
    filter(mobile %in% ids) %>% 
    spread(key = mobile, value = rmsd) %>%
    remove_rownames() %>%
    column_to_rownames(var = "fixed") %>%
    as.matrix(.[-1]) -> pbps_rmsd_mat

pbps_rmsd_mat_ids <- rownames(pbps_rmsd_mat)
pbps_rmsd_mat_dist <- as.dist(pbps_rmsd_mat)

pbps_rmsd_mat_dist <- get_dist(pbps_rmsd_mat)
pbps_rmsd_mat_dist_pearson <- get_dist(pbps_rmsd_mat, method = "pearson")


auto_clust <- function(
    data, 
    methods  = c("average", "single", "complete", "ward"), 
    metric   = "euclidean",
    fundist  = factoextra::get_dist,
    funclust = function(x, y) cluster::agnes(x = x, method = y), 
    ...
){
    # avoid NA values, use method gower to distances
    if(any(class(data) == "dist")){
        dh_dist <- data
    }else{
        dh_dist <- fundist(data, metric) #daisy(data , metric="euclidean")
    }
    # methods to try
    names(methods) <- methods
    # get all clusters methods
    methods %>%
        map( function(x) funclust(x = dh_dist, y = x) ) -> myclusters
    # select best method
    myclusters %>%
        map(function(x) list(x = x$ac)) %>%
        enframe() %>%
        unnest(value) %>%
        unnest(value) %>%
        filter(value == max(value)) %>%
        pull(name) -> mymethod
    myclusters[[mymethod]]
}


# hc.nma <- hclust(pbps_rmsd_mat_dist^2, method = "ward.D")
hc.nma <- auto_clust(pbps_rmsd_mat, metric  = "euclidean")


#hclustplot(hc.nma, k = 3, labels = ids)

grps.nma1 <- cutree(hc.nma, k=2)
grps.nma2 <- cutree(hc.nma, k=3)
grps.nma3 <- cutree(hc.nma, k=6)
# grps.nma4 <- cutree(hc.nma, k=7)
# grps.nma4 <- cutree(hc.nma, k=8)
vec <- cbind(
    ids           = pbps_rmsd_mat_ids,
    `k2`          = as.character(grps.nma1),
    # Cluster       = as.character(annot[match(pbps_rmsd_mat_ids, annot$acc),]$Clusternum),
    `k3`          = as.character(grps.nma2),
    `k6`          = as.character(grps.nma3),
    # `k7`          = as.character(grps.nma4),
    Clustername   = as.character(annot[match(pbps_rmsd_mat_ids, annot$acc),]$Cluster),
    full          = as.character(annot[match(pbps_rmsd_mat_ids, annot$acc),]$full)
    # `m1e-1`       = as.character(annot[match(pbps_rmsd_mat_ids, annot$acc),]$e0.1 ),
    # `m1e-5`       = as.character(annot[match(pbps_rmsd_mat_ids, annot$acc),]$`e-5` )
    )

vec %>% 
    as_tibble() %>%
    mutate(
        k2 = case_when(
            k2 == "1" ~ "green", #PBPII
            k2 == "2" ~ "red"    #PBPI/III
        ), 
        k3 = case_when(
            k3 == '1' ~ 'orange', #PBPII
            k3 == '2' ~ 'firebrick',   # PBPI
            k3 == '3' ~ 'midnightblue'   # PBPIII
        ),
        k6 = case_when(
            k6 == '1' ~ 'darkgreen',    # PBPII D/G
            k6 == '2' ~ 'firebrick',    # PBPI B
            k6 == '3' ~ 'darkorange',   # PBPII F/E 
            k6 == '4' ~ 'purple',       # PBPIII C
            k6 == '5' ~ 'midnightblue', # PBPIII A
            k6 == '6' ~ 'orange',   # PBPII E
            # k7 == '7' ~ 'orange'        # PBPII E/F
        ),
        Cluster = case_when(
            Clustername == 'A' ~ 'midnightblue',     # PBPIII A
            Clustername == 'B' ~ 'firebrick',      # PBPI   B
            Clustername == 'C' ~ 'purple',   # PBPII  C
            Clustername == 'D' ~ 'darkgreen',    # PBPII  D
            Clustername == 'E' ~ 'orange',   # PBPII  E
            Clustername == 'F' ~ 'darkorange',     # PBPII  F
            Clustername == 'G' ~ 'lightgreen', # PBPII  G
        ),
        full = case_when(
            full == '1' ~ 'darkgreen',      # PBPII E/F
            full == '2' ~ 'darkorange',     # PBPII D
            full == '3' ~ 'firebrick',       # PBPI/III A/B
            full == '4' ~ 'green',  # PBPII D/G
            full == '5' ~ 'purple',    # PBPII C
            full == '6' ~ 'midnightblue'    # PBPII C
        ) 
        # Cluster = case_when(
        #     Cluster == '1' ~ 'blue',      # PBPIII A
        #     Cluster == '2' ~ 'red',       # PBPI   B
        #     Cluster == '3' ~ 'purple',    # PBPII  C
        #     Cluster == '4' ~ 'green',     # PBPII  D
        #     Cluster == '5' ~ 'orange',    # PBPII  E
        #     Cluster == '6' ~ 'darkorange' # PBPII  F
        # ),
        # `m1e-1` = case_when(
        #     `m1e-1` == '1' ~ 'orange', # PBPII E/F
        #     `m1e-1` == '2' ~ 'green', # PBPII D
        #     `m1e-1` == '3' ~ 'red', # PBPI
        #     `m1e-1` == '4' ~ 'darkgreen', # PBPII D/G
        #     `m1e-1` == '5' ~ 'purple', # PBPII C
        #     `m1e-1` == '6' ~ 'blue', # PBPIII
        # ),
        # `m1e-5` = case_when(
        #     `m1e-5` == '1'  ~ 'darkorange', # PBPII F
        #     `m1e-5` == '2'  ~ 'green', # PBPII D
        #     `m1e-5` == '3'  ~ 'darkgreen', # PBPII D/G
        #     `m1e-5` == '4'  ~ 'red', # PBPI
        #     `m1e-5` == '5'  ~ 'darkred', # PBPI
        #     `m1e-5` == '6'  ~ 'purple', # PBPII C
        #     `m1e-5` == '7'  ~ 'blue', # PBPIII Peripla_BP_2
        #     `m1e-5` == '8'  ~ 'lightgreen', # PBPII D
        #     `m1e-5` == '9'  ~ 'orange', # PBPII E
        #     `m1e-5` == '10' ~ 'darkblue', # PBPIII ZnuA
        # ),
    ) %>% 
    select(ids, Cluster, k2, k3, k6, full ) -> vec_tmp

vec_tmp %>% 
    as.data.frame %>% 
    column_to_rownames(var = ids) -> vec_tmp

vec_tmp %>%
    select(-ids) %>% 
    # select(-Clustername) %>%
    as.matrix() -> vec

#    `k19`         = as.character(grps.nma4),
#    `k19`         = as.character(grps.nma4),
#    `e0.001`      = as.character(pull(annot, e0.001))
#    `Cluster`     = as.character(pull(annot, Cluster))
#`Subcluster`  = as.character(pull(annot, Subcluster))
tvec <- t(vec)

as.ggplot(
    heatmap.3(pbps_rmsd_mat, distfun = daisy, labRow = ids, labCol = ids,
              hclustfun = auto_clust,
              # ColSideColors= vec,
              RowSideColors= tvec,
              # ColSideColorsSize=as.integer(ncol(vec)),
              RowSideColorsSize=as.integer(nrow(tvec)),
              density.info = "density" , breaks = seq(0,7, by = 0.5), sepcolor = "rainbow"
    )
) 


#as.data.frame(pbps_rmsd_mat) %>% rownames
pp <- as.dendrogram(hc.nma)
ggtree(pp) +
    geom_tiplab() -> my_tree

vec_tmp %>% 
    gather(key = key, value = value, -ids) -> vec_tmp_td

my_tree +
    geom_fruit(
        data = vec_tmp_td,
        mapping = aes(key, ids, col = value),
        geom = geom_tile
    )

as.data.frame(pbps_rmsd_mat) -> tt
rownames(tt) <- vec_tmp$ids
colnames(tt) <- vec_tmp$ids
gheatmap(p = my_tree, tt, low = "red", high = "yellow")

ggheatmap(pbps_rmsd_mat)
heatmap.3(pbps_rmsd_mat, distfun = daisy, labRow = ids, labCol = ids,
          hclustfun = auto_clust,
          # ColSideColors= vec,
          RowSideColors= tvec,
          # ColSideColorsSize=as.integer(ncol(vec)),
          RowSideColorsSize=as.integer(nrow(tvec)),
          density.info = "density" , breaks = seq(0,7, by = 0.5), sepcolor = "rainbow"
          ) 

tvec <- t(vec[, -1])
as.ggplot(
    function() heatmap.3(pbps_rmsd_mat, distfun = daisy, labRow = ids, labCol = ids,
                         hclustfun = auto_clust,
                         # dendrogram = as.dendrogram(hc.nma),
                         # ColSideColors= vec,
                         RowSideColors= tvec,
                         # ColSideColorsSize=as.integer(ncol(vec)),
                         RowSideColorsSize=as.integer(nrow(tvec)),
                         density.info = "density" , breaks = seq(0,7, by = 0.5), sepcolor = "rainbow"
    ) 
) + #Cluster
    labs(title = "a) SBPs PDB relationship") -> hm


q1 <- as.ggplot( function() plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,1])) + #Cluster
    labs(title = "b) Cluster Sheepers et al")
q2 <- as.ggplot( function() plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,3])) + #k3
    labs(title = "c) 3 kmers RMSD") 
q3 <- as.ggplot( function() plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,4])) + #k6
    labs(title = "d) 6 kmers RMSD")
q4 <- as.ggplot( function() plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,5])) + #full
    labs(title = "e) community")

 hm / ((q1 + q2) / (q4 + q3)) + patchwork::plot_annotation() 
    ylab("PC2 (22.27%)")


svg("/home/acpguedes/Pictures/Cluster_PCA.svg", width = 1127, height = 579)
plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,1]) #Cluster
dev.off()
svg("/home/acpguedes/Pictures/k3_PCA.svg", width = 1127, height = 579)
plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,3]) #k3
dev.off()
svg("/home/acpguedes/Pictures/full_PCA.svg", width = 1127, height = 579)
plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,5]) #full
dev.off()
svg("/home/acpguedes/Pictures/m1e1_PCA.svg", width = 1127, height = 579)
plot(pca(pbps_rmsd_mat), pc.axes = 1:2, col = vec[,6]) #m1e-1
dev.off()



vec_tmp %>%
    distinct() %>%
    arrange(Clustername,k7) %>%
    write_tsv("/home/acpguedes/Documents/colortab.tsv")
