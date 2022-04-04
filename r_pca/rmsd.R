library(bio3d)
library(ape)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

devtools::install_github('coolbutuseless/minilexer')
library(parallel)
library(dendextend)
library(pvclust)
#devtools::install_github('coolbutuseless/minilexer')
library(minilexer)
library(plot3D)
library(rgl)
library(car)
library(future.apply)
plan(multiprocess) ## Run in parallel on local computer

library(tidyverse)
library(doParallel)
tmalign <- function(fixed, mobile){
  capture.output(
    cat(
      system(
        paste(
          "TMalign", 
          fixed, 
          mobile
        ), 
        intern = TRUE
      ),
      collapse = "\n"
    )
  ) 
}

tmalign_parser <- function(data){
  require(magrittr)
  tokens <- c(
    fixed      = "Name\\s+of\\s+Chain_1:\\s+\\S+",
    mobile     = "Name\\s+of\\s+Chain_2:\\s+\\S+",
    fixed_len  = "Length\\s+of\\s+Chain_1:\\s+\\d+",
    mobile_len = "Length\\s+of\\s+Chain_2:\\s+\\d+",
    ali_len    = "Aligned\\s+length=\\s+\\d+",
    RMSD       = "RMSD=\\s+\\d*\\.?\\d+",
    Seq_ID     = "n_aligned=\\s+\\d*\\.?\\d+", 
    TM_score_1 = "\\d*\\.?\\d+\\s+\\(if\\s+normalized\\s+by\\s+length\\s+of\\s+Chain_1", 
    TM_score_2 = "\\d*\\.?\\d+\\s+\\(if\\s+normalized\\s+by\\s+length\\s+of\\s+Chain_2"
  ) 
  
  data %>%
    stringr::str_extract_all(tokens) %>% 
    tidyr::as_tibble(.name_repair = ~names(tokens)) %>%
    dplyr::mutate(
      fixed      = stringr::str_extract(fixed, "\\S+$"),
      mobile     = stringr::str_extract(mobile, "\\S+$"),
      fixed_len  = as.integer(stringr::str_extract(fixed_len, "\\d+$")),
      mobile_len = as.integer(stringr::str_extract(mobile_len, "\\d+$")),
      ali_len    = as.integer(stringr::str_extract(ali_len, "\\d+$")),
      RMSD       = as.double(stringr::str_extract(RMSD, "\\d*\\.?\\d+$")),
      Seq_ID     = as.double(stringr::str_extract(Seq_ID, "\\d*\\.?\\d+$")),
      TM_score_1 = as.double(stringr::str_extract(TM_score_1, "^\\d*\\.?\\d+")),
      TM_score_2 = as.double(stringr::str_extract(TM_score_2, "^\\d*\\.?\\d+"))
    ) -> data
  data
}

pdb2tmalign <- function(f,m) tmalign_parser(tmalign(f,m))


pdb2tmalign_df <- function(data, col_index_mobile, col_index_fixed, cpu = 1){
  require(doParallel)
  cl <- makeCluster(cpu)
  clusterExport(cl, c("data","foo", "tmalign", "tmalign_parser"))
  res <- foreach(i = 1:NROW(data)) %dopar% {
    pdb2tmalign(data[i,col_index_mobile],data[i,col_index_fixed]) 
  }
  stopCluster(cl)
  Reduce("rbind", res)
}

#pdb2tmalign_df(dd[c(1:10),], 1, 2, 12)

tmalign_list <- function(fixedList, mobileList, cpu = 1){
  data <- tidyr::crossing( fixed=fixedList, mobile=mobileList)
  pdb2tmalign_df(data=data, 1, 2, cpu=cpu)
}

tmalign_allvsall <- function(l, cpu = 1){
  tmalign_list(fixedList = l, mobileList = l, cpu = cpu)
}

#setwd("/home/acpguedes/projects/sig_trans/work/SBP_3/")
setwd("/home/acpguedes/projects/sig_trans/work/SBP_3/pdb_ali/")
pdbs <- read_tsv("/home/acpguedes/projects/sig_trans/work/SBP_3/pbps_pdb.rawtable")
ids <- pdbs %>% group_by(id_0.3cov_0.8evalue_0.001) %>% top_n(1) %>% pull(ID)
chey <- c("4WXM_A", "5CHY_A", "4QIC_A", "3C3M_A", "5M7N", "4TMY_A", "1QMP",   "1ORD_A", "1C4K_A", "5FKX_A", "2AYX_A", "2KX7_A", "1SR2_A", "4LDA_A", "4JC0_A")
ids <- c(ids, "5UB6_A", chey)
mark <- c(rep(1, 223), rep(2, 15))
dat <- tibble(id = ids, fold = mark)

raw.files <- get.pdb(dat$id)
files <-pdbsplit(raw.files, dat$id)

tidyr::crossing( fixed=files, mobile=files) %>% write_tsv("sbps_chey.tsv")
read_file("sbps_chey.tmalign") -> tmalign 

# TM-align

t1 <- Sys.time()
  ts <- tmalign_allvsall(files, 12)
t2 <- Sys.time()
t2-t1

          foo(.[1,1], .[1,2])

    tidyr::crossing( fixed=files, mobile=files) %>%
  head() %>%
  tmalign(fixed = fixed, mobile = mobile)

get_rmsd <- function(data){
  tmalign("split_chain/2FQX_A.pdb", "split_chain/2FQX_A.pdb") %>%
    tmalign_parser() 
}

%>%
  mutate(
    fixed = str_extract(fixed, "[A-Za-z0-9]{4}_\\w(?=\\.)"),
    mobile = str_extract(mobile, "[A-Za-z0-9]{4}_\\w(?=\\.)")
  ) -> t


t %>%
  select(fixed,mobile,RMSD) %>%
  spread(key = mobile, value = RMSD) -> t.rmsd

rownames(t.rmsd) <- t.rmsd$fixed
ids <- rownames(t.rmsd) 
t.rmsd %>%
  select(-fixed) %>%
  as.matrix() -> t.rmsd
  
p <- pvclust(t.rmsd, method.dist="euclidian", 
        method.hclust="ward.D", nboot=100)
plot(p)
pvrect(p)

p.pca <- pca(t.rmsd)  
scatter3d(x=p.pca$U[,1],y=p.pca$U[,2],z=p.pca$U[,3], col.var = dat[ order(dat$id, rownames(t.rmsd)) ,]$fold )

  plot(p.pca)
  
  
  
mysum <-function(a, b){
    a + b
}
d <- tidyr::tibble(a = seq(1,10,1), b = seq(11,20,1))
d %>%
  head(n=1L) %>%
  group_by(a,b) %>%
  mysum(.$a,.$b)
