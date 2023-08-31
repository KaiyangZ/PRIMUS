#!/usr/bin/env Rscript

# simulations were performed using splatter-1.16 

source("splatPopSim.R")

n = 20

#' simulation scenario 1 
dir.create("simulation_1")

# generate seed 
set.seed(123)
sim1.seeds = sample(seq.int(100, 9999), n)

for (i in seq.int(n) ){
    s = sim1.seeds[i]

    sim = splatPopSim(vcf.n.samples = 6,
                      vcf.n.snps = 20000,
                      bulk.n.genes = 5000,
                      bulk.n.samples = 100,
                      group.prob = rep(1/5, 5),
                      nGenes = 5000,
                      batchCells = 500,
                      seed = s
                      )

    saveRDS(sim, paste0( "simulation_1/simData_", i, ".RDS") )
}
 
#' simulation scenario 2 
dir.create("simulation_2")

# scenario design

#       sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
#Group1        1        1        0        1        0        0
#Group2        1        0        1        0        1        1
#Group3        0        1        0        0        0        1
#Group4        1        0        1        1        1        0
#Group5        0        0        1        0        0        0

# generate seed 
set.seed(345)
sim2.seeds = sample(seq.int(100, 9999), n)

group2samples = lapply( list(c(1, 2, 4), c(1, 3, 5, 6), c(2, 6), c(1, 3, 4, 5), c( 3) ), function(x) paste0("sample_", x) )

names(group2samples) = paste0("Group", seq.int(5))

for (i in seq.int(n) ){
    s = sim2.seeds[i]

    sim = splatPopSim(vcf.n.samples = 6,
                      vcf.n.snps = 20000,
                      bulk.n.genes = 5000,
                      bulk.n.samples = 100,
                      group.prob = rep(1/5, 5),
                      nGenes = 5000,
                      batchCells = 500,
                      group2samples = group2samples,
                      seed = s
                      )

    saveRDS(sim, paste0( "simulation_2/simData_", i, ".RDS") )
}

#' simulation scenario 3
dir.create("simulation_3")

# generate seed 
set.seed(678)
sim3.seeds = sample(seq.int(100, 9999), n)

# define size of each sample 
# 2000, 1000, 500, 200, 100, 20

sampling.n = c("sample_1" = 2000, "sample_2" = 1000,
               "sample_3" = 500, "sample_4" = 200,
               "sample_5" = 100, "sample_6" = 20
               )

for (i in seq.int(n) ){
    print(i)
    s = sim3.seeds[i]

    sim = splatPopSim(vcf.n.samples = 6,
                      vcf.n.snps = 20000,
                      bulk.n.genes = 5000,
                      bulk.n.samples = 100,
                      group.prob = rep(1/5, 5),
                      nGenes = 5000,
                      batchCells = 5000,
                      group2samples = group2samples,
                      sampling.n = sampling.n,
                      seed = s
                      )

    saveRDS(sim, paste0( "simulation_3/simData_", i, ".RDS") )
}
