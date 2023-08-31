#!/usr/bin/env Rscript

splatPopSim <- function(vcf.n.samples = 5, vcf.n.snps = 2000, 
			bulk.n.genes = 2000, bulk.n.samples = 50,
			group.prob, nGenes = 2000, batchCells = 500, 
			group2samples = NULL, sampling.n = NULL, 
			seed = 123, ...){

    require(splatter)
    require(scater)
    require(VariantAnnotation)

    vcf <- mockVCF(n.samples = vcf.n.samples, n.snps = vcf.n.snps, seed = seed)
    bulk.eqtl <- mockBulkeQTL(n.genes = bulk.n.genes, seed = seed)
    bulk.means <- mockBulkMatrix(n.genes = bulk.n.genes, n.samples = bulk.n.samples, seed = seed)
    
    params.est <- splatPopEstimate(means = bulk.means,
                                   eqtl = bulk.eqtl,
                                   counts = mockSCE()
    				   )
    
    params.est <- setParams(params.est, nGenes = nGenes, group.prob = group.prob, batchCells = batchCells, similarity.scale = 6)
    
    sim.sc.est <- splatPopSimulate(vcf = vcf, params = params.est, seed = seed)

    if (!is.null(group2samples ) ){
	keep = sapply(names(group2samples), function(x) sim.sc.est$Group == x & sim.sc.est$Sample %in% group2samples[[x]] ) 	    	 
        
    	sim.sc.est = sim.sc.est[, rowSums(keep) == 1]    
    }

    sim.annot = colData(sim.sc.est) 

    if(!is.null(sampling.n) ){
	    require(sampling)
	    set.seed(seed)
	    sampling.n = sampling.n[unique(sim.annot$Sample) ]

 	    sim.annot.ds = sampling::strata(sim.annot, stratanames = "Sample", size = sampling.n, method = "srswor")	
        sim.sc.est = sim.sc.est[, sim.annot.ds$ID_unit]    
    }

    return(sim.sc.est)
}

