calculate_new_worker_infections <- function(c.inf, w.inf, params){
  #c.inf <- split_arrays(c.inf)
  w <- worker_calculate_infectiousness(c.inf[1], c.inf[2], c.inf[3], 
                                       w.inf[1], w.inf[2], w.inf[3], 
                                       params)
  return(1-exp(-w))
}



#calculate for all workers
worker_calculate_infectiousness <- function(community.asx, community.sx, community.n,
                                            workers.asx, workers.sx, workers.n,
                                            params){
  with(params,{
    #assume 75% of interactions come from patients and 25% comes from other workers
    asx <- sum(c(0.75*lambda.cw*community.asx/community.n, 
                    0.25*lambda.ww*workers.asx/workers.n), na.rm = TRUE)
    sx <- sum(c(0.75*lambda.cw*community.sx/community.n,
                0.25*lambda.ww*(workers.sx)/workers.n), na.rm = TRUE)
    p <- alpha.a*asx + sx
    return(p)
  })
}
