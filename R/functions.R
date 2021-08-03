require('geomorph')
require('rgl')

getCoresNumber = function()
{
  require('doParallel')
  return(detectCores())
}

vcrossp <- function(a, b, normalize=F) {
  if (is.null(nrow(a)) || is.null(nrow(b))) {
    a = matrix(a,nrow=1)
    b = matrix(b,nrow=1)
  }
  result <- matrix( NA, nrow(a), 3 )
  result[,1] <- a[,2] * b[,3] - a[,3] * b[,2]
  result[,2] <- a[,3] * b[,1] - a[,1] * b[,3]
  result[,3] <- a[,1] * b[,2] - a[,2] * b[,1]
  if (normalize) {
    result = result / norm(result)
  }
  return (result)
}

re_orient_normals = function(scene_full)
{
  require('Rvcg')
  # align the scene according to the normal median
  
  scene_full = vcgUpdateNormals(scene_full)
  normals = scene_full$normals[1:3,]
  meanvec = apply(normals,1,median)
  
  current_orientation=rbind(meanvec)
  desired_orientation=rbind(c(0,0,1))
  
  scene_full = align_cloud(desired_orientation, current_orientation, scene_full)
  return (scene_full)
}

score_rotation = function(a,ref_vec='x', angle=0)
{
  if (ref_vec=='x') {
    ref_vec = c(1,0,0)
  } else if (ref_vec=='y') {
    ref_vec = c(0,1,0)
  } else {
    cat('not a good value for ref_vec, should be "x" or "y"\n')
    return(-1)
  }
  cuttedv = cut((rotate3d(a, 
                          angle,
                          ref_vec[1], 
                          ref_vec[2], 
                          ref_vec[3]))$vb[3,],
                include.lowest = T, 
                labels = F,
                breaks=seq(min(a$vb[3,],na.rm = T), 
                           max(a$vb[3,],na.rm = T), 
                           le=101))
  return(max(table(cuttedv)))
}

rotate_best_angle = function(a, maxangle, le=21)
{
  best_score = 0
  best_ax = NA
  best_angle = NA
  pb = txtProgressBar(min = 0, max = le*2, style=3); pbi = 0
  for (ax in c('x','y')) {
    for (angle in seq(-maxangle,maxangle,le=le)) {
      score = score_rotation(a, ref_vec=ax, angle)
      if (score > best_score) {
        best_score = score
        best_ax = ax
        best_angle = angle
      }
      pbi = pbi + 1
      setTxtProgressBar(pb,pbi)
    }
  }
  cat('\nThe best angle found was:',best_angle,'around axis',best_ax,'\n')
  return(list(angle=best_angle,ax=best_ax))
}

align_cloud = function(desired_orientation, current_orientation, mesh)
{
  # ensure the right dimensionality of the input
  if (is.null(dim(desired_orientation))) {
    desired_orientation = rbind(desired_orientation)
  }
  if (is.null(dim(current_orientation))) {
    current_orientation = rbind(current_orientation)
  }
  angle = acos( sum(desired_orientation*current_orientation) / ( sqrt(sum(desired_orientation * desired_orientation)) * sqrt(sum(current_orientation * current_orientation)) ) )  
  
  ref_vec = vcrossp(desired_orientation,current_orientation)
  return(rotate3d(mesh, angle, ref_vec[1], ref_vec[2], ref_vec[3]))
}

optimize_horizontal = function(a, maxangle = pi/4, factors = c(1,10,100,1000), le=10,verbose=F)
{
  b = a
  a$normals = NULL
  a$material = NULL
  a$vb[1,] = a$vb[1,] - median(a$vb[1,])
  a$vb[2,] = a$vb[2,] - median(a$vb[2,])
  a$vb[3,] = a$vb[3,] - median(a$vb[3,])
  pb = txtProgressBar(min = 0, max = 2*le*length(factors)*4, style=3); pbi = 0
  for (factor in factors) {
    for (ax in c('x','y','x','y')) {
      for (direction in c(1, -1)) {
        lagg = -1 # we allow score regression only 2 times 
        best_score = max(table(cut(a$vb[3,],
                                   include.lowest = T, 
                                   labels = F,
                                   breaks=seq(min(a$vb[3,],na.rm = T), 
                                              max(a$vb[3,],na.rm = T), 
                                              le=10*factor+1))))
        best_angle = 0
        if (direction == 1) {
          seq_angle = seq(0,maxangle/factor,le=le+1)[2:(le+1)]
        } else {
          seq_angle = seq(0,-maxangle/factor,le=le+1)[2:(le+1)]
        }
        for (angle_i in seq_along(seq_angle)) {
          angle = seq_angle[angle_i]
          pbi = pbi + 1
          setTxtProgressBar(pb,pbi)
          score = max(table(cut((rotate3d(a, 
                                          angle,
                                          ax=='x', 
                                          ax=='y', 
                                          1))$vb[3,],
                                include.lowest = T, 
                                labels = F,
                                breaks=seq(min(a$vb[3,],na.rm = T), 
                                           max(a$vb[3,],na.rm = T), 
                                           le=10*factor+1))))
          if (score > best_score) {
            best_score = score
            best_angle = angle
          } else {
            if (lagg >= 0) {
              break
            } else {
              lagg = lagg + 1
            }
          }
        }
        if (best_angle != 0) {
          if (verbose) {
            cat('performed a rotation of',best_angle,'around axis',ax,'\n\n')
          }
          a = rotate3d(a, best_angle, ax=='x', ax=='y', 0)
          a$vb[1,] = a$vb[1,] - median(a$vb[1,])
          a$vb[2,] = a$vb[2,] - median(a$vb[2,])
          a$vb[3,] = a$vb[3,] - median(a$vb[3,])
        }
        pbi = pbi + (le - angle_i)
        setTxtProgressBar(pb,pbi)
      }
    }
  }  
  b$vb = a$vb
  return (b)
}

do.resample = function(ground, le=11)
{
  coords = apply(ground$vb[1:3,],1,range)
  box = sapply(1:3,function(i) {
    cut(ground$vb[i,],breaks=seq(coords[1,i],coords[2,i],le=le), labels = F)
  })
  ground$vb = ground$vb[,which(!duplicated(box))]
  return(ground)
}

do.scale = function(aligned, area_width=1)
{
  aligned$vb = aligned$vb[c(1,3,2),]
  for (i in 1:3) {
    aligned$vb[i,] = aligned$vb[i,] - min(aligned$vb[i,])
  }
  # reverse z:
  aligned$vb[3,] = max(aligned$vb[3,]) - aligned$vb[3,]
  # find the largest diagonal to perform the last rotation in the ground plane
  le = 3600
  d = 0
  d.a = NA
  convex_hull = aligned
  convex_hull$vb = convex_hull$vb[,chull(aligned$vb[1,],aligned$vb[2,])]
  convex_hull$normals = matrix(NaN,nrow=4,ncol=ncol(convex_hull$vb))
  cat('\nFinding the optimal rotation in the ground (XY) plane:\n')
  pb = txtProgressBar(min = 1, max = le, style=3)
  i = 0
  for (a in seq(0,2*pi,le=le+2)[2:(le+1)]) {
    i = i + 1
    setTxtProgressBar(pb, i)
    rotated = transform3d(convex_hull, rotationMatrix(a,0,0,1))
    p1.i = which.max(rotated$vb[1,])
    p2.i = which.min(rotated$vb[1,])
    p3.i = which.max(rotated$vb[2,])
    p4.i = which.min(rotated$vb[2,])
    d.new = max(c(
      rotated$vb[1,p1.i] - rotated$vb[1,p2.i],
      rotated$vb[2,p3.i] - rotated$vb[2,p4.i]
    ))
    if (d.new > d) {
      d.a = a + pi/4
      d = d.new
    }
  }
  rotated = transform3d(aligned, rotationMatrix(d.a,0,0,1))
  s = d/sqrt(2)
  rotated$vb = rotated$vb / s
  #finally, ensure that all axes are correctly oriented
  rotated$vb = rotated$vb * ((rowMeans(rotated$vb) > 0)*2-1)
  
  rotated$vb = rotated$vb * area_width
  
  return(rotated)
}

do.slice = function(scaled, breaks=9)
{
  cutted = cut(scaled$vb[3,], breaks=breaks,inc=T)
  cuts = sort(as.numeric(na.omit(unique(cutted))))
  slices = lapply(as.numeric(cuts), 
                  function(i) {cbind(scaled$vb[1,which(as.numeric(cutted)==i)],
                                     scaled$vb[2,which(as.numeric(cutted)==i)])})
  return(list(levels(cutted), slices))
}

do.cluster = function(slice, method = 'BIC')
{
  if (method != 'BIC') {
    cat('Only BIC method is implemented for now. exiting.')
    return(NaN)
  }
  require(mclust)
  slice_clust = Mclust(as.matrix(slice), G=1:50)
  return(slice_clust$classification)
}

do.plot.classif = function(slice, classif, colmap=rainbow, ...) {
  plot(slice, col=colmap(length(unique(classif)))[classif], ...)
}

do.circle = function(slice, classif, maxshoot = 0.1) {
  p = t(sapply(1:length(unique(classif)), function(i) {
    ind = which(classif == i)
    mx = mean(slice[ind,1])
    my = mean(slice[ind,2])
    rmedian = median(sqrt((slice[ind,1]-mx)^2 + (slice[ind,2]-my)^2))
    rmean = mean(sqrt((slice[ind,1]-mx)^2 + (slice[ind,2]-my)^2))
    r90 = quantile(sqrt((slice[ind,1]-mx)^2 + (slice[ind,2]-my)^2),probs = 0.9)
    rmax = max(sqrt((slice[ind,1]-mx)^2 + (slice[ind,2]-my)^2))
    return (c(mx, my, rmedian, rmean, r90, rmax))
  }))
  # clean the cluster reconstruction
  p = p[complete.cases(p),]
  if (any(p[,3] > maxshoot)) {
    p = p[-which(p[,3]>maxshoot),]
  }
  return(p)
}



generate_centers = function()
{
  n = sample.int(n = 200,size = 1)
  return(matrix(runif(2*n),ncol=2,byrow = T))
}

cost = function(slice, centers, diameters=NA)
{
  require('FNN')
  clust = get.knnx(centers, slice, k=1)
  # first try at cost function include looking at:
  # - number of cluster (the lower the better)
  # - average distance of points to the nearest cluster (the lower the better)
  c1 = nrow(centers)
  c2 = median(sapply(clust$nn.index,function(i) {mean(clust$nn.dist[clust$nn.index==i])} ))
  return(c(c1,c2))
}


find_solution = function(slice,tries=100)
{
  require('FNN')
  require('emoa')
  sol_costs = NULL
  sols = NULL
  pbi = 0; pb = txtProgressBar(min = 1, max = tries, style=3)
  for (i in 1:tries) {
    pbi = pbi + 1; setTxtProgressBar(pb, pbi)
    
    # we pick a set of points from the original set
    sols[[i]] = matrix(slice[sample.int(nrow(slice),sample(50:300,size = 1)),],ncol=2,byrow = F)
    
    # and we remove points that are less than 5 mm apart (0.005)
    tst = get.knn(sols[[i]], k=1)
    good = which(tst$nn.dist < 0.005)
    sols[[i]] = sols[[i]][-unique(t(apply(cbind(good,tst$nn.index[good]),1,sort)))[,1], ]
    
    c1 = nrow(sols[[i]])
    
    clust = get.knnx(sols[[i]], slice, k=1)
    c2 = mean(sapply(1:c1,function(iii) {quantile(clust$nn.dist[clust$nn.index==iii], probs=0.9)} ))
    
    scost = c(c1,c2)
    sol_costs = rbind(sol_costs,scost)
  }
  poptimal = !is_dominated(t(sol_costs))
  plot(sol_costs,col=rainbow(2)[1+poptimal])
  optimal_sols = list()
  optimal_costs = list()
  for (i in 1:tries) {
    if (poptimal[i]) {
      optimal_sols[[length(optimal_sols)+1]] = sols[[i]]
      optimal_costs[[length(optimal_costs)+1]] = sol_costs[i,]
    }
  }
  return(list(optimal_sols,optimal_costs))
}


guess_optimal_bandwidth = function(slices, resample.n = Inf, iter = 100, cores=NA)
{
  if (is.na(cores)) {
    cores = getCoresNumber()
  }
  require('doParallel')
  registerDoParallel(cores=cores)
  require('LPCM')
  ns = sapply(1:length(slices),function(i.slice)nrow(slices[[i.slice]]))
  hs = seq(0.0025,0.0075,by=0.00025)
  mss =
    foreach (h = hs,
             .packages = 'LPCM') %dopar% 
             {   i.slice=1
             return(c(cluster.slice=i.slice,
                      h=h,
                      ms(slices[[i.slice]], h=h, thr=0.0001, 
                         subset = sample.int(ns[i.slice],min(resample.n,ns[i.slice])),
                         iter=iter, scaled=F,
                         plotms=F)))
             }
  
  # re-order the list according to the bandwidth
  for (i in length(mss):1) {
    if (is.null(mss[[i]])) {
      mss[i] = NULL
    }
  }
  mss = mss[order(unlist(lapply(mss,function(m)m$h)))]
  
  mmm2 = assign.sector.width(mss,slices,s = 8)
  
  symm_ratios = unlist(lapply(mmm2,function(m){
    if (nrow(m$cluster.ws) == 1) {
      diams = m$cluster.ws[1:4] + m$cluster.ws[5:8]
      r = diff(range(diams)) / max(diams,na.rm=T)
    } else {
      diams = m$cluster.ws[,1:4] + m$cluster.ws[,5:8]
      r = median(c(diff(apply(diams,1,range))) / apply(diams,1,max),na.rm=T)     
    }
    return(r)
  }))
  
  bh = which.max(symm_ratios)
  
  plot(hs,symm_ratios,
       las=1,
       xlab='bandwidth',
       ylab='Mean root symmetry ratio',
       main='Bandwidth yielding optimal root circularity\n(highest symmetry ratio)'
  )
  
  hss = seq(min(hs),max(hs),by=0.00001)
  
  con = sapply(hss,
               function(x){
                 sum(sapply(1:length(hs), function(i){
                   symm_ratios[i]*dnorm(x = x-hs[i], sd = 0.00025)
                 }))
               })
  hsbh = mean(hss[order(con, decreasing = T)[1:10]])
  
  abline(v=hsbh,col='red')
  text(hsbh,1,'Max likelihood',col='red',xpd=T,adj=0)
  
  par(new=TRUE)
  
  plot(hss, con, xlab="", ylab="", 
       axes=FALSE, type="l", col="orange")
  text(mean(hss),mean(con),"Likelihood",col="orange") 
  
  return(hsbh)
}


approximate = function(slices, breaks=seq(0.1,0.6,by=0.02), resample.n = Inf, iter = 100, h = 0.005, cores=NA)
{
  ns = sapply(1:length(slices),function(i.slice)nrow(slices[[i.slice]]))
  maxnb = length(slices)
  
  if (is.na(cores)) {
    cores = getCoresNumber()
  }
  
  if (cores==1) {
    # that's the mono-thread version
    require('LPCM')
    require('tcltk')
    pb = tkProgressBar('Progress', min=0, max=maxnb)
    nb = 0
    
    cat(paste('\nApproximation took',system.time({
      mss =
        lapply (c(sample(1:length(slices))), function(i.slice) {
          m = c(cluster.slice=i.slice,
                cluster.hm=breaks[i.slice],
                cluster.hM=breaks[i.slice+1],
                ms(slices[[i.slice]], h=h, thr=0.0001, 
                   subset = sample.int(ns[i.slice],min(resample.n,ns[i.slice])),
                   iter=iter, scaled=F,
                   plotms=F))
          nb = nb + 1;
          setTkProgressBar(pb, nb)
          return(m)
        })
    })[3],'seconds'))
    
  } else {
    require('doParallel')
    registerDoParallel(cores=cores)
    
    token = as.character(round(runif(1,min=10000,max=99999)))
    tmpdir = tempdir()
    
    cat(paste('\nApproximation took',system.time({
      mss =
        foreach (i.slice = c(0,sample(1:length(slices))),
                 .packages = 'LPCM') %dopar% 
                 {
                   if (i.slice == 0) {
                     ## MASTER
                     require('tcltk')
                     # require('tcltk2')
                     pb = tkProgressBar('Progress', min=0, max=maxnb-cores)
                     nb = 0
                     while (nb < maxnb-1) {
                       nb = length(list.files(tmpdir,pattern=paste0('^',token)))
                       nb2 = max(0,nb-cores)
                       Sys.sleep(1)
                       setTkProgressBar(pb, nb2)
                     }
                     close(pb)
                     Sys.sleep(1)
                     return(NULL)
                   } else {
                     file.create(tempfile(pattern=token, tmpdir = tmpdir))
                     return(c(cluster.slice=i.slice,
                              cluster.hm=breaks[i.slice],
                              cluster.hM=breaks[i.slice+1],
                              ms(slices[[i.slice]], h=h, thr=0.0001, 
                                 subset = sample.int(ns[i.slice],min(resample.n,ns[i.slice])),
                                 iter=iter, scaled=F,
                                 plotms=F)))
                   }
                 }
    })[3],'seconds'))
  }
  
  # re-order the list according to the slice numbers
  for (i in length(mss):1) {
    if (is.null(mss[[i]])) {
      mss[i] = NULL
    }
  }
  mss = mss[order(unlist(lapply(mss,function(m)m$cluster.slice)))]
  return(mss)
}


plotcyl = function(mss,breaks=seq(0.1,0.6,by=0.02),maxn=NULL,
                   colmap = rainbow, offset=0)
{
  # maxn is the number of roots (i.e. the number of unique colors)
  # offset is an additional translation over the X axis
  if (!is.null(maxn)) {
    colmap = sample(colmap(maxn))
  } else {
    colmap = colmap
  }
  cat(paste('\n3D plot took',system.time({
    for (ii in 1:(length(mss))) {
      ii_breaks = ii
      m = mss[[ii]]
      n = nrow(m$cluster.center)
      hm = breaks[ii_breaks]
      hM = breaks[ii_breaks+1]
      for (i in 1:n) {
        p = m$cluster.center[i,]
        if (!is.null(m$cluster.id) & !is.null(maxn)) {
          col = colmap[m$cluster.id[i]]
        } else {
          col=colmap(n)[i]
        }
        if (!is.null(m$cluster.ws)) {
          s = ncol(m$cluster.ws) # maximal number of cluster
          middle_s = s/2-1
          section=rep(m$cluster.ws[i,],each=2)[c(2:(2*s),1)]*
            cbind(cos(seq(0, 2 * pi, len = s + 1)[-1]),
                  sin( seq(0, 2 * pi, len = s + 1)[-1]))[rep(c(middle_s:1,s:(middle_s+1)),each=2),] # sector parameterization
          radius = 1
        } else if (!is.null(m$cluster.dia)) {
          section = NULL
          radius = m$cluster.dia[i]
        } else {
          section = NULL
          radius = 0.01
        }
        plot3d(cylinder3d(rbind(c(p[1]+offset, p[2], hm),
                                c(p[1]+offset, p[2], hM)),
                          radius=radius, 
                          section=section,
                          close=-2), 
               col=col, add=T)
      }
    }
  })[3],'seconds'))
  
}


merge.clusters = function(mss, threshold = 0.01)
{
  require('FNN')
  unique.id = 1
  new.mss = NULL
  for (i.slice in 1:length(mss)) {
    m = c(cluster.id = NA, mss[[i.slice]])
    n = nrow(m$cluster.center)
    m$cluster.id = rep(NA,n)
    if (i.slice > 1) {
      mp = new.mss[[i.slice-1]]
      correspondance = get.knnx(mp$cluster.center, m$cluster.center, k=1)
      cluster.close = correspondance$nn.dist < threshold
      m$cluster.id[cluster.close] = mp$cluster.id[correspondance$nn.index[cluster.close]]
      tofill = sum(m$cluster.id[cluster.close])
      n = sum(!cluster.close)
      # browser()
      if (n > 0) {
        m$cluster.id[!cluster.close] = unique.id:(unique.id+n-1)
      }
    } else {
      if (n > 0) {
        m$cluster.id = unique.id:(unique.id+n)
      }
    }
    unique.id = unique.id + n
    new.mss[[i.slice]] = m
  }
  return (list(mss=new.mss,unique.id=unique.id))
}

get_sector = function(P, maxs=8)
{
  # this function returns the sector number of points P, when discretizing over maxs sectors
  s = floor(acos(P[,1]/sqrt(rowSums(P^2)))/pi*maxs/2)
  s[which(P[,2]<0)] = (-s[which(P[,2]<0)]-1) %% maxs
  return(s)
}

get_sector_distance = function(P, maxs=8)
{
  # this function returns the triangular sector distance, when discretizing over maxs sectors
  a = sapply(1:nrow(P),function(i) {
    aaa = acos(P[i,1]/sqrt(sum(P[i,]^2)))
    aaa = ifelse(P[i,2] < 0, 2*pi-aaa, aaa)
    bbb = floor(acos(P[i,1]/sqrt(sum(P[i,]^2)))/pi*maxs/2)
    bbb = ifelse(P[i,2] < 0, ((-bbb) %% maxs) *pi*2/maxs, bbb * pi*2/maxs)
    ccc = ceiling(acos(P[i,1]/sqrt(sum(P[i,]^2)))/pi*maxs/2)
    ccc = ifelse(P[i,2] < 0, ((-ccc) %% maxs) *pi*2/maxs, ccc * pi*2/maxs)
    # there are easier ways...
    abs(c(aaa-(bbb+ccc)/2))
  })
  mlength = cos(a)*sqrt(rowSums(P^2))
  return(abs(mlength/cos(pi*2/maxs/2)))
}


assign.width = function(mss, slices, threshold=0.06)
{
  require('FNN')
  unique.id = 1
  new.mss = NULL
  pb = txtProgressBar(min = 1, max = length(mss), style=3); pbi = 0
  for (i.slice in 1:length(mss)) {
    pbi = pbi + 1; setTxtProgressBar(pb, pbi)
    m = c(cluster.dia = NA, mss[[i.slice]])
    n = nrow(m$cluster.center)
    m$cluster.dia = rep(NA,n)
    for (i.cluster in 1:n) {
      correspondance = get.knnx(m$cluster.center, slices[[i.slice]], k=1)
      cluster.close = correspondance$nn.dist < threshold
      correspondance$nn.dist[!cluster.close] = NA
      m$cluster.dia = sapply(1:n, function(i)quantile(correspondance$nn.dist[correspondance$nn.index==i],probs=0.95,na.rm=T))
    }
    new.mss[[i.slice]] = m
  }
  return (new.mss)
}

assign.sector.width = function(mss, slices, s=8, threshold=0.06, fun=max, ...)
{
  require('FNN')
  unique.id = 1
  new.mss = NULL
  pb = txtProgressBar(min = 1, max = length(mss), style=3); pbi = 0
  for (i.slice in 1:length(mss)) {
    pbi = pbi + 1; setTxtProgressBar(pb, pbi)
    m = c(cluster.ws = NA, mss[[i.slice]])
    n = nrow(m$cluster.center)
    m$cluster.ws = matrix(NA,nrow=n,ncol=s)
    correspondance = get.knnx(m$cluster.center, slices[[i.slice]], k=1)
    # this gives the sector for all points, relative to their cluster:
    sectors = get_sector(slices[[i.slice]] - m$cluster.center[correspondance$nn.index,],maxs=s)
    # browser()
    sectors_dist = get_sector_distance(
      slices[[i.slice]] - m$cluster.center[correspondance$nn.index,],maxs=s)
    for (i.cluster in 1:n) {
      # for each cluster within a slice, get the quantiles per sector
      m$cluster.ws[i.cluster,] = sapply(1:s, function(i){
        iii = which((sectors==i-1) & (correspondance$nn.index == i.cluster))
        indicator = max(sectors_dist[iii],na.rm=T, ...)
        if (is.infinite(indicator) | is.na(indicator)) {
          return (0)
        } else {
          return (indicator)
        }
      })
    }
    new.mss[[i.slice]] = m
  }
  return (new.mss)
}


