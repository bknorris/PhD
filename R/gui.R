source('functions.R')

library(gWidgets)
library(gWidgetsRGtk2)
require(RGtk2)
options("guiToolkit"="RGtk2")

myplot3d = function(...){
  if (length(rgl.dev.list()))
    rgl.set(rgl.dev.list())
  par3d(windowRect=c(10,10,800,600))
  plot3d(...,xlab='X',ylab='Y',zlab='Z')
  rgl.pop("lights")
  light3d(specular="white") 
}


read_ply = function(h, var, ...)
{
  tryCatch(
    {
      f = svalue(h$obj)
      if (substring(f,1,1) %in% c('"',"'")) {
        # remove quotes if present
        f = substring(f,2,nchar(f)-1)
      }
      assign(var, NULL, envir = globalenv())
      assign(var, read.ply(f, ShowSpecimen = F), envir = globalenv())
    },
    error = function(e) {
      gmessage('Could not open file! Is it an ASCII ply file? (BINARY ply files can not be opened)')
    }
  )
}




panel1 = function()
{
  win <- gwindow("3D Pneumatophores   -   Data Import")
  
  allwin = ggroup(cont=win, horizontal=T)
  
  BigDataGroup <- ggroup(cont=allwin, horizontal=F)
  DataGroup <- gframe("Data", container = BigDataGroup, horizontal=FALSE)
  
  grp.file <- ggroup(horizontal=FALSE, container = DataGroup)
  lbl.file <- glabel("3D model: ", container = grp.file)
  browse.file <- gfilebrowse(text    = "Select ASCII ply file",
                             container = grp.file,
                             filter = list("ply files" = list(patterns = c("*.ply"))),
                             handler = function(h, ...){
                               enabled(browse.file) = F;
                               enabled(grp.file) = F;
                               read_ply(h, var='tst', ...);
                               enabled(txt_data_frame_name) = T;
                               enabled(chk_ground) = T;
                               enabled(btn_next) = T;
                             })
  
  grp.file2 <- ggroup(horizontal=FALSE, container = DataGroup)
  lbl.file2 <- glabel("Horizontal optimization: ", container = grp.file2)
  chk_ground <- gcheckbox(
    text      = "Skip alignment with ground",
    handler = function(h,...){
    },
    container = grp.file2
  )
  enabled(chk_ground) = F
  
  SizeGroup <- gframe("Size", container = BigDataGroup, horizontal=FALSE)
  grp_name <- ggroup(horizontal=FALSE, container = SizeGroup)
  lbl_data_frame_name <- glabel(
    "Enter the size of the area (in meters) : ",
    container = grp_name
  )
  txt_data_frame_name <- gedit("1.0", container = grp_name)
  enabled(txt_data_frame_name) = F
  
  RightGroup <- ggroup(cont=allwin, horizontal=F)
  addSpring(RightGroup)
  
  btn_next <- gbutton(
    text      = "Go to next step",
    container = RightGroup,
    handler   = function(h, ...)
    {
      enabled(btn_next) = F
      assign('area_width', as.numeric(svalue(txt_data_frame_name)), envir = globalenv())
      
      if (svalue(chk_ground) == 1) {
        aligned = re_orient_normals(tst, plot=F)
        aligned = optimize_horizontal(aligned)
      } else {
        aligned = tst
      }
      
      
      assign('scaled', do.scale(aligned, area_width), envir = globalenv())
      panel3()
      dispose(win)
    }
  )
  enabled(btn_next) = F
  
}

panel3 = function() {
  win <- gwindow("3D Pneumatophores   -   Ground Removal")
  
  grp_name <- ggroup(container = win,horizontal=F)

  lbl_data_frame_name2 <- glabel(
    "Offset to align the ground with 0 (in meters): ",
    container = grp_name
  )
  txt_data_frame_name2 <- gedit("0", container = grp_name)

  lbl_data_frame_name <- glabel(
    "Height of the first slice (relative to the ground, in meters): ",
    container = grp_name
  )
  txt_data_frame_name <- gedit("0.1", container = grp_name)

  btn_flip <- gbutton(
    text      = "Flip it yeah!",
    container = grp_name,
    handler   = function(h, ...)
    {
      scaled$vb[3,] = -scaled$vb[3,]
      assign('scaled',scaled, envir = globalenv())    }
  )

    
  btn_v <- gbutton(
    text      = "Update",
    container = win,
    handler   = function(h, ...)
    {
      tryCatch(
        {
          scaled$vb[3,] = scaled$vb[3,] - as.numeric(svalue(txt_data_frame_name2))
          assign('scaled',scaled, envir = globalenv())
          svalue(txt_data_frame_name2) = '0'
          assign('baseline',as.numeric(svalue(txt_data_frame_name)), envir = globalenv())
          cat(paste('baseline is',baseline,'\n'))
          thr = 0.3
          hist(scaled$vb[3,scaled$vb[3,]<thr],breaks=seq(min(scaled$vb[3,]),thr+0.01,by=0.01), include.lowest = T, xlab='height', main='Distribution of point heights')
          abline(v=baseline,col='red')
          myplot3d(t(scaled$vb))
          planes3d(a=0,b=0,c=-1,d=baseline,col='red')
          planes3d(a=0,b=0,c=1,d=0,col='green')
          enabled(btn_next) = T
        },
        error = function(e) {
          enabled(btn_next) = F
          gmessage('Could not read the relative height of the first slice or the ground level. Are they numerical value?')
        })
    }
  )
  
  grp_next <- ggroup(container = win)
  btn_next <- gbutton(
    text      = "Go to next step",
    container = grp_next,
    handler   = function(h, ...)
    {
      panel4()
      dispose(win)
    }
  )
  enabled(btn_next) = F
  
}






panel4 = function() {
  win <- gwindow("3D Pneumatophores   -   Cylindrical Approximation")
  
  allwin = ggroup(cont=win, horizontal=T)
  
  BigDataGroup <- ggroup(cont=allwin, horizontal=F)
  DataGroup <- gframe("Clustering Parameters", container = BigDataGroup, horizontal=FALSE)
  
  grp_name <- ggroup(horizontal=FALSE, container = DataGroup)
  lbl_data_frame_name <- glabel(
    "Enter the height of a slice (in meters)\nValues like 0.005, 0.01 or 0.02 make sense:",
    container = grp_name
  )
  txt_data_frame_name <- gedit("0.01", container = grp_name)
  
  lbl_data_frame_name2 <- glabel(
    "Enter the max number of points within a slice (lower = faster, but less accurate)\nValues like 1000, 10000 or Inf make sense:",
    container = grp_name
  )
  txt_data_frame_name2 <- gedit("Inf", container = grp_name)
  
  lbl_data_frame_name3 <- glabel(
    "Enter the number of iterations for each step of approximation (lower = faster, but less accurate)\nValues like 100, 200 or 300 make sense:",
    container = grp_name
  )
  txt_data_frame_name3 <- gedit("300", container = grp_name)
  
  grp_bandwidth <- ggroup(cont=grp_name, horizontal=T)
  
  lbl_data_frame_name4 <- glabel(
    "Enter the bandwidth of the mean-shift kernel (lower = more small clusters)\nValues like 0.0025, 0.005 or 0.0075 make sense:",
    container = grp_bandwidth
  )
  txt_data_frame_name4 <- gedit("0.005", container = grp_bandwidth)
  btn_bandwidth <- gbutton(
    text      = "Guess optimal",
    container = grp_bandwidth,
    handler   = function(h, ...)
    {
      tryCatch(
        {
          step = as.numeric(svalue(txt_data_frame_name))
          cat(paste('step is',step,'\n'))
        },
        error = function(e) {
          step = 0
          gmessage('Could not read the slice width. Is it a numerical value?')
        })
      tryCatch(
        {
          resample.n = as.numeric(svalue(txt_data_frame_name2))
          cat(paste('resample.n is',resample.n,'\n'))
        },
        error = function(e) {
          resample.n = 0
          gmessage('Could not read the max number of points. Is it a numerical value (or, "Inf")?')
        })
      tryCatch(
        {
          iter = as.numeric(svalue(txt_data_frame_name3))
          cat(paste('iter is',iter,'\n'))
        },
        error = function(e) {
          iter = 0
          gmessage('Could not read the iter number. Is it a numerical value?')
        })
      if ((step != 0) & (iter != 0) & (resample.n != 0)) {
#         tryCatch(
#           {
            assign('resample.n', resample.n, envir = globalenv())
            assign('iter', iter, envir = globalenv())
            assign('step', step, envir = globalenv())
            assign('breaks', seq(baseline,max(scaled$vb[3,])+step,by=step), envir = globalenv())
            tmp_list = do.slice(scaled, breaks)
            assign('breaks_label', tmp_list[[1]], envir = globalenv())
            assign('slices', tmp_list[[2]], envir = globalenv())
            svalue(txt_data_frame_name4) = guess_optimal_bandwidth(slices,resample.n = resample.n, iter = iter)
#           }, error = function(e) {
#             enabled(txt_data_frame_name4) = T
#             gmessage('Could not perform automatic bandwidth detection. Falling back into debug mode..')
#             browser()
#           }
        # )
      }
    }
  )
  
  RightGroup <- ggroup(cont=allwin, horizontal=F)
  addSpring(RightGroup)
  
  lbl_data_frame_name_lab <- glabel(
    "Expect processing time up to 30 minutes",
    container = RightGroup
  )
  lbl_data_frame_name_lab2 <- glabel(
    "using the default settings...",
    container = RightGroup
  )
  
  btn_scale <- gbutton(
    text      = "Approximate",
    container = RightGroup,
    handler   = function(h, ...)
    {
      tryCatch(
        {
          step = as.numeric(svalue(txt_data_frame_name))
          cat(paste('step is',step,'\n'))
        },
        error = function(e) {
          step = 0
          gmessage('Could not read the slice width. Is it a numerical value?')
        })
      tryCatch(
        {
          resample.n = as.numeric(svalue(txt_data_frame_name2))
          cat(paste('resample.n is',resample.n,'\n'))
        },
        error = function(e) {
          resample.n = 0
          gmessage('Could not read the max number of points. Is it a numerical value (or, "Inf")?')
        })
      tryCatch(
        {
          iter = as.numeric(svalue(txt_data_frame_name3))
          cat(paste('iter is',iter,'\n'))
        },
        error = function(e) {
          iter = 0
          gmessage('Could not read the iter number. Is it a numerical value?')
        })
      tryCatch(
        {
          b = as.numeric(svalue(txt_data_frame_name4))
          cat(paste('bandwidth is',b,'\n'))
        },
        error = function(e) {
          b = 0
          gmessage('Could not read the bandwidth parameter. Is it a numerical value?')
        })
      if ((step != 0) & (iter != 0) & (resample.n != 0) & (b != 0)) {
        tryCatch(
          {
            assign('resample.n', resample.n, envir = globalenv())
            assign('iter', iter, envir = globalenv())
            assign('step', step, envir = globalenv())
            assign('b', b, envir = globalenv())
            assign('breaks', seq(baseline,max(scaled$vb[3,])+step,by=step), envir = globalenv())
            tmp_list = do.slice(scaled, breaks)
            assign('breaks_label', tmp_list[[1]], envir = globalenv())
            assign('slices', tmp_list[[2]], envir = globalenv())
            assign('mss',
                   approximate(slices,breaks=breaks,iter=iter,resample.n=resample.n,h=b,cores=8),
                   envir = globalenv())
            tmp_list = merge.clusters(mss, threshold=0.1); new.mss=tmp_list[[1]]; 
            assign('unique.n',tmp_list[[2]], envir = globalenv())
            s = 8
            new.mss2 = assign.sector.width(new.mss, slices, threshold=0.06, s=s)
            assign('mss',new.mss2, envir = globalenv())
            
            enabled(btn_upload2) = T
            
            myplot3d(t(scaled$vb));
            plotcyl(mss,breaks=breaks,maxn=unique.n)
            
            # a few statistics
            pneu_nb = sapply(1:length(mss), function(i.slice){
              nrow(mss[[i.slice]]$cluster.center)
            })
            mean_rad = sapply(1:length(mss), function(i.slice){
              mean(mss[[i.slice]]$cluster.ws)
            })
            sem_rad = sapply(1:length(mss), function(i.slice){
              sqrt(var(c(mss[[i.slice]]$cluster.ws))/pneu_nb[i.slice]/s)
            })
            mean_diam = sapply(1:length(mss), function(i.slice){
              mean(c(mss[[i.slice]]$cluster.ws[,1:(s/2)]+mss[[i.slice]]$cluster.ws[,(1+s/2):s]))
            })
            sem_diam = sapply(1:length(mss), function(i.slice){
              sqrt(var(c(
                mss[[i.slice]]$cluster.ws[,1:(s/2)]+mss[[i.slice]]$cluster.ws[,(1+s/2):s]
                ))/pneu_nb[i.slice]/(s/2))
            })

            min_diam = sapply(1:length(mss), function(i.slice){
              if(nrow(mss[[i.slice]]$cluster.ws) == 1) {
                min(mss[[i.slice]]$cluster.ws[1:(s/2)]+mss[[i.slice]]$cluster.ws[(1+s/2):s])
              } else {
                mean(apply(mss[[i.slice]]$cluster.ws[,1:(s/2)]+mss[[i.slice]]$cluster.ws[,(1+s/2):s],1,min))
              }
            })
            max_diam = sapply(1:length(mss), function(i.slice){
              if(nrow(mss[[i.slice]]$cluster.ws) == 1) {
                max(mss[[i.slice]]$cluster.ws[1:(s/2)]+mss[[i.slice]]$cluster.ws[(1+s/2):s])
              } else {
                mean(apply(mss[[i.slice]]$cluster.ws[,1:(s/2)]+mss[[i.slice]]$cluster.ws[,(1+s/2):s],1,max))
              }
            })
            median_diam = sapply(1:length(mss), function(i.slice){
              if(nrow(mss[[i.slice]]$cluster.ws) == 1) {
                median(mss[[i.slice]]$cluster.ws[1:(s/2)]+mss[[i.slice]]$cluster.ws[(1+s/2):s])
              } else {
                mean(apply(mss[[i.slice]]$cluster.ws[,1:(s/2)]+mss[[i.slice]]$cluster.ws[,(1+s/2):s],1,median))
              }
            })
            
            ex_p = c(T,F)
            bp = barplot(mean_rad*100, horiz=T,main='Mean radius as a function of height',xlab='Radius (cm)')
            require('plotrix')
            plotCI(mean_rad*100,bp,uiw=sem_rad*100, err='x',sfrac=0.005,add=T,pch=NA,xpd=T)
            text(x=par('usr')[1],bp[ex_p],breaks_label[ex_p],adj=1, cex=0.75,xpd=T)
            text(x=par('usr')[1],bp[length(bp)]+2*diff(bp[1:2]),'Height (m)',adj=1, cex=1,xpd=T)
            bp = barplot(mean_diam*100, horiz=T,main='Mean diameter as a function of height',xlab='Diameter (cm)')
            plotCI(mean_diam*100,bp,uiw=sem_diam*100, err='x',sfrac=0.005,add=T,pch=NA,xpd=T)
            text(x=par('usr')[1],bp[ex_p],breaks_label[ex_p],adj=1, cex=0.75,xpd=T)
            text(x=par('usr')[1],bp[length(bp)]+2*diff(bp[1:2]),'Height (m)',adj=1, cex=1,xpd=T)
            bp = barplot(pneu_nb, horiz=T,main='Number of pneumatophores per height',xlab='Count')
            text(x=par('usr')[1],bp[ex_p],breaks_label[ex_p],adj=1, cex=0.75,xpd=T)
            text(x=par('usr')[1],bp[length(bp)]+2*diff(bp[1:2]),'Height (m)',adj=1, cex=1,xpd=T)
            frontal = pneu_nb * mean_diam
            sem_frontal = sem_diam * pneu_nb 
            tot_frontal = sum(frontal)
            bp = barplot(frontal, horiz=T,
                         main=paste0('Mean +/- SEM frontal area densities\n(sum from ',breaks[1],' cm up to ',breaks[length(breaks)],' cm: ',round(tot_frontal,digits=3),' per meter)'),
                         xlab='Frontal area density (per meter and per slice)')
            plotCI(frontal,bp,uiw=sem_frontal, err='x',sfrac=0.005,add=T,pch=NA,xpd=T)
            text(x=par('usr')[1],bp[ex_p],breaks_label[ex_p],adj=1, cex=0.75,xpd=T)
            text(x=par('usr')[1],bp[length(bp)]+2*diff(bp[1:2]),'Height (m)',adj=1, cex=1,xpd=T)
            minfrontal = pneu_nb * min_diam
            tot_minfrontal = sum(minfrontal)
            bp = barplot(minfrontal, horiz=T,
                         main=paste0('Min frontal area densities\n(sum from ',breaks[1],' cm up to ',breaks[length(breaks)],' cm: ',round(tot_minfrontal,digits=3),' per meter)'),
                         xlab='Frontal area density (per meter and per slice)')
            text(x=par('usr')[1],bp[ex_p],breaks_label[ex_p],adj=1, cex=0.75,xpd=T)
            text(x=par('usr')[1],bp[length(bp)]+2*diff(bp[1:2]),'Height (m)',adj=1, cex=1,xpd=T)
            maxfrontal = pneu_nb * max_diam
            tot_maxfrontal = sum(maxfrontal)
            bp = barplot(maxfrontal, horiz=T,
                         main=paste0('Max frontal area densities\n(sum from ',breaks[1],' cm up to ',breaks[length(breaks)],' cm: ',round(tot_maxfrontal,digits=3),' per meter)'),
                         xlab='Frontal area density (per meter and per slice)')
            text(x=par('usr')[1],bp[ex_p],breaks_label[ex_p],adj=1, cex=0.75,xpd=T)
            text(x=par('usr')[1],bp[length(bp)]+2*diff(bp[1:2]),'Height (m)',adj=1, cex=1,xpd=T)
            medianfrontal = pneu_nb * median_diam
            tot_medianfrontal = sum(medianfrontal)
            bp = barplot(medianfrontal, horiz=T,
                         main=paste0('Median frontal area densities\n(sum from ',breaks[1],' cm up to ',breaks[length(breaks)],' cm: ',round(tot_medianfrontal,digits=3),' per meter)'),
                         xlab='Frontal area density (per meter and per slice)')
            text(x=par('usr')[1],bp[ex_p],breaks_label[ex_p],adj=1, cex=0.75,xpd=T)
            text(x=par('usr')[1],bp[length(bp)]+2*diff(bp[1:2]),'Height (m)',adj=1, cex=1,xpd=T)

            ## assign the newly computed statistics into the main environment
            assign('pneu_nb', pneu_nb,  envir = globalenv())
            assign('mean_rad', mean_rad,  envir = globalenv())
            assign('sem_rad', sem_rad,  envir = globalenv())
            assign('mean_diam', mean_diam,  envir = globalenv())
            assign('sem_diam', sem_diam,  envir = globalenv())
            assign('min_diam', min_diam,  envir = globalenv())
            assign('max_diam', max_diam,  envir = globalenv())
            assign('median_diam', median_diam,  envir = globalenv())
            assign('frontal', frontal,  envir = globalenv())
            assign('sem_frontal', sem_frontal,  envir = globalenv())
            assign('minfrontal', minfrontal,  envir = globalenv())
            assign('maxfrontal', maxfrontal,  envir = globalenv())
            
            cat('\n\nSTART OUTPUT TO CSV FILE\n\nhm,hM,nb,aveDIA,semDIA,minDIA,maxDIA,aveFDA,semFDA,minFDA,maxFDA\n')
            for (i.slice in 1:length(mss)){
              m = mss[[i.slice]]
              cat(paste0(m$cluster.hm,
                         ',',
                         m$cluster.hM,
                         ',',
                         nrow(m$cluster.ws),
                         ',',
                         mean_diam[i.slice],
                         ',',
                         sem_diam[i.slice],
                         ',',
                         min_diam[i.slice],
                         ',',
                         max_diam[i.slice],
                         ',',
                         frontal[i.slice],
                         ',',
                         sem_frontal[i.slice],
                         ',',
                         minfrontal[i.slice],
                         ',',
                         maxfrontal[i.slice],
                         '\n'
              ))
            }
            cat('\n\nEND OUTPUT TO CSV FILE\n\n')
            
          },
          error = function(e) {
            gmessage('Could not perform the approximation :(')
          })
      }
    }
  )
  
  btn_upload2 <- gbutton(
    text      = "Save Workspace",
    container = RightGroup,
    handler   = function(h, ...)
    {
      gfile(
        text    = "Where to save the workspace?",
        type    = "open",
        # action  = function(){},
        handler = function(h, ...)
        {
          tryCatch(
            {
              save.image(file=h$file)
            },
            error = function(e) {
              gmessage('Could not save the workspace!')
            }
          )
        },
        filter = list(
          "RData files" = list(patterns = c("*.RData"))
        )
      )
    }
  )
  enabled(btn_upload2) = F
  
}


tst = NULL
ground = NULL
scaled = NULL
area_width = 1.0
baseline = 0
mss=NULL
breaks=NULL
breaks_label=NULL
iter=NULL
resample.n=NULL
step=NULL
slices = NULL
b=NULL
slices=NULL
pneu_nb = NULL
mean_rad = NULL
sem_rad = NULL
mean_diam = NULL
sem_diam = NULL
min_diam = NULL
max_diam = NULL
median_diam = NULL
frontal = NULL
sem_frontal = NULL
minfrontal = NULL
maxfrontal = NULL
panel1()
