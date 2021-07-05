  library('Cairo')
	library('R2HTML')
  library('waveslim')
	library('brainwaver')
######### modif: igraph or igrahp0
	library('igraph') # commence a 1
	library('ape')
  library('corrplot')
  source('Figure.R')
  
  compute_Graph_annotated <- function(rs_path,
                                      atlas_file,
                                      graph_path,
                                      file_coord,
                                      graphs = TRUE, 
                                      regions_selected=NULL,
                                      percentage_selected=4, 
                                      num.levels,
                                      num.nb.edges,
                                      length_time_series=NULL)
  {

    templBaseName <- file_path_sans_ext(basename(atlas_file))
    name.ts <- file_path_sans_ext(basename(atlas_file))
    
    length.nb.edges <- length(num.nb.edges)
    serie_temp <- paste(graph_path,'/Functional/data/',templBaseName,'/',name.ts,'_ts.txt',sep='')
    data.roi <- read.table(serie_temp)
    proc.length <- dim(data.roi)[1]
    
    if(is.null(length_time_series)){
      length_time_series<-1:proc.length
    }
    
    data.roi<-data.roi[length_time_series,]
    proc.length<-dim(data.roi)[1]
    original.data<-read.table(serie_temp)
    n.regions<-dim(data.roi)[2] # /!\ a changer selon le nb de regions souhaitees
    cat("nb_pt_tps ", proc.length, " nb_regions ", n.regions,"\n")
    
    if(is.null(regions_selected)){
      regions_selected = seq_len(n.regions)
    }
    
    ######### DOSSIER DE SORTIE
    dir.create(file.path(graph_path, "Graph_Measures"), showWarnings = FALSE)
    dir.create(file.path(paste(graph_path,'Graph_Measures',sep='/'), templBaseName), showWarnings = FALSE)
    name.dir<-paste(graph_path,'Graph_Measures', templBaseName,sep='/')
    
    plot_mvt(name.dir, rs_path)#output of movement parameters
    
  # ----------------------------------
    data.roi<-as.matrix(data.roi)
    data.roi<-data.roi[,regions_selected] 
    n.regions<-length(regions_selected)

    cat('serie_temp ',serie_temp,'n regions ',n.regions,' proc length ',proc.length,'\n')
    n.levels<-4 ## maximum number of wavelet scale to extract
    cor.list <- const.cor.list(data.roi, method = 'modwt', wf = 'la8',
                               n.levels = n.levels, boundary = 'periodic',
                               save.wave = FALSE, export.data = FALSE) # calcul des mat
    
    for(j in 1:n.regions){
      if(sum(data.roi[,j])==0){
        warning(paste('problem region',j,'equal to 0'))
        if(j ==1){
          for(k in 1:n.levels){
            cor.list[[k]][,j]<-c(1,rep(0,length((j+1):n.regions)))
            cor.list[[k]][j,]<-c(1,rep(0,length((j+1):n.regions)))
          }
        }
        else if(j==n.regions){
          for(k in 1:n.levels){
            cor.list[[k]][,j]<-c(rep(0,length(1:(j-1))),1)
            cor.list[[k]][j,]<-c(rep(0,length(1:(j-1))),1)
          }
        }
        else{
          for(k in 1:n.levels){
            cor.list[[k]][,j]<-c(rep(0,length(1:(j-1))),1,rep(0,length((j+1):n.regions)))
            cor.list[[k]][j,]<-c(rep(0,length(1:(j-1))),1,rep(0,length((j+1):n.regions)))
          }
        }				
      }
    }
    print('Extraction Correlation...')
    
    for(i in 1:n.levels){
      
      write.table(cor.list[[i]],paste(name.dir,'/','/','wave.cor.mat_n.levels_',i,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      write.table(cor.list[[i+n.levels]],paste(name.dir,'/','/','wave.cor.mat.lower_n.levels_',i,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      write.table(cor.list[[i+2*n.levels]],paste(name.dir,'/','/','wave.cor.mat.upper_n.levels_',i,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      #mtc: modifications to save correlation plot
      CairoPNG(file=paste(name.dir,'/','/','Corr.mat_n.levels_',i,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.grey.matter.png',sep=''),width = 1000, height = 800,family='times')
      col1 <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "yellow",  "#FF7F00", "red", "#7F0000"))
      corrplot(cor.list[[i]], method="color", col = col1(100))
      dev.off()
    }
    
    print('Graphs...')	
    
    if(graphs){
      ##browser()
      cor.mat<-read.table(paste(name.dir,'/','wave.cor.mat_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.grey.matter.txt',sep=''))
      cat(paste(name.dir,'/','wave.cor.mat_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.grey.matter.txt',sep=''),'\n') 
      
      cor.mat<-as.matrix(cor.mat)
      
      Eglob<-array(0,dim=c(length.nb.edges)) # efficacite global= plus court chemin 1 valeur par nb aretes
      # moyenne sur le graphe de l'efficacite globale
      Eloc<-array(0,dim=c(length.nb.edges))# moyenne sur le graphe de l'efficacite locale = clustering
      tot.Eglob<-array(0,dim=c(n.regions,length.nb.edges)) # sur le total des regions 78X40
      tot.Eloc<-array(0,dim=c(n.regions,length.nb.edges))
      in.degree<-array(0,dim=c(n.regions,length.nb.edges))# degres
      
      tot.nbedges<-array(0,dim=c(length.nb.edges))# nb d'aretes selectionnes en tout car pas exactement le nb demandé au seuil (si 400 peut etre 395, 402 etc)
      tot.sup<-array(0,dim=c(length.nb.edges)) # valeur de correlation utilisee pr mettre le seuil (voir si pas tp faible ou pas tp fort)
      
      tot.BC<-array(0,dim=c(n.regions,length.nb.edges)) 
      SmallWorld<-array(0,dim=c(4,length.nb.edges)) 
      
      ## parametres de graph	    
      Percolation<-array(0,dim=c(length.nb.edges))
      total_distance_connections<-array(0,dim=c(length.nb.edges))
      ratio_long_dist_connections<-array(0,dim=c(length.nb.edges))
      max.nodes.removal<-n.regions
      attack_robustness<-array(0,dim=c(max.nodes.removal,length.nb.edges))
      attack_robustness_nodes<-array(0,dim=c(max.nodes.removal,length.nb.edges))
      
      total_distance_connections<-array(0,dim=c(length.nb.edges))
      ratio_long_dist_connections<-array(0,dim=c(length.nb.edges))
      max.nodes.removal<-n.regions
      attack_robustness<-array(0,dim=c(max.nodes.removal,length.nb.edges))
      attack_robustness_nodes<-array(0,dim=c(max.nodes.removal,length.nb.edges))
      total_distance_connections_nodes<-array(0,dim=c(n.regions,length.nb.edges))
      
      fastgreedy.merges<-vector('list',length.nb.edges) # modularite
      fastgreedy.modularity<-vector('list',length.nb.edges)
      # matrice d'adjacence (meme taille mat corr) => meme VARIANCE pr ttes les paires, calculee a partir du nb de pts ds la serie temporelle
      adj.mat<-const.adj.mat(cor.mat, thresh = 0.05, sup = 0.000000001, proc.length = proc.length,num.levels=num.levels) #on teste si c sup a zero
      
      ################### /!\ /!\ /!\ A mettre en relation avec le nb d'aretes qu'on a choisi d'extraire au debut		
      corr_zero<-sum(adj.mat)/2 # nb de correlation signif != 0 ( ne depend plus du seuil, nb de correlation qu'on a le dt d'extraire pr une mat donnee)
      cat('nb correlation signif diff zero ', c(sum(adj.mat)/2),'\n') # affiche le nb de correlation signif
      
      abs.cor.mat<-abs(cor.mat)
      abs.cor.mat[adj.mat==1]->min_signif
      min(min_signif)
      cat(min(min_signif),'\n') #affiche valeur min des corr parmi le nb significatif diff de 0
      ######################## boucle sur le nb d'aretes
      count<-1
      
      max.num.edges <- factorial(n.regions)/(2*factorial(n.regions-2)) # on calcule avec les combinaisons de 2 parmi le nombre de régions n!/(k!(n-k)!)
      cost<-round(num.nb.edges/max.num.edges*100*(n.regions*(n.regions-1))/2/100) # nb d'aretes en % par rapport au nb total possible tous les 2.5% on commence a 5% car sinon 75 aretes pas suffisant pr minimum spaning tree
      
      coord<-paste(file_coord)
      set1 <- read.table(coord)
      set1 <- set1[regions_selected,]
      #browser()
      
      euclid <- 2 * dist(set1, method = "euclidean") 
      x.euclid <- as.matrix(euclid)
      
      # on ne garde pas les regions deconnectees donc on cree le minimum spanning tree    
      for(i in cost){
        if (i<=corr_zero){ #MTC: Remove this line to allow graphs without significance constrain
          # OK
          #### MST
          diag(abs.cor.mat)<-rep(0,n.regions)
          MST<-mst(sqrt(2*(1-abs.cor.mat)))>0 ### ajacency matrix with MST connextions
          MST_adj_mat<-MST
          MST_adj_mat[lower.tri(MST_adj_mat)]<-0
          adj.mat<-matrix(0,n.regions,n.regions)
          pvalue.cor <- abs(cor.mat[upper.tri(cor.mat)])
          MST_no_edges<-(MST_adj_mat[upper.tri(MST_adj_mat)]==0)
          cor_wo_MST<-pvalue.cor[MST_no_edges]
          
          pvalue.thresh <- sort(cor_wo_MST,decreasing=T)[i-(n.regions-1)] ### only valid for number of edges greater than n.regions
          n.sup<-pvalue.thresh # plus petite valeur de corr pr obtenir i aretes
          test.sign <- (pvalue.cor >= pvalue.thresh)
          l <- 1
          for (k in 2:(n.regions)) {
            for (q in 1:(k - 1)) {
              if ((test.sign[l]) == TRUE) {
                adj.mat[k, q] <- 1
                #cat('no MST edges',k,' ',q,' ',abs.cor.mat[k,q],'\n')
              }
              else{
                if(MST[k,q]==1){ 
                  adj.mat[k, q] <- 1
                  #cat('MST edges',k,' ',q,' ',abs.cor.mat[k,q],'\n')
                }
              }
              l <- l + 1
            }
          }
          
          adj.mat<-adj.mat+t(adj.mat) # matrice d'adjacence avec i aretes premier a 5%=150 aretes
          ##browser()
          
          if(count %in% c(percentage_selected-1, percentage_selected, percentage_selected +1)){
            write.table(adj.mat,paste(name.dir,'/','Adj_mat_n.levels_',num.levels,'_n.regions_',n.regions,'proc_length',proc.length, '.cost_', i,'.txt',sep=''),col.names=F,row.names=F,quote=F)
          }
          
          print(c(i,"*****",n.sup,"*****",sum(adj.mat)/2))
          
          tot.sup[count]<-n.sup
          
          tot.nbedges[count]<-sum(adj.mat)/2
          in.degree[,count]<-rowSums(adj.mat)
          
          tmp<-global.efficiency(adj.mat,weight.mat=matrix(1,n.regions,n.regions))
          Eglob[count]<-tmp$eff
          tot.Eglob[,count]<-tmp$nodal.eff
          
          tmp1<-local.efficiency(adj.mat,weight.mat=matrix(1,n.regions,n.regions))
          
          Eloc[count]<-tmp1$eff
          
          tot.Eloc[,count]<-tmp1$loc.eff
          
          tmp2<-small.world(adj.mat, dat = "reduced", distance = "norm", coord = 0, export.data = FALSE)
          
          Percolation[count]<-tmp2$size.large.connex
          
          total_distance_connections[count]<-sum(adj.mat*x.euclid)
          
          tmp3<-sum(adj.mat*(x.euclid>85))
          
          ratio_long_dist_connections[count]<-100*tmp3/sum(adj.mat)/2
          
          tmp4<-targeted.attack(adj.mat,max.nodes.removal)
          
          attack_robustness[,count]<-tmp4$size.large.connex
          attack_robustness_nodes[,count]<-tmp4$rem.nodes
          
          total_distance_connections_nodes[,count]<-rowSums(adj.mat*x.euclid)
          tmp.graph<-graph.adjacency(adjmatrix=adj.mat,mode='undirected',weighted=NULL,diag=F)
          
          res1<-fastgreedy.community(tmp.graph)
          fastgreedy.merges[[count]]<-res1$merges
          fastgreedy.modularity[[count]]<-res1$modularity
          tot.BC[,count]<-betweenness(tmp.graph, directed=FALSE) #MTC: node betweenness centrality computation
          SmallWorld[,count]<-unlist(small.world(adj.mat)) #MTC: node betweenness centrality computation
          print(count)
          
          if(corr_zero>cost[percentage_selected]){
            
            if(count==percentage_selected){	        
              community<-as.vector(membership(res1))
              dim.edges<-cost[count]
              name.png<-paste('n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,sep='')
              plot_graphs(adj.mat,name.dir,file_coord,regions_selected,dim.edges,name.png)
              
              nodal_metrics<-in.degree[,count]
              if(length(nodal_metrics)!=dim(original.data)[2]){
                nodal_metrics<-c(nodal_metrics,rep(0,(dim(original.data)[2]-length(nodal_metrics))))
              }
              
              name.nii<-paste(name.dir,'/degree_',name.png,sep='')
              plot_nifti(atlas_file,nodal_metrics,name.nii)
              
              nodal_metrics<-tot.Eglob[,count]*100
              if(length(nodal_metrics)!=dim(original.data)[2]){
                nodal_metrics<-c(nodal_metrics,rep(0,(dim(original.data)[2]-length(nodal_metrics))))
              }
              
              name.nii<-paste(name.dir,'/Eglob_',name.png,sep='')
              plot_nifti(atlas_file,nodal_metrics,name.nii)
              
              nodal_metrics<-community
              if(length(nodal_metrics)!=dim(original.data)[2]){
                nodal_metrics<-c(nodal_metrics,rep(0,(dim(original.data)[2]-length(nodal_metrics))))
              }
              
              name.nii<-paste(name.dir,'/modules_',name.png,sep='')
              plot_nifti(atlas_file,nodal_metrics,name.nii)
              
              nodal_metrics<-tot.Eloc[,count]*100
              if(length(nodal_metrics)!=dim(original.data)[2]){
                nodal_metrics<-c(nodal_metrics,rep(0,(dim(original.data)[2]-length(nodal_metrics))))
              }
              
              name.nii<-paste(name.dir,'/Eloc_',name.png,sep='')
              plot_nifti(atlas_file,nodal_metrics,name.nii)
              
              betweenness_centrality<-betweenness(tmp.graph, directed=FALSE) #MTC: node betweenness centrality computation
              param.sw.brain<-small.world(adj.mat,dat="all")          		
            }
          }
          count<-count+1
        }
      }
      
      
      write.table(Eglob,paste(name.dir,'/','Eglob_mean_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(Eloc,paste(name.dir,'/','Eloc_mean_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(tot.sup,paste(name.dir,'/','Thresh_mean_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(tot.nbedges,paste(name.dir,'/','Nb.edges_mean_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(in.degree,paste(name.dir,'/','In.degree_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      
      write.table(tot.Eglob,paste(name.dir,'/','Eglob_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(tot.Eloc,paste(name.dir,'/','Eloc_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      
      write.table(c(corr_zero,min_signif),paste(name.dir,'/','corr_zero_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(Percolation,paste(name.dir,'/','Percolation_n.levels_',num.levels,'_n.regions_', n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(total_distance_connections,paste(name.dir,'/','total_distance_connections_n.levels_', num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(ratio_long_dist_connections,paste(name.dir,'/','ratio_long_dist_connections_n.levels_' ,num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(attack_robustness,paste(name.dir,'/','attack_robustness_n.levels_',num.levels, '_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(attack_robustness_nodes,paste(name.dir,'/','attack_robustness_nodes_n.levels_', num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(total_distance_connections_nodes,paste(name.dir,'/', 'total_distance_connections_nodes_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F)
      
      write.table(tot.BC,paste(name.dir,'/','BC_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter.txt',sep=''),col.names=F,row.names=F,quote=F) #MTC
      if(corr_zero>400){
        write.table(betweenness_centrality,paste(name.dir,'/','BC_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length, '.mst.grey.matter___column4.txt',sep=''),col.names=F,row.names=F,quote=F) #MTC
      }
      save(fastgreedy.merges,file=paste(name.dir,'/','fastgreedy.merges_n.levels_',num.levels, '_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''))
      save(fastgreedy.modularity,file=paste(name.dir,'/','fastgreedy.modularity_n.levels_',num.levels, '_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,'.mst.grey.matter.txt',sep=''))
      
      ### html outputs
      
      print(getwd())
      
      dir.create(file.path(paste(name.dir), paste('HTML_results_n.levels_',num.levels, '_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,sep='')), showWarnings = FALSE)
      
      html.dir<-paste(name.dir,paste('HTML_results_n.levels_',num.levels, '_n.regions_',n.regions,'_start_time_series_',length_time_series[1],'proc_length',proc.length,sep=''),sep='/')
      
      
      HTMLStart(outdir=html.dir, file="average_metrics",extension="html", echo=FALSE, HTMLframe=TRUE)
      target<-paste(html.dir,'/average_metrics_main.html',sep='')
      HTML.title("Average graph metrics with brainwaver", HR=1)
      
      HTML.title("Significant number of edges", HR=3)
      HTMLhr(file=target)
      HTML(paste('Taking into account the amount of possible number of edges and point in the time series, for this dataset the maximum number of possible edges in the graph is equal to ',corr_zero),file=target)
      HTMLbr(x=1,file=target)
      HTMLhr(file=target)
      HTML('see Achard et al. J Neurosci. 2006',file=target)
      HTMLbr(x=1,file=target)
      
      modularity<-rep(0,length.nb.edges)
      for(m in 1:length.nb.edges){
        modularity[m]<-max(fastgreedy.modularity[[m]])
      }
      
      results<-data.frame(NUM_EDGES=tot.nbedges,E_GLOB=Eglob,E_LOC=Eloc,MODULARITY=modularity, CLUSTERING=SmallWorld[2,], MIN_PATH=SmallWorld[3,])
      
      HTML(results,file=target)
      HTMLStop()
      
      if(cost[percentage_selected]<=corr_zero){
        
        HTMLStart(outdir=html.dir, file="10_percent",extension="html", echo=FALSE, HTMLframe=TRUE)
        target<-paste(html.dir,'/10_percent_main.html',sep='')
        HTML.title("Graph metrics for 10% cost with brainwaver", HR=1)
        
        HTML.title("Table including all regions of template", HR=3)
        
        results<-data.frame(DEGREE=in.degree[,percentage_selected],E_GLOB=tot.Eglob[,percentage_selected],E_LOC=tot.Eloc[,percentage_selected],MODULES=community, BC=tot.BC[,percentage_selected], CLUSTERING=param.sw.brain$Cp,MIN_PATH_LENGTH=param.sw.brain$Lp)
        
        
        HTML(results,file=target)
        HTMLStop()
      }
      
    }
    print('The graph computation is finished!!!')
  }
  
  #############################################################################################################################################
  
  #compute_Graph(graph_path,templBaseName,rs_path,file_coord,1, 1, 90,2,seq(100,400,100))
  #compute_Graph <- function(graph_path,templBaseName,rs_path,file_coord,compute_cor, graphs, regions_selected=90, num.levels,num.nb.edges)
  