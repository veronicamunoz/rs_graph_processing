
num.levels <- 3
n.regions <- 100
n.proc <- 500
control.dir <- '/media/veronica/DATAPART2/EmoPark/Graphs/Controls/' 
patient.dir <- '/media/veronica/DATAPART2/EmoPark/Graphs/Patients/'
atlas.name <- 'wAAL3_mod_v1'
num.nb.edges = seq(200,3000,200)


list.control <- dir(control.dir)
list.patient <- dir(patient.dir)
list.patient <- list.patient[list.patient != "PAT15"]

eglob.file.name <- paste0('Eglob_mean_n.levels_',num.levels,'_n.regions_',n.regions,'_start_time_series_1proc_length',n.proc,'.mst.grey.matter.txt')
eloc.file.name <- gsub('Eglob','Eloc',eglob.file.name)

max.num.edges <- factorial(n.regions)/(2*factorial(n.regions-2))
cost <- round(num.nb.edges/max.num.edges,3) 

eglob.control <- data.frame(matrix(NA,nrow = 15, ncol = length(list.control)),row.names = cost)
eglob.patient <- data.frame(matrix(NA,nrow = 15, ncol = length(list.patient)),row.names = cost)

eloc.control <- data.frame(matrix(NA,nrow = 15, ncol = length(list.control)),row.names = cost)
eloc.patient <- data.frame(matrix(NA,nrow = 15, ncol = length(list.patient)),row.names = cost)

colnames(eglob.control) <- colnames(eloc.control)  <- list.control
colnames(eglob.patient) <- colnames(eloc.patient)  <- list.patient


for(subject in list.control){
  eglob.control[subject] <- read.table(file.path(control.dir,subject,'Graph_Measures',atlas.name,eglob.file.name))
  eloc.control[subject] <- read.table(file.path(control.dir,subject,'Graph_Measures',atlas.name,eloc.file.name))
}

for(subject in list.patient){
  eglob.patient[subject] <- read.table(file.path(patient.dir,subject,'Graph_Measures',atlas.name,eglob.file.name))
  eloc.patient[subject] <- read.table(file.path(patient.dir,subject,'Graph_Measures',atlas.name,eloc.file.name))
}



library(ggplot2)
library(reshape)
eglob.control[ "cost" ] <- as.numeric(rownames(eglob.control))
control.molten <- melt( eglob.control, id.vars="cost", value.name="eglob", variable.name="name" )
colnames(control.molten) <- c("cost", "name", "eglob")
eglob.control[eglob.control==0] = NA

eglob.patient[ "cost" ] <- as.numeric(rownames(eglob.patient))
patient.molten <- melt( eglob.patient, id.vars="cost", value.name="eglob", variable.name="name" )
colnames(patient.molten) <- c("cost", "name", "eglob")
eglob.patient[eglob.patient==0] = NA


eglob <- rbind(control.molten,patient.molten)
eglob$type <- c(rep("control", nrow(control.molten)), rep("patient", nrow(patient.molten)))
eglob$type <- as.factor(eglob$type)

eglob$eglob[eglob$eglob==0] = NA

ggplot(eglob,aes(x=reorder(cost,-name), y=eglob, group=name, color = type)) + 
  geom_line() +
  labs(x="Cost",
       y="Eglob") +
  theme_classic()

mean.eglob <- data.frame(cost = rep(cost,3), 
                         group = c(rep("control", length(cost)), rep("patient", length(cost)), rep("diff", length(cost))),
                         eglob = c(apply(X = eglob.control, MARGIN = 1, FUN = mean, na.rm=TRUE), apply(X = eglob.patient, MARGIN = 1, FUN = mean,na.rm=TRUE), (apply(X = eglob.control, MARGIN = 1, FUN = mean,na.rm=TRUE) + apply(X = eglob.patient, MARGIN = 1, FUN = mean,na.rm=TRUE))/2)
)

ggplot(mean.eglob, aes(x=cost, y=eglob, group=group, color=group)) +
  geom_line() +
  scale_x_continuous(breaks=cost) +
  labs(x="Cost",
       y="Eglob") +
  geom_vline(xintercept = cost[3]) +
  geom_vline(xintercept = cost[4]) +
  theme_classic()




mtcars$cyl <- as.factor(mtcars$cyl)
head(mtcars)
library(ggplot2)
# Nuage de points simples
ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point()
# Changer la taille et la forme
ggplot(mtcars, aes(x=wt, y=mpg)) +
  geom_point(size=2, shape=23)