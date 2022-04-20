#install.packages("igraph")
#install.packages("parallel")

library(igraph) # utile pour simuler des réseaux et les représenter
#source("fonctionsReseaux_cours.R") # fonctions R codées  pour simuler et agir sur des réseaux
source("fonctionsReseaux.R")
source("fonctionsSimulations.R") # fonctions R codées pour simuler des modèles SIR
source("fonctionsModDeterministes.R") # fonctions R codées pour les approximations déterministes

# Descripteurs du réseau
N = 200 # effectif total
Nenf = 100 # effectif des enfants
nb_classes = 4
probaAA = 5
pdetect = 0.2
tpsisolement=10 # temps d'isolement
tpsfermC=3 # temps de fermeture des classes
tpsfermE=21 # temps de fermeture de l'école

classes = genere_ecole(Nenf, nb_classes)
AD = genere_pop(N=N, Nenf=Nenf, classes=classes, probaAA=probaAA)
AD

# Paramètres de l’épidémie
m = 0.1# taux de contamination 
e =.28 # taux de guérison
val = e/(m*getvap(AD))
print(val)
nT = 150 # nombre de pas de temps
nb_infec_init = 5
Xini = rep(0,N) # vecteur de N éléments valant 0
# choix aléatoire de 5 indices pour les 5 premiers individus contaminés :
prem_cont = sample(x = 1:N, size = nb_infec_init)
Xini[prem_cont] = 1

# Simulations
nrep = 20


## Fonctions de constitution de la matrice des contacts

# Contacts adulte-adulte :
relationsAA = function(AD, N,  Nenf, probaAA = 5) {
  for (i in seq(Nenf+1:N)) {
    for (j in seq(Nenf+1:N)){
      if (i != j) {
        if (runif(1,0,100) < probaAA) {
          AD[i,j] = 1 # l'adulte a 40% de proba d'etre en contact avec un autre adulte
          AD[j,i] = 1 
        }
      }
    }
  }
  return(AD)
}


# Familles : 2 frères/soeurs + 2 parents, tous en contacts entre eux

familles = function(AD, N, Nenf) { # fonctionne ssi Nenf <= N)
  for (i in seq(1,Nenf,by=2)) { # on parcourt 1 enfant sur 2, pour créer des fratries de 2 enfants
    AD[i,i+1] = 1 # l'enfant i est en contact avec l'enfant i+1 car ils sont de la même famille
    AD[i+1,i] = 1 # et vice versa
    
    # Liens du parent 1 (i+Nenf) avec les enfants
    AD[i+Nenf, i] = 1 # parent 1 - enfant 1
    AD[i, i+Nenf] = 1 # réciproque
    AD[i+Nenf, i+1] = 1 # parent 1 - enfant 2 
    AD[i+1, i+Nenf] = 1 # réciproque 
    
    # Liens du parent 2 (i+1+Nenf) avec les enfants
    AD[i+1+Nenf, i] = 1 # parent 2 - enfant 1
    AD[i, i+1+Nenf] = 1 
    AD[i+1+Nenf, i+1] = 1 # parent 2 - enfant 2 
    AD[i+1, i+1+Nenf] = 1 
    
    AD[i+Nenf, i+1+Nenf] = 1 # lien entre parent 1 et parent 2
    AD[i+1+Nenf, i+Nenf] = 1 
  }
  return(AD)
}


## Fonction génératrice d'une population modulable (choix du nombre de personnes, d'enfants/adultes, des probas de contact entre chaque groupe)

# Constitution de l'école : 4 classes d'enfants
genere_ecole = function(Nenf, nb_classes=4) 
{
  cl = rep(1:nb_classes,each=Nenf/nb_classes) # marche que pour nb_classes et Nenf multiples
  G = Communaute(N=Nenf,nedge=2*N,cluster=cl,rapport=20)
  classes = getAdj(G) # matrice d'adjacence du réseau d'enfants
  plot(G)
  return(classes)
}

# Constitution de la matrice des contacts globale
genere_pop = function(N, Nenf, classes, probaAA=5) 
{
  AD = matrix(data = 0, nrow = N, ncol = N) # matrice nulle de taille N*N
  AD[1:Nenf,1:Nenf] = classes
  AD = familles(AD, N=N, Nenf=Nenf)
  AD = relationsAA(AD, N=N, Nenf=Nenf, probaAA=probaAA)
  return(AD)
}


## MODELES ##

# Modèle SIR classique 
sim_SIR = SimuSIR(m,e,AD,Xini,nT)

# Modèle SIR test/isolement 
sim_test = SimuSIR_isol(m,e,AD,Xini,nT,pdetect,tpsisolement=10)

# Modèle SIR fermeture des classes
sim_classe = SimuSIR_classe(m,e,AD,Xini,nT,pdetect,Nenf,nb_classes,tpsfermC=3)

# Modèle SIR fermeture de l'école
sim_ecole = SimuSIR_ecole(m,e,AD,Xini,nT,pdetect)#,tpsfermE=21,tpsisolement=10)

# Plot du modèle SIR
plot(rowSums(sim_SIR==0), main='S(t), I(t) et R(t) en modèle SIR', type ='l', ylim=c(0,200),xlab="temps",ylab="population") # ajout des nrep S(t) selon les simulations
lines(rowSums(sim_SIR==1),col=2)
lines(rowSums(sim_SIR==2),col=3)
legend("right",lty=1,col=c(1,2,3),legend=c("S","I","R")) #,"Contacts coupés"))

# Plot du modèle test/isolement
plot(rowSums(sim_test==0),main='S(t), I(t) et R(t) en modèle SIR avec tests et isolement', type ='l', ylim=c(0,200),xlab="temps",ylab="population")
lines(rowSums(sim_test==1),col=2)
lines(rowSums(sim_test==2),col=3)
legend("right",lty=1,col=c(1,2,3),legend=c("S {isolement}","I {isolement}","R {isolement}"))

# Plot du modèle fermeture de classe
plot(rowSums(sim_classe==0), main='S(t), I(t) et R(t) en modèle SIR avec fermeture de classes', type ='l', ylim=c(0,200),xlab="temps",ylab="population") 
lines(rowSums(sim_classe==1),col=2)
lines(rowSums(sim_classe==2),col=3)
legend("right",lty=1,col=c(1,2,3),legend=c("S {ferm. classe}","I {ferm. classe}","R {ferm. classe}"))

# Plot du modèle fermeture de l'école 
plot(rowSums(sim_ecole==0), main="S(t), I(t) et R(t) en modèle SIR avec fermeture d'école", type ='l', ylim=c(0,200),xlab="temps",ylab="population") 
lines(rowSums(sim_ecole==1),col=2)
lines(rowSums(sim_ecole==2),col=3)
legend("right",lty=1,col=c(1,2,3),legend=c("S {ferm. école}","I {ferm. école}","R {ferm. école}"))

# Comparaison du nombre d'infectés/de résistants entre les 2 modèles en fonction de différents paramètres m et e

comparaison_I_selon_m = function(AD,Xini,nT,pdetect=0.2,tpsisolement,Nenf,nb_classes = 4,tpsfermC=3)
{
  for (e in seq(0.05,0.5,by=0.05)) 
  { # de 1% à 10%, tous les 5% (=1,6,11% juste)
    nT = 90
    print(c('e',e))
    couleurs = c(1,'yellow2','orange3',2,6) # liste des couleurs de chaque valeur de m
    i = 1 # compteur des couleurs
    plot(NULL,xlim=c(0,nT),main="Comparaison de I(t) pour plusieurs m, avec e fixé",ylim=c(0,N),xlab="temps",ylab="Infectés")
    for (m in c(0.01,0.05,0.1,0.15))
    { # de 5% à 50%, tous les 5%
      sim_SIR = SimuSIR(m,e,AD,Xini,nT)
      sim_test = SimuSIR_isol(m,e,A,Xini,nT,pdetect,tpsisolement=10)
      sim_classe = SimuSIR_classe(m,e,A,Xini,nT,pdetect=0.2,Nenf,nb_classes = 4,tpsfermC=3)
      sim_ecole = SimuSIR_ecole(m,e,A,Xini,nT,pdetect=,tpsfermE=21,tpsisolement=10)
      
      lines(rowSums(sim_SIR==1),col=couleurs[i], type = "l", lty = 1) # SIR en trait plein
      lines(rowSums(sim_test==1),col=couleurs[i], type = "l", lty = 2) # avec tests/isolement en segments longs
      lines(rowSums(sim_classe==1),col=couleurs[i], type = "l", lty = 3) # avec fermeture de classe en pointillés
      lines(rowSums(sim_ecole==1),col=couleurs[i], type = "l", lty = 5) # avec fermeture d'école en points
      
      legend("right", lty=1, title="Valeurs de m (en %)",
             col=couleurs,
             legend=c(0.01,0.05,0.1,0.15)*100)
      i <- i+1
    }
  }
}
comIm = comparaison_I_selon_m(AD,Xini,nT) #on constate plus faible influence des interventions avec augmentation de m

comparaison_I_selon_e = function(AD,Xini,nT)
{
  for (m in c(0.05,0.1,0.15,0.2,0.3,0.4)) {
    nT=90
    print(c('m',m))
    couleurs = c(1,'yellow2','orange3',2,6)
    i = 1 
    plot(NULL,xlim=c(0,nT),main="Comparaison de I(t) pour plusieurs e, avec m fixé",ylim=c(0,N),xlab="temps",ylab="Infectés")
    for (e in c(0.05,0.1,0.15,0.2))
    { # de 5% à 50%, tous les 5%
      sim_SIR = SimuSIR(m,e,AD,Xini,nT)
      sim_test = SimuSIR_corrige(m,e,AD,Xini,nT,pdetect)
      sim_classe = SimuSIR_action(m,e,AD,Xini,nT,pdetect=0.2)
      sim_ecole = SimuSIR_ecole(m,e,AD,Xini,nT,pdetect=0.2)
      
      lines(rowSums(sim_SIR==1),col=couleurs[i], type = "l", lty = 1)
      lines(rowSums(sim_test==1),col=couleurs[i], type = "l", lty = 3)
      lines(rowSums(sim_classe==1),col=couleurs[i], type = "l", lty = 3) 
      lines(rowSums(sim_ecole==1),col=couleurs[i], type = "l", lty = 5) 
      
      legend("right", lty=1, title="Valeurs de e (en %)",
             col=couleurs,
             legend=c(0.05,0.1,0.15,0.2)*100)
      i <- i+1
    }
  }
}
comIe = comparaison_I_selon_e(AD,Xini,nT)


comparaison_I = function(AD,Xini,nT) {
  
  for (m in seq(0.01,0.1,by=0.05)) {
    for (e in seq(0.05,0.5,by=0.1)) {
      # print(c('m',m))
      # print(c('e',e))
      plot(rowSums(sim_SIR==1), type='l', ylim = c(0,200),xlab="temps",ylab="Infectés")
      lines(rowSums(sim_test==1),col='yellow2')
      lines(rowSums(sim_classe==1),col='orange3')
      lines(rowSums(sim_ecole==1),col=2)
      
      legend("right", lty=1,
             col=c(1,'yellow2','orange3',2),
             legend=c("I", "I {isolement}", "I {ferm. classe}", "I {ferm. école}"))
    }
  }
}

comI = comparaison_I(AD,Xini,nT)

comparaison_R = function(AD,Xini,nT) {
  for (m in seq(0.01,0.1,by=0.05)) {
    for (e in seq(0.05,0.5,by=0.05)) {
      plot(rowSums(sim_SIR==2), type='l', col=1, ylim = c(0,200),xlab="temps",ylab="immunisés contre Réinfection")
      lines(rowSums(sim_test==2),col='yellow2')
      lines(rowSums(sim_classe==2),col='orange3')
      lines(rowSums(sim_ecole==2),col=2)
      
      legend("right", lty=1,
             col=c(1,'yellow2','orange3',2),
             legend=c("R", "R {isolement}", "R {ferm. classe}", "R {ferm. école}"))
    }
  }
}

comR = comparaison_R(AD,Xini,nT)

# Comparaison du nombre d'infectés dans les 3 situations 
plot(NULL,xlim=c(0,nT),main="Comparaison de I(t) avec et sans tests / 
     fermeture de classes / d'école",ylim=c(0,N),xlab="temps",ylab="Infectés")
for (i in 1:10) {
  lines(rowSums(sim_SIR==1),col=6)
  lines(rowSums(sim_test==1),col=4)
  lines(rowSums(sim_classe==1),col='yellow2')
  lines(rowSums(sim_ecole==1),col='orange3')
}
legend("right",lty=1,col=c(6,4,'yellow2', 'orange3'),legend=c("I", "I {isolement}", "I  {ferm. classes}","I {ferm. école}"))

# Plot de comparation des modèles 

comp_modeles_graph <- function(){
  couleurs = c('hotpink',4,5,6)
  for (m in (1:5)){
    plot(NULL,xlim=c(0,nT),main="S(t), R(t) et I(t) pour les trois types d'interventions et le modèle SIR avec m=%",
         ylim=c(0,N),xlab="temps",ylab="population")
    for (i in 1:5){
      sim_sir = SimuSIR(m/100,e,AD,Xini,nT)
      sim_test = SimuSIR_isol(m/100,e,AD,Xini,nT,0.2)
      sim_classe = SimuSIR_classe(m/100,e,AD,Xini,nT,pdetect=0.2)
      sim_ecole = SimuSIR_ecole(m/100,e,AD,Xini,nT,0.2)
      lines(rowSums(sim_sir==1),col=couleurs[1])
      lines(rowSums(sim_test==1),col=couleurs[2])
      lines(rowSums(sim_ecole==1),col=couleurs[3])
      lines(rowSums(sim_classe==1),col=couleurs[4])
    }
  }
  return(tps)
}

legend("right",lty=1,col=c('hotpink',4,5,6),
       title = ("Modèle"),
       legend=c("SIR","Test/isolements","Classe","Ecole"),box.lty=0)

comparaison = comp_modeles()

comp_modeles_inf <- function(){
  stock = data.frame('obs' = 1:104, 'moy'=0,'tps'=0) # 104 = nombre total d'observations
  for (m in (1:5)){
    tps = matrix(0,nrow=20,ncol=150) 
    for (i in (1:5)){
      sim_sir = SimuSIR(m/100,e,AD,Xini,nT)
      sim_test = SimuSIR_isol(m/100,e,AD,Xini,nT,0.2)
      sim_classe = SimuSIR_classe(m/100,e,AD,Xini,nT,pdetect=0.2)
      sim_ecole = SimuSIR_ecole(m/100,e,AD,Xini,nT,0.2)
      
      stock$tps[m+4*(i-1)+(m-1)*20] = temps_moit(rowSums(sim_sir==1)) # calcul du temps d'atteinte de 50% de la pop pour une valeur de m et une simulation i
      stock$tps[m+1+4*(i-1)+(m-1)*20] = temps_moit(rowSums(sim_test==1))
      stock$tps[m+2+4*(i-1)+(m-1)*20] = temps_moit(rowSums(sim_classe==1))
      stock$tps[m+3+4*(i-1)+(m-1)*20] = temps_moit(rowSums(sim_ecole==1))
      
      stock$moy[m+4*(i-1)+(m-1)*20] = mean(rowSums(sim_sir==1)) # Calcul de la moyenne d'infectés
      stock$moy[m+1+4*(i-1)+(m-1)*20] = mean(rowSums(sim_test==1))
      stock$moy[m+2+4*(i-1)+(m-1)*20] = mean(rowSums(sim_classe==1))
      stock$moy[m+3+4*(i-1)+(m-1)*20] = mean(rowSums(sim_ecole==1))
    }
  }
  return(stock)
}

res = comp_modeles_inf()

# Fonction permettant de calculer le temps de contamination de 50% de la population

temps_moit <- function(vec){
  t = 1
  i = 0
  while (i < N/2 & t < 150) {
    i = i + vec[t]
    t = t+1
  }
  if (t<150){
    return(t)
  }
  else{
    return(NA)
  }
}

## REPETITIONS de simulations
# Plots du modèle SIR

infecSIR = data.frame(obs=1:nrep, S=0, I=0, R=0)
#plot(NULL,xlim=c(0,nT),ylim=c(0,N),main='S(t), I(t) et R(t)',xlab="temps",ylab="population")
for (i in 1:nrep){
  sim_SIR = SimuSIR(m,e,AD,Xini,nT)
  infecSIR$S[i] = (rowSums(sim_SIR==0))[nT]
  infecSIR$I[i] = sum(rowSums(sim_SIR==1))/nT # moyenne de personnes infectées par jour 
  infecSIR$R[i] = (rowSums(sim_SIR==2))[nT]
}
#legend("right",lty=1,col=c(1,2,3),legend=c("S", "I", "R"))
sains_SIR = sum(infecSIR$S)/(nrep)
moy_SIR = sum(infecSIR$I)/(nrep)
immu_SIR = sum(infecSIR$R)/(nrep)

# Plots du modèle test/isolement
infecI = data.frame(obs=1:nrep, S=0, I=0, R=0)
#plot(NULL,xlim=c(0,nT),ylim=c(0,N),main='S(t), I(t) et R(t) avec tests et isolement', xlab="temps",ylab="population")
for (i in 1:nrep){
  sim_test = SimuSIR_isol(m,e,AD,Xini,nT,pdetect)
  # lines(rowSums(sim_test==0),type="l",ylim=c(0,200))
  # lines(rowSums(sim_test==1),col=2)
  # lines(rowSums(sim_test==2),col=3)
  infecI$S[i] = (rowSums(sim_test==0))[nT]
  infecI$I[i] = sum(rowSums(sim_test==1))/nT
  infecI$R[i] = (rowSums(sim_test==2))[nT]
}
#legend("right",lty=1,col=c(1,2,3),legend=c("S {isolement}", "I {isolement}", "R {isolement}"))

sains_I = sum(infecI$S)/nrep
moy_I = sum(infecI$I)/nrep
immu_I = sum(infecI$R)/nrep

# Plots du modèle fermeture de classe
infecC = data.frame(obs=1:nrep, S=0, I=0, R=0)
#plot(NULL,xlim=c(0,nT),ylim=c(0,N),main='S(t), I(t) et R(t) avec fermeture de classes', xlab="temps",ylab="population")
for (i in 1:nrep){
  sim_classe = SimuSIR_classe(m,e,AD,Xini,nT,pdetect)
  # lines(rowSums(sim_classe==0),type="l",ylim=c(0,200))
  # lines(rowSums(sim_classe==1),col=2)
  # lines(rowSums(sim_classe==2),col=3)
  infecC$S[i] = (rowSums(sim_classe==0))[nT]
  infecC$I[i] = sum(rowSums(sim_classe==1))/nT
  infecC$R[i] = (rowSums(sim_classe==2))[nT]
}
#legend("right",lty=1,col=c(1,2,3),legend=c("S {ferm. classes}", "I {ferm. classes}", "R {ferm. classes}"))

immu_C = sum(infecC$R)/nrep

# Plot du modèle fermeture de l'école 
infecE = data.frame(obs=1:nrep, S=0, I=0, R=0)
#plot(NULL,xlim=c(0,nT),ylim=c(0,N),main="S(t), I(t) et R(t) avec fermeture d'école",xlab="temps",ylab="population")
for (i in 1:20){
  sim_ecole = SimuSIR_ecole(m,e,AD,Xini,nT,pdetect=0.2)
  # lines(rowSums(sim_ecole==0),type="l",ylim=c(0,200))
  # lines(rowSums(sim_ecole==1),col=2)
  # lines(rowSums(sim_ecole==2),col=3)
  infecE$S[i] = (rowSums(sim_ecole==0))[nT]
  infecE$I[i] = sum(rowSums(sim_ecole==1))/nT
  infecE$R[i] = (rowSums(sim_ecole==2))[nT]
}
#legend("right",lty=1,col=c(1,2,3),legend=c("S {ferm. école}", "I {ferm. école}", "R {ferm. école}"))

immu_E = sum(infecE$R)/nrep

# Lien probaAA et m/e :
mt=0.05
et=0.1
i=1
couleurs = c(1,'yellow2','orange3')
plot(NULL,xlim=c(0,nT),main="Comparaison de R(t) avec et sans fermeture de classe,
     avec m=0.05, e=0.1, pour plusieurs probaAA",ylim=c(0,N),xlab="temps",ylab="Immunisés")
for (pA in c(5,20,40))
{ 
  probaAA=pA
  sim_SIR = SimuSIR(mt,et,AD,Xini,nT)
  sim_classe = SimuSIR_classe(mt,et,AD,Xini,nT,pdetect=0.2,Nenf,nb_classes = 4,tpsfermC=3)
  #sim_test = SimuSIR_isol(mt,et,AD,Xini,nT,pdetect,tpsisolement=10)
  #sim_ecole = SimuSIR_ecole(mt,et,AD,Xini,nT,pdetect,tpsfermE=21,tpsisolement=10)
  
  lines(rowSums(sim_SIR==2),col=couleurs[i], type = "l", lty = 1)
  #lines(rowSums(sim_test==2),col=couleurs[i], type = "l", lty = 3)
  lines(rowSums(sim_classe==2),col=couleurs[i], type = "l", lty = 3) 
  #lines(rowSums(sim_ecole==2),col=couleurs[i], type = "l", lty = 3) 
  
  legend("topright", lty=1, title="Valeurs de probaAA",
         col=couleurs,
         legend=c(5,20,40))
  
  i <- i+1
}



## Détermination des paramètres m et e les plus intéressants
determ_e_m = function(A)#, Xini, nT, e, m, Rlim, contactlim)
{
  couleurs = 1:length(c(0.01,0.05,0.1,0.15))
  #e_m_ok = c() # initialisation liste de couples (e, m) qui fonctionnent selon nos critères
  for (e in seq(0.01,0.5,by=0.01))#(liste d’intervalles 0.005 entre 0.001 et 0.25))
  {
    i = 1 
    plot(NULL,xlim=c(0,nT),main=" ",ylim=c(0,N),xlab="temps",ylab="population")
    
    for (m in seq(0.01,0.1,by=0.01))#(liste d’intervalles 0.005 entre 0.001 et 0.25))
    {
      val = e/(m*getvap(A))
      #print(c(m,e))
      if (val > 0.9 & val < 1.1)
      {
        #SimuSIR(m,e,A,Xini,nT)
        print(c(m,e,val))
      }
      i = i+1
    }
  }
}
d = determ_e_m(AD) 
#renvoie par ex : couples tq m = 1% et e= 0.235 0.24 0.245 0.25


# renvoie la proba que l'épidémie ait disparu et le I(T) moyen, R(t) moyen
repetSimuSIR(m,e,AD,Xini,nT,nrep=20) 


