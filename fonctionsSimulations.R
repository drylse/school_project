#Script : fonctionsSimulations.R
#Regroupe les fonctions de simulation en modèle SIR sans et avec interventions dépendantes de l'état des individus.

## Paramètres :

#m taux de migration/colonisation : contaminations
#e taux d'extinction : guérisons
#A matrice d'adjacence du graphe de connection
#Xini liste de l'état initial de chaque individu
#nT nombre de pas de temps
#pdetect probabilité de détection d'une contamination
#tpsisolement nombre de pas durant lesquels un individu détecté est isolé
#tfermC nombre de pas durant lesquels une classe est fermée quand un cas est détecté 
#tfermE nombre de pas durant lesquels l'école est fermée


# Modèle SIR sans intervention
SimuSIR=function(m,e,A,Xini,nT)
  #fonction du modele retournant les noeuds infectes au cours du temps : S=0 I=1 R=2
{
  N=length(Xini)# nombre d'individus
  Xprec=Xini# stocke les individus infectés au pas de temps précédent
  # sortie une ligne = une generation, une colonne = un individu
  statut=matrix(0,nrow=nT,ncol=N) # matrice des statuts de chaque individu à chaque pas de temps
  # 0=Susceptible ; 1=Infecté ; 2=immunisé contre Réinfection
  for (t in 1:nT)
  {
    Ared = A
    Xnew = Xprec # stockage de l'état actuel des individus qui va changer
    # cure
    gueri = rbinom(N,1,e)
    Xnew[which(gueri==1 & Xprec==1)] = 2 # ceux qui ont guÃ©ri et qui Ã©taient malades au temps précédent sont immunisés
    
    # infection
    infection = Ared*matrix(rbinom(N*N,1,m),
                            nrow=N)*matrix(Xnew==1,nrow=N,ncol=N)
    # infection (taille N*N) : susceptibilité pour chaque individu d'infecter chacun, sous reserve m (taux de migration) et selon quels individus sont contaminés
    vect= Xnew 
    subsel=which((vect==0)&colSums(infection)>0) # stock des individus S qui sont susceptibles d'être contaminés
    Xnew[subsel] = 1 # on modifie leur état dans Xnew
    statut[t,] = Xnew # on entre le nouveau statut de la population à t dans la matrice statut
    Xprec = Xnew
    
  }
  return(statut)
}

#Modèle SIR avec tests & isolement
SimuSIR_isol=function(m,e,A,Xini,nT,pdetect,tpsisolement=10)
  #fonction du modele retournant les noeuds infectes au cours du temps (S=0 I=1 R=2), avec une intervention de détection et isolement
{
  N=length(Xini)
  Xprec=Xini
  statut=matrix(0,nrow=nT,ncol=N)
  detection = data.frame('identifiant' = 1:N, 
                         'infectdetect' = 0, 
                         'tempsdet' = 0) #1 = identifiant, 2 = infect_detect, 3 = temps_detect
  # avec infectdetect = 0 par défaut, 1 si I puis 2 apres tpsisolement 
  for (t in 1:nT)
  {
    Ared = A # on reprend la matrice initiale à chaque tour de boucle
    Xnew = Xprec
    
    # Réinsertion des infectés détectés isolés 10j dans la matrice 
    for (i in which(detection$infectdetect == 1)) { # pour chaque individu ayant été détecté
      if (t - detection$tempsdet[i] >= tpsisolement) { # s'il a été isolé pendant au moins 10 jours
        detection$infectdetect[i] = 2 # on le considère guéri et apte à retourner à l'école
      }
    }
    
    # réisolement de ceux détectés à isoler
    isol = which(detection$infectdetect == 1 & t-detection$tempsdet < tpsisolement)
    Ared[isol,] = 0
    Ared[,isol] = 0 # les individus ayant été détectés comme malades sont isolés
    
    Ared = familles(Ared, N, Nenf) # on maintient le contact dans la famille
    
    # cure
    gueri = rbinom(N,1,e)
    Xnew[which(gueri==1 & Xprec==1)] = 2 # ceux qui ont guéri et qui étaient malades au temps précédent sont immunisés
    
    # infection
    infection = Ared*matrix(rbinom(N*N,1,m),
                            nrow=N)*matrix(Xnew==1,nrow=N,ncol=N)
    vect= Xnew 
    subsel=which((vect==0)&colSums(infection)>0) 
    Xnew[subsel] = 1
    statut[t,] = Xnew 
    
    # détection des infectés : ils seront isolés à partir du temps suivant, avant la prochaine infection
    nvxdetect =  which(Xprec == 1 & rbinom(N,1,pdetect)==1 & detection$infectdetect == 0) 
    detection$infectdetect[nvxdetect] = 1
    detection$tempsdet[nvxdetect] = t 
    
    Xprec = Xnew
    
  }
  return(statut)
}


# Modèle SIR avec fermeture des classes
SimuSIR_classe=function(m,e,A,Xini,nT,pdetect=0.2,Nenf,nb_classes = 4,tpsfermC=3)
{
  N=length(Xini)
  Xprec=Xini
  statut=matrix(0,nrow=nT,ncol=N)
  detection = data.frame('identifiant' = 1:N, 
                         'infectdetect' = 0, 
                         'tempsdet' = 0)
  etatclasses = data.frame('etat'=rep(0,nb_classes),
                           'tferm'=rep(0,nb_classes)) 
  # $etat : 0=classe ouverte, 1=classe fermée ; $tferm: pas de temps où la classe ferme
  contacts_coupes=numeric(nT) #  initialisation de la liste des contacts coupés à chaque t
  Nenf  = 100
  tpsfermC = 3
  tpsfermE = 10
  for (t in 1:nT)
  {
    Ared = A
    
    ## REIMPLEMENTATION ETAT CLASSES
    # ouverture des classes après tpsfermC
    for (c in 1:nb_classes)
    {
      if (etatclasses$tferm[c]<t-tpsfermC) # si tpsfermC est écoulé
      {
        etatclasses$etat[c] = 0 # on note la classe comme ouverte
        etatclasses$tferm[c] = 0 # on réinitialise tferm
      }
    }
    
    # refermeture des classes fermées restantes
    for (c in 1:nb_classes)
    {
      if (etatclasses$etat[c] == 1) # si la classe est fermée
      {
        Ared[Nenf*(c-1)/nb_classes +1:Nenf*c/nb_classes, Nenf*(c-1)/nb_classes +1:Nenf*c/nb_classes] = 0 # on isole la classe
      }
    }
    
    # réisolement de ceux à isoler de manière générale
    isol = which(detection$infectdetect == 1 & t-detection$tempsdet < tpsisolement)
    Ared[isol,] = 0
    Ared[,isol] = 0 # les individus ayant été détectés comme malades sont isolés
    
    Ared = familles(Ared, N, Nenf) # on maintient le contact dans la famille
    
    ## CHANGEMENTS D'ETAT DES INDIVIDUS
    Xnew = Xprec 
    
    # CURE
    gueri = rbinom(N,1,e)
    Xnew[which(gueri==1 & Xprec==1)] = 2
    # si le guéri etait détecté, et donc isolé, alors on retire les guéris des détectés -> donc enfants guéris sont plus isolés qu'au tour d'après /!\
    gueridetec = which(detection$infectdetect[gueris]==1) # liste des guéris qui étaient détectés
    detection$infectdetect[gueridetec] = 2
    
    # INFECTION
    infection = Ared*matrix(rbinom(N*N,1,m),
                            nrow=N)*matrix(Xnew==1,nrow=N,ncol=N)
    vect = Xnew 
    subsel = which((vect==0)&colSums(infection)>0)
    Xnew[subsel] = 1
    statut[t,] = Xnew
    Xprec = Xnew
    
    ### ACTION sur le reseau
    
    # détection des infectés
    nvxdetect =  which(statut[t,] == 1 & rbinom(N,1,pdetect)==1 & detection$infectdetect == 0) # ind malades et détectés à t (et pas avant)
    detection$infectdetect[nvxdetect] = 1
    detection$tempsdet[nvxdetect] = t
    
    # fermeture de la classe des enfants nouvellement détectés
    nvxenfantsdet = which(detection$tempsdet == t & detection$identifiant<=Nenf)
    cl = rep(1:nb_classes,each=Nenf/nb_classes) # Nenf-liste de la classe de chaque enfant
    cl_nvx_enf_det = cl[nvxenfantsdet] # liste de la classe de chacun des nouveaux enfants détectés
    for (c in cl_nvx_enf_det)
    {
      if (etatclasses$etat[c] == 0) # si la classe i n'est pas déjà fermée
      {
        etatclasses$etat[c] = 1 # on note la classe comme fermée
        etatclasses$tferm[c] = t
      }
    }
    #print(etatclasses)
    
    contacts_coupes[t] <- sum(A-Ared)/2 # stock de la diff du nb de contacts (initial - à t)
    #print (contacts_coupes)
    
  }
  return(statut)#(list(statut=statut,contacts_coupes=contacts_coupes))
}


SimuSIR_ecole=function(m,e,A,Xini,nT,pdetect=0.2,tpsfermE=21,tpsisolement=10) 
  #fonction du modele retournant les noeuds infectes au cours des generations (S=0 I=1 R=2)
{
  N=length(Xini)
  Xprec=Xini
  statut=matrix(0,nrow=nT,ncol=N) 
  detection = data.frame('identifiant' = 1:N, 
                         'infectdetect' = 0, 
                         'tempsdet' = 0)
  statut_ecole = data.frame('detect' = 0, 'tpsdet' = 0)
  for (t in 1:nT)
  {
    Ared = A
    Xnew = Xprec
    
    # réisolement de ceux à isoler
    isol = which(detection$infectdetect == 1 & t-detection$tempsdet < tpsisolement)
    Ared[isol,] = 0
    Ared[,isol] = 0 # les individus ayant été détectés comme malades sont isolés
    
    # refermeture de l'école si elle est toujours fermée
    if (statut_ecole$detect == 1 & t - statut_ecole$tpsdet < tpsfermE) { # si l'école est fermée depuis moins de tpsfermE
      Ared[1:Nenf,1:Nenf] = 0
    }
    
    Ared = familles(Ared,N,Nenf)
    
    # cure
    gueri = rbinom(N,1,e)
    Xnew[which(gueri==1 & Xprec==1)] = 2
    
    # infection
    infection = Ared*matrix(rbinom(N*N,1,m),
                            nrow=N)*matrix(Xnew==1,nrow=N,ncol=N)
    vect= Xnew 
    subsel=which((vect==0)&colSums(infection)>0) 
    Xnew[subsel] = 1 
    statut[t,] = Xnew 
    
    # détection des infectés
    det = rbinom(N,1,pdetect)
    nvxdetect =  which(statut[t,] == 1 & det==1 & detection$infectdetect == 0) # ind malades et détectés à t (et pas avant)
    detection$tempsdet[nvxdetect] = t
    detection$infectdetect[nvxdetect] = 1
    # si on détecte de nouveaux enfants 
    inf_enf = which(statut[t,] == 1 & det==1 & detection$infectdetect==0 & detection$identifiant<Nenf) # enfants infectés détectés
    if (length(isol) >= 5 & statut_ecole$tpsdet==0) {  # si l'école est ouverte et compte au moins 5 cas
      statut_ecole$detect = 1
      statut_ecole$tpsdet = t
    }
    
    Xprec = Xnew
    
  }
  return(statut)
}


# Valeurs obtenues sur un type d'épidémie :
Rlim = 50 # nb max d'immunisés souhaités sur la pop (quart)
# détermination de probaext, moyI et moyR sur nrep simulations: 
repetSimuSIR = function(m,e,A,Xini,nT,nrep=100)
{
  RES = sapply(1:nrep, function(i){
    sim = SimuSIR_corrige(m,e,A,Xini,nT,pdetect)
    return(c(sum(sim[nT,]==1),sum(sim[nT,]==2))) 
  })
  vec = c(mean(RES[1,]==0),mean(RES[1,]),mean(RES[2,]))
  names(vec) = c("probaext","moyI","moyR")
  # if (mean(RES[2,])>Rlim) {print ("Trop d'immunisés)}
  #else {print("Suffisamment peu d'immunisés")}
  return(vec)
}


