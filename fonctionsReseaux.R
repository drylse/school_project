require(igraph)


### fonctions utilitaires
getAdj = function(G)
{
  as.matrix(get.adjacency(G))
}

# faire un graphe a partir d'une matrice d'adjacence
plotAdj = function(A)
{
  dir = isSymmetric(A)
  if (dir==TRUE) plot(graph.adjacency(A,mode="undirected"))
  else plot(graph.adjacency(A,mode="directed"))
}


Communaute=function (N ,nedge,cluster, rapport)   
{
  
  #frequency in clusters
  alpha=as.numeric(unclass(table(cluster)))
  nclust=length(alpha)
  totedgeintra=sum(choose(alpha,2))
  totedgeinter=choose(N,2)-totedgeintra
  
  #number of edge inter and intra
  
  if (length(rapport)==1)
  {
    ninter=rbinom(1,nedge,1/(rapport+1))
    nintra= nedge-ninter 
    cptwhile=0
    
    while ((ninter<(nclust-1))|(nintra>totedgeintra))
    {
      ninter=rbinom(1,nedge,1/(rapport+1))
      nintra= nedge-ninter 
      cptwhile=cptwhile+1
      if (cptwhile>50) {break}
    }
  }
  else 
  {
    nintra=rapport[1]
    ninter=rapport[2]
    
  }
  
  if (nintra>totedgeintra | ninter>totedgeinter | ninter+nintra>nedge) 
  {
    print("nombres de liens intra/inter bloques par defaut")
    nintra=totedgeintra
    ninter=nedge-nintra
    
  }
  
  
  #tirage des aretes intra et inter
  lintra=sort(sample(1:totedgeintra,nintra,replace=FALSE))
  linter=sort(sample(1:totedgeinter,ninter,replace=FALSE))
  if (length(linter)<1){linter=0}
  
  A=matrix(0,nrow=N,ncol=N)
  
  cptintra=0
  cptinter=0
  
  for (i in 1:(N-1))
  {
    for (j in (i+1):N)
    {
      if (cluster[i]==cluster[j])
      {
        cptintra=cptintra+1
        if (lintra[1]==cptintra)
        {
          lintra=lintra[-1]
          if (length(lintra)<1){lintra=0}
          A[i,j]<-1->A[j,i]
        }
      }
      else 
      {
        cptinter=cptinter+1
        if (linter[1]==cptinter)
        {
          linter=linter[-1]
          if (length(linter)<1){linter=0}
          A[i,j]<-1->A[j,i]
        } 
      }
    } 
  }
  G=graph.adjacency(A,mode='undirected')
  return(G)
}

## Intervention sur le réseau :

# ferme classe si cas contaminé :

ferme_classes = function(A, t, etatclasses, enfantsdetect, Nenf, nb_classes)
{
  for (enfant in enfantsdetect)
  {
    for (c in 1:nb_classes)
    {
      if (Nenf*(c-1)/nb_classes +1 < enfant & enfant < Nenf*c/nb_classes) #s'il appartient à la c-ième classe
      {
        # alors on note que cette classe est fermée
        etatclasses[1][c]=1
        etatclasses[2][c]=t
      }
    }
  }
  return(A)
}

maj_ferm_cl = function(A,classes,etatclasses,Nenf,nb_classes=4,tpsfermC=3)
{
  for (c in 1:nb_classes)
  {
    if (etatclasses[1][c] == 1 & etatclasses[2][c]>t-tpsfermC) # si la classe est notée fermée et le tps de ferm n'est pas écoulé
    {
      # alors on enlève les relations au sein de cette classe
      A[Nenf*(c-1)/nb_classes +1:Nenf*c/nb_classes, Nenf*(c-1)/nb_classes +1:Nenf*c/nb_classes] = 0
    }  
    if (etatclasses[1][c] == 1 & etatclasses[2][c]<t-tpsfermC) # si le tps de fermeture est écoulé 
    {
      etatclasses[1][c]=0 
      etatclasses[2][c]=0
    }
  }
  return(A)
}



# pour voir l'effet de la ferm. des écoles, on construit une nouvelle matrice d'adjacence, identique mais sans les contacts enfants-enfants

ferme_ecole = function(A, N, Nenf)
{
  ecole_fermee = matrix(0, nrow = N, ncol = N) # il suffit de rappeler la fonction familles sur une matrice nulle
  familles(ecole_fermee, N, Nenf)
  relationsAA(ecole_fermee, N, Nenf)
  return (ecole_fermee) # matrice avec uniquement les contacts familiaux et AA
}
# pas besoin de fct rouvre_ecole : on retournejuste à matrice A

