calculSISDeterministe = function(m,e,A,Xini,nT)
{
  N = length(Xini)
  Icour = mean(Xini)
  Scour = 1 - Icour
  # info A nb voisin moyen
  tau = m * sum(A)/N
  
  vI = numeric(nT)
  vS = numeric(nT)
  
  for (t in 1:nT)
  {
    Inew = Icour + Icour * Scour * tau - e * Icour  
    Snew = 1 - Inew
    
    vI[t] = Inew
    vS[t] = Snew      
      
    Icour = Inew
    Scour = Snew
  }
  
  return(list(nI=vI*N,nS=vS*N))
    
}

calculSISProba = function(m,e,A,Xini,nT)
{
  probacour = Xini
  N = length(Xini)
  mproba = matrix(0,nT,N)
  
  for (t in 1:nT)
  {
    probnoncont =  apply(1-A*m*matrix(probacour,N,N,byrow=F),2,prod) 
    probacour = 1 - ( (1-probacour)*probnoncont + e * probacour  )  
    mproba[t,] = probacour
  }
  
  return(mproba)
  
}
