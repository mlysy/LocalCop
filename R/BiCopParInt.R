################################################
# Permisible Intervals for Copula Parameter
################################################

BiCopParInt = function(family) {
  if(family==1 || family==2){ intv  <-  c(-1,1)}
  if(family==3 || family==13){ intv  <-  c(0,100)}
  if(family==4 || family==14 || family== 6){ intv  <-  c(1,100)}
  if(family==5){ intv  <-  c(-100,100)}

  if(family==7 || family==17){intv  <-  c(0,7)}
  if(family==8 || family==9 || family==18 || family==19){intv  <-  c(1,6)}
  if(family==10 || family==20){intv  <-  c(1,8)}

  if(family==23 || family==33){intv  <-  c(-100,0)}
  if(family==24 || family==34  || family==26 || family==36){intv  <-  c(-100,-1)}

  if(family==27 || family==37){intv  <-  c(-7,0)}
  if(family==28 || family==29 || family==38 || family==39){intv  <-  c(-6,-1)}
  if(family==30 || family==40){ intv  <-  c(-8,-1) }

  if(family==104 || family==114 || family==204 || family==214 ){intv  <-  c(1,500)}
  if(family==124 || family==134 || family==224 || family==234){intv  <-  c(-500,-1)}
  return(intv)
}





