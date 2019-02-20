#--- local likelihood for frank copula -----------------------------------------

llfrank <- function(u1, u2, z, wgt, eta) {
  theta <- BiCopEta2Par(family = 5, eta = eta[1] + eta[2]*z)
  ds <- VineCopula::BiCopPDF(u1, u2, family = 5,
                             par = theta$par, par2 = theta$par2)
  sum(wgt * log(ds))
}
