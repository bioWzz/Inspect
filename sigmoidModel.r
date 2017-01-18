sigmoidModel <- function(x, par) 
{
  # h0= par[1]; h1=par[2]; h2=par[3]; t1=par[4]; t2=par[5]; b=par[6]
  par[1]+(par[2]-par[1])*(1/(1+exp((-1*par[4])*(x-par[3]))))
}