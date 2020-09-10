############################################################################
########################## Running all models  #############################
############################################################################

# In this script we are going to compute model 0-3. While in  Model 0 and 
# Model 1 non-exposure effects are considered, in model 2 and model 3 
# the exposure effect is included in three different ways (linear effect, 
# Random walk and SPDE of 1 dimension). We consider four confuding factors
# in all the four models.

# NOTE 1: The following models include the four confounding factors 
# considered in this study. In case, you want to compute the same models
# without considering this confounding factors, just remove it from
# the definition of the formula. For examples, removing lines 43-46, 
# model 0 will be computed without these confounding factors).

#Loading necessary packages
library(sp)
library(INLA)
library(maptools)
library(spatstat)
library(rgeos)
library(RColorBrewer)
library(ggmap)
library(rgdal)


########################################################
# Model 0: No covariates, no specific spatial effects
########################################################

#Load stack data and others
load("Data/Stack-LinRW-models.RData")


t1 <-Sys.time()
E <- inla.stack.data(join.stack)$e

#formula
form <- as.formula( 
  paste0("y ~ ", 
    '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
    f(spatial.field1, model = spde) + 
    f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
    f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
    f(cont.copy4, copy = "spatial.field1", fixed = TRUE) ' ,
    " + ",
    paste(paste0("PARO2059.", 2:4), collapse = " + " ),
    " + ",
    paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
    " + ",
    paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
    " + ",
    paste(paste0("ESTPRE.", 2:4), collapse = " + " )
  )
)

#Estimation
Sys.time()
pp.res <- inla(formula = form, 
  family = rep('poisson', 4), data = inla.stack.data(join.stack),
  control.predictor = list(A = inla.stack.A(join.stack), link = 1,
    compute = TRUE),
  E = inla.stack.data(join.stack)$e, 
  control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE))
Sys.time()

summary(pp.res)
t2 <-Sys.time()

ti <- t2-t1; 
print(
  paste0("Model 0 has used ", ti, " cpu time. And have finished at ",
    Sys.time()))

rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, controls_ppp, 
   kid_ppp, lung_ppp, sto_ppp)


save.image("Results/Res-M0-ConFact.RData")


rm(list = ls())

########################################################

########################################################
#Model 1: No covariates, with specific spatial effects
########################################################

#Load stack data and others
load("Data/Stack-LinRW-models.RData")

t1 <-Sys.time()
E <- inla.stack.data(join.stack)$e

#formula
form <- as.formula( 
  paste0("y ~ ", 
    '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
    f(spatial.field1, model = spde) + 
    f(spatial.field2, model = spde) + 
    f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
    f(spatial.field3, model = spde) + 
    f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
    f(spatial.field4, model = spde) + 
    f(cont.copy4, copy = "spatial.field1", fixed = TRUE) ' , 
    " + ",
    paste(paste0("PARO2059.", 2:4), collapse = " + " ),
    " + ",
    paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
    " + ",
    paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
    " + ",
    paste(paste0("ESTPRE.", 2:4), collapse = " + " )
  )
)

## To reduce the computation time we can state as inicial values of the 
# modes of the different hyperparameters the posterior estimates of these
# parameters obtained from other fitted models. In the case, you want to
# use this option, you should be sure the results of the previous models
# are not remove it. 
#hyper.mode.M1 <- pp.res.M1$summary.hyperpar$mode


#Estimation
Sys.time()
pp.res <- inla(formula = form, 
  family = rep('poisson', 4), 
  data = inla.stack.data(join.stack),
  control.predictor = list(A = inla.stack.A(join.stack), link = 1,
    compute = TRUE),
  E = inla.stack.data(join.stack)$e, 
  control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE)
  #, control.mode = list(theta=hyper.mode.M1)
)
Sys.time()


summary(pp.res)
t2 <-Sys.time()

ti <- t2-t1; 
print(
  paste0("Model 1 has used ", ti, " cpu time. And have finished at ", 
    Sys.time()))

rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, 
   controls_ppp, kid_ppp,lung_ppp, sto_ppp)

save.image("Results/Res-M1-ConFact.RData")

rm(list = ls())

########################################################

########################################################
#Model 2: Lin Eff covariates, no specific spatial effects
########################################################

for(p in 1:13){
  
  #Load stack data and others
  load("Data/Stack-LinRW-models.RData")
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- as.formula( 
    paste0("y ~ ", 
      '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
      f(spatial.field1, model = spde) + 
      f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
      f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
      f(cont.copy4, copy = "spatial.field1", fixed = TRUE) ' , 
      " + ",
      paste(paste0("var", p, ".", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("PARO2059.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTPRE.", 2:4), collapse = " + " )
    )
  )
  
  ## To reduce the computation time we can state as inicial values of the 
  # modes of the different hyperparameters the posterior estimates of these
  # parameters obtained from other fitted models. In the case, you want to
  # use this option, you should be sure the results of the previous models
  # are not remove it. 
  #hyper.mode.M0 <- pp.res.M0$summary.hyperpar$mode
  
  #Estimation
  Sys.time()
  pp.res <- inla(formula = form, 
    family = rep('poisson', 4), data = inla.stack.data(join.stack),
    control.predictor = list(A = inla.stack.A(join.stack), link = 1,
      compute = TRUE),
    E = inla.stack.data(join.stack)$e, 
    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE)
    #,control.mode = list(theta=hyper.mode.M0)
  )
  Sys.time()
  
  summary(pp.res)
  
  t2 <-Sys.time()
  
  ti <- t2-t1; 
  print(
    paste0("Model with the covariable ", p, " has used ", ti, " cpu time. 
      And have finished at ", Sys.time()))
  
  rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, pp.res.M1, 
     pp.res.M0, controls_ppp, kid_ppp, lung_ppp, sto_ppp)
  
  guardar <- paste0("Results/M2-Var", p , "_LinEffect.RData" )
  
  save.image(guardar)
  
  rm(list = ls())
  
}

########################################################

########################################################
#Model 2: Scal. Rw1 covariates, no specific spatial effects
########################################################

for(p in 1:13){
  
  #Load stack data and others
  load("Data/Stack-LinRW-models.RData")
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- as.formula( 
    paste0("y ~ ", 
      '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
      f(spatial.field1, model = spde) + 
      f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
      f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
      f(cont.copy4, copy = "spatial.field1", fixed = TRUE) ' , 
      " + ",
      paste(paste0("PARO2059.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTPRE.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("f(inla.group(", paste0("var", p, ".", 2:4), ", n=50), 
        model='rw1', scale.model=TRUE)"), collapse = " + " )
    )
  )
  
  ## To reduce the computation time we can state as inicial values of the 
  # modes of the different hyperparameters the posterior estimates of these
  # parameters obtained from other fitted models. In the case, you want to
  # use this option, you should be sure the results of the previous models
  # are not remove it. 
  # hyper.mode.M0 <- pp.res.M0$summary.hyperpar$mode
  
  
  #Estimation
  pp.res <- inla(formula = form, 
    family = rep('poisson', 4), data = inla.stack.data(join.stack),
    control.predictor = list(A = inla.stack.A(join.stack), link = 1,
      compute = TRUE),
    E = inla.stack.data(join.stack)$e, 
    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE)
  )#,control.mode = list(theta=hyper.mode.M0))
  summary(pp.res)
  
  t2 <-Sys.time()
  
  ti <- t2-t1; 
  print(
    paste0("Model with the covariable ", p, " has used ", ti,
      " cpu time. And have finished at ", Sys.time()))
  
  rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, pp.res.M1, 
     pp.res.M0, controls_ppp, kid_ppp, lung_ppp, sto_ppp)
  
  guardar <- paste0("Results/M2_Var", p , "_RW1.RData" )
  
  save.image(guardar)
  
  rm(list = ls())
  
}

########################################################

########################################################
#Model 2: SPDE 1D covariates, no specific spatial effects
########################################################

for(p in 1:13){
  
  #Load stack data and others
  load("Data/Stack-SPDE1-Models.RData")
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- as.formula( 
    paste0("y ~ ", 
      '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
      f(spatial.field1, model = spde) + 
      f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
      f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
      f(cont.copy4, copy = "spatial.field1", fixed = TRUE) ' , 
      " + ",
      paste(paste0("f(", paste0("ind_var", p, ".", 2:4),
        ", model=spde_cov)"), collapse = " + " ),
      " + ",
      paste(paste0("PARO2059.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTPRE.", 2:4), collapse = " + " )
    )
  )
  
  ## To reduce the computation time we can state as inicial values of the 
  # modes of the different hyperparameters the posterior estimates of these
  # parameters obtained from other fitted models. In the case, you want to
  # use this option, you should be sure the results of the previous models
  # are not remove it. 
  # hyper.mode.M0 <- pp.res.M0$summary.hyperpar$mode
  
  #Estimation
  pp.res <- inla(formula = form, 
    family = rep('poisson', 4), data = inla.stack.data(join.stack),
    control.predictor = list(A = inla.stack.A(join.stack), link = 1,
      compute = TRUE),
    E = inla.stack.data(join.stack)$e, 
    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE)
  )#,control.mode = list(theta=hyper.mode.M0))
  Sys.time()
  summary(pp.res)
  
  t2 <-Sys.time()
  
  ti <- t2-t1 
  print(
    paste0("Model with the covariable ", p, " has used ", ti,
           " cpu time. And have finished at ", Sys.time()))
  
  
  rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, 
     pp.res.M1, pp.res.M0, controls_ppp, kid_ppp, lung_ppp, sto_ppp)
  
  guardar <- paste0("Results/M2_Var", p , "_SPDE.RData" )
  
  save.image(guardar)
  
  rm(list = ls())
  
}

########################################################

########################################################
#Model 3: Lin Eff covariates, with specific spatial effects
########################################################

for(p in 1:13){
  
  #Load stack data and others
  load("Data/Stack-LinRW-models.RData")
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- as.formula( 
    paste0("y ~ ", 
      '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
      f(spatial.field1, model = spde) + 
      f(spatial.field2, model = spde) + 
      f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
      f(spatial.field3, model = spde) +  
      f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
      f(spatial.field4, model = spde) + 
      f(cont.copy4, copy = "spatial.field1", fixed = TRUE)  ' , 
      " + ",
      paste(paste0("var", p, ".", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("PARO2059.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTPRE.", 2:4), collapse = " + " )
    )
  )
  
  ## To reduce the computation time we can state as inicial values of the 
  # modes of the different hyperparameters the posterior estimates of these
  # parameters obtained from other fitted models. In the case, you want to
  # use this option, you should be sure the results of the previous models
  # are not remove it. 
  #hyper.mode.M1 <- pp.res.M1$summary.hyperpar$mode
  
  #Estimation
  pp.res <- inla(formula = form, 
    family = rep('poisson', 4), data = inla.stack.data(join.stack),
    control.predictor = list(A = inla.stack.A(join.stack), link = 1,
      compute = TRUE),
    E = inla.stack.data(join.stack)$e, 
    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE),
  )#,control.mode = list(theta=hyper.mode.M1)

  summary(pp.res)
  
  t2 <-Sys.time()
  
  ti <- t2-t1; 
  print(
    paste0("Model with the covariable ", p, " has used ", ti, 
           " cpu time. And have finished at ", Sys.time()))
  
  
  rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, pp.res.M1,
     pp.res.M0, controls_ppp, kid_ppp, lung_ppp, sto_ppp)
  
  guardar <- paste0("Results/M3-Var", p , "_LinEffect.RData" )
  
  save.image(guardar)
  
  rm(list = ls())
  
}

########################################################

########################################################
#Model 3: Scal. Rw1 covariates, with specific spatial effects
########################################################

for(p in 1:13){
  
  #Load stack data and others
  load("Data/Stack-LinRW-models.RData")
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- as.formula( 
    paste0("y ~ ", 
      '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
      f(spatial.field1, model = spde) + 
      f(spatial.field2, model = spde) + 
      f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
      f(spatial.field3, model = spde) +  
      f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
      f(spatial.field4, model = spde) + 
      f(cont.copy4, copy = "spatial.field1", fixed = TRUE)  ' , 
      " + ",
      paste(paste0("PARO2059.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTPRE.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("f(inla.group(", paste0("var", p, ".", 2:4), 
        ", n=50), model='rw1', scale.model=TRUE)"),  collapse = " + " )
    )
  )
  
  ## To reduce the computation time we can state as inicial values of the 
  # modes of the different hyperparameters the posterior estimates of these
  # parameters obtained from other fitted models. In the case, you want to
  # use this option, you should be sure the results of the previous models
  # are not remove it. 
  #hyper.mode.M1 <- pp.res.M1$summary.hyperpar$mode
  
  
  #Estimation
  pp.res <- inla(formula = form, 
    family = rep('poisson', 4), data = inla.stack.data(join.stack),
    control.predictor = list(A = inla.stack.A(join.stack), link = 1,
      compute = TRUE),
    E = inla.stack.data(join.stack)$e, 
    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE)
  )#,control.mode = list(theta=hyper.mode.M1))
  summary(pp.res)
  
  t2 <-Sys.time()
  
  ti <- t2-t1; 
  print(
    paste0("Model with the covariable ", p, " has used ", ti,
      " cpu time. And have finished at ", Sys.time()))
  
  rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, pp.res.M1, 
     pp.res.M0, controls_ppp, kid_ppp, lung_ppp, sto_ppp)
  
  guardar <- paste0("Results/M3_Var", p , "_RW1.RData" )
  
  save.image(guardar)
  
  rm(list = ls())
  
}

########################################################

########################################################
#Model 3: SPDE 1D covariates, with specific spatial effects
########################################################

for(p in 1:13){
  
  #Load stack data and others
  load("Data/Stack-SPDE1-Models.RData")
  
  t1 <-Sys.time()
  E <- inla.stack.data(join.stack)$e
  
  #formula
  form <- as.formula( 
    paste0("y ~ ", 
      '0 + Intercept.1 + Intercept.2 + Intercept.3 + Intercept.4 + 
      f(spatial.field1, model = spde) + 
      f(spatial.field2, model = spde) + 
      f(cont.copy2, copy = "spatial.field1", fixed = TRUE)  + 
      f(spatial.field3, model = spde) +  
      f(cont.copy3, copy = "spatial.field1", fixed = TRUE)  + 
      f(spatial.field4, model = spde) + 
      f(cont.copy4, copy = "spatial.field1", fixed = TRUE)  ' , 
      " + ",
      paste(paste0("f(", paste0("ind_var", p, ".", 2:4),
        ", model=spde_cov)"), collapse = " + " ),
      " + ",
      paste(paste0("PARO2059.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("SOCMEDIA.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTMEDIO.", 2:4), collapse = " + " ),
      " + ",
      paste(paste0("ESTPRE.", 2:4), collapse = " + " )
    )
  )
  
  ## To reduce the computation time we can state as inicial values of the 
  # modes of the different hyperparameters the posterior estimates of these
  # parameters obtained from other fitted models. In the case, you want to
  # use this option, you should be sure the results of the previous models
  # are not remove it. 
  # hyper.mode.M1 <- pp.res.M1$summary.hyperpar$mode
  
  
  #Estimation
  pp.res <- inla(formula = form, verbose = F,
    family = rep('poisson', 4), data = inla.stack.data(join.stack),
    control.predictor = list(A = inla.stack.A(join.stack), link = 1,
      compute = TRUE),
    E = inla.stack.data(join.stack)$e, 
    control.compute = list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE, po=TRUE)
  )#,control.mode = list(theta=hyper.mode.M1))
  summary(pp.res)
  
  t2 <-Sys.time()
  
  ti <- t2-t1; 
  print(
    paste0("Model with the covariable ", p, " has used ", ti,
      " cpu time. And have finished at ", Sys.time()))
  
  
  rm(t1, t2, E, s.index1, form, s.index2, s.index3, s.index4, pp.res.M1,
     pp.res.M0, controls_ppp, kid_ppp, lung_ppp, sto_ppp)
  
  guardar <- paste0("Results/M3_Var", p , "_SPDE.RData" )
  
  save.image(guardar)
  
  rm(list = ls())
  
}

########################################################




