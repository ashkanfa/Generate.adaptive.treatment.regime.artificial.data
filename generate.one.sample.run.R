no.rows=case$n[2]
no.bin =(case$covar[2]) - 1
no.nor = 1
prop.vec.bin = round(runif(no.bin,min=0.2,max=0.8),digits = 1)
mean.vec.nor=0
var.nor=1
#make this constant, this is not the focus of the study

corr.vec = cor.vec.generate(p=case$covar[2], true.p=0.5, cor1=case$cor3[2],cor2=case$cor2[2],cor3=case$cor1[2])
no.trt=case$trt[2]
no.stage=3
bincorr.vec = cor.vec.generate(p=case$trt[2], true.p=case$true.trt.prop[2], cor1=case$cor1[2],cor2=case$cor2[2],cor3=case$cor3[2])
true.trt.prop=case$true.trt.prop[2]
true.cov.prop=case$true.cov.prop[2]
min.coef=case$min.coef[2]
max.coef=case$max.coef[2]
no.int.co=case$no.int.co[2]
no.int.cot=case$no.int.cot[2]
no.int.t=case$no.int.t[2] 

set.seed(10)
sim.dataaa<-NULL
sim.dataaa<-generate.adaptive.data(no.rows,no.bin, no.nor, prop.vec.bin, mean.vec.nor,var.no, 
                                      corr.vec, no.trt,no.stage,bincorr.vec, true.trt.prop,true.cov.prop,
                                      min.coef, max.coef, no.int.co, no.int.cot, no.int.t)
