### R code from vignette source 'trigger.Rnw'

###################################################
### code chunk number 1: import
###################################################
library(trigger)
 data(yeast)
 names(yeast)

#reduce data size for vignette run time
set.seed(123)
#select subset of 400 traits
gidx = sort(sample(1:6216, size = 400))
yeast$exp = yeast$exp[gidx,]
yeast$exp.pos = yeast$exp.pos[gidx,]
#select subset of markers
midx = sort(sample(1:3244, size = 500))
yeast$marker = yeast$marker[midx,]
yeast$marker.pos = yeast$marker.pos[midx,]

attach(yeast)
dim(exp)


###################################################
### code chunk number 2: build
###################################################
trig.obj <- trigger.build(marker=marker, exp=exp, 
	marker.pos=marker.pos, exp.pos=exp.pos)
trig.obj
detach(yeast)


###################################################
### code chunk number 3: link
###################################################
trig.obj = trigger.link(trig.obj, norm = TRUE)


###################################################
### code chunk number 4: fig1
###################################################
plot(trig.obj, type = "link", cutoff = 1e-5)


###################################################
### code chunk number 5: trigger.Rnw:163-164
###################################################
plot(trig.obj, type = "link", cutoff = 1e-5)


###################################################
### code chunk number 6: mlink
###################################################
trig.obj = trigger.loclink(trig.obj)
trig.obj = trigger.mlink(trig.obj, B = 10, idx = NULL)


###################################################
### code chunk number 7: fig2
###################################################
plot(trig.obj, type = "mlink", qcut = 0.2, bin.size = 50000)


###################################################
### code chunk number 8: trigger.Rnw:191-192
###################################################
plot(trig.obj, type = "mlink", qcut = 0.2, bin.size = 50000)


###################################################
### code chunk number 9: eigenR2
###################################################
trig.obj = trigger.eigenR2(trig.obj, adjust = FALSE)


###################################################
### code chunk number 10: fig3
###################################################
plot(trig.obj, type = "eigenR2")


###################################################
### code chunk number 11: trigger.Rnw:212-213
###################################################
plot(trig.obj, type = "eigenR2")


###################################################
### code chunk number 12: trigger.net
###################################################
trig.obj = trigger.loclink(trig.obj, window.size = 10000)
trig.prob = trigger.net(trig.obj, Bsec = 100, idx = NULL) 
 dim(trig.prob)


###################################################
### code chunk number 13: detach
###################################################
# Re-attach dataset to re-include all traits for analysis
data(yeast); attach(yeast)
dim(exp)
trig.obj <- trigger.build(marker=marker, exp=exp, 
	marker.pos=marker.pos, exp.pos=exp.pos)


###################################################
### code chunk number 14: export2cross
###################################################
cross = trigger.export2cross(trig.obj, plotarg = FALSE, verbose = FALSE)


###################################################
### code chunk number 15: trait-trigger
###################################################
causreg = trigger.trait(trig.obj, trait = "DSE1", cross = cross, addplot = TRUE, thr = 3)
causreg


