library(ggplot2)
library(ComplexUpset)
library("hgu133a2.db")
library(UpSetR)
library(Seurat)
library(readxl)

########## Controlled Ovarian Stimulation ###########

## GnRH-a & GnRH-anta

Devjak <- read_excel("Data COH.xlsx", 
                     sheet = "Devjak 2012")

## HP-hMG & rFSH CC

Adriaenssens <- read_excel("Data COH.xlsx", 
                           sheet = "Adriaenssens 2010")
AdriaenssensUP=Adriaenssens[which(Adriaenssens$Direction=="Up HP-hMG"),]
AdriaenssensDOWN=Adriaenssens[which(Adriaenssens$Direction=="Down HP-hMG"),]

Assou <- read_excel("Data COH.xlsx", 
                    sheet = "Assou 2013")
AssouUP=Assou[which(Assou$Direction=="Up HP-hMG"),]
AssouDOWN=Assou[which(Assou$Direction=="Down HP-hMG"),]

Cruz1 <- read_excel("Data COH.xlsx", 
                    sheet = "Cruz 2017 rFSH hMG")
CruzUP1=Cruz1[which(Cruz1$Direction=="Up HP-hMG"),]
CruzDOWN1=Cruz1[which(Cruz1$Direction=="Down HP-hMG"),]

Gatta1 <- read_excel("Data COH.xlsx", 
                     sheet = "Gatta 2013 HP-hMG rFSH")
GattaUP1=Gatta1[which(Gatta1$Direction=="Up HP-hMG"),]
GattaDOWN1=Gatta1[which(Gatta1$Direction=="Down HP-hMG"),]

## HP-hMG & urinary FSH CC

Cruz2 <- read_excel("Data COH.xlsx", 
                    sheet = "Cruz 2017 urinary FSH hMG")

## urinary FSH & rFSH CC

Cruz3 <- read_excel("Data COH.xlsx", 
                    sheet = "Cruz 2017 urinary FSH rFSH")

## rLH+rFSH & rFSH CC

Gatta2 <- read_excel("Data COH.xlsx", 
                     sheet = "Gatta 2013 rLH+rFSH rFSH")

Barbieri <- read_excel("Data COH.xlsx", 
                       sheet = "Barbieri 2012")

## FSH & rFSH CC

Gurgan1 <- read_excel("Data COH.xlsx", 
                      sheet = "Gurgan 2014 FSH rFSH")

## FSH+rFSH & FSH CC

Gurgan2 <- read_excel("Data COH.xlsx", 
                      sheet = "Gurgan 2014 FSH+rFSH FSH")

## FSH+rFSH & rFSH CC

Gurgan3 <- read_excel("Data COH.xlsx", 
                      sheet = "Gurgan 2014 FSH+rFSH rFSH")

## HP-hMG & rFSH MGC

Brannian <- read_excel("Data COH.xlsx", 
                       sheet = "Brannian 2010")
BrannianUP=Brannian[which(Brannian$Direction=="Up HP-hMG"),]
BrannianDOWN=Brannian[which(Brannian$Direction=="Down HP-hMG"),]

Grondahl <- read_excel("Data COH.xlsx", 
                       sheet = "Grondahl 2009")
k <- Grondahl$`Probe set`
Genes=select(hgu133a2.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="PROBEID")
Grondahl$Genes=Genes$SYMBOL

GrondahlUP=Grondahl[which(Grondahl$Direction=="Up HP-hMG"),]
GrondahlDOWN=Grondahl[which(Grondahl$Direction=="Down HP-hMG"),]



## Mild vs COS MGC

Liu <- read_excel("Data COH.xlsx", 
                  sheet = "Liu 2022")

### ALL

"Devjak et al. (2012)"=na.omit(unique(Devjak$Genes))
"Adriaenssens et al. (2010) UP"=UpdateSymbolList(na.omit(unique(AdriaenssensUP$Genes)),timeout = 10000)
"Adriaenssens et al. (2010) DOWN"=UpdateSymbolList(na.omit(unique(AdriaenssensDOWN$Genes)),timeout = 10000)

"Assou et al. (2013) UP"=UpdateSymbolList(na.omit(unique(AssouUP$Genes)),timeout = 10000)
"Assou et al. (2013) DOWN"=UpdateSymbolList(na.omit(unique(AssouDOWN$Genes)),timeout = 10000)

"Cruz et al. (2017) UP a"=UpdateSymbolList(na.omit(unique(CruzUP1$Genes)),timeout = 10000)
"Cruz et al. (2017) DOWN a"=UpdateSymbolList(na.omit(unique(CruzDOWN1$Genes)),timeout = 10000)

"Gatta et al. (2013) UP a"=UpdateSymbolList(na.omit(unique(GattaUP1$Genes)),timeout = 10000)
"Gatta et al. (2013) DOWN a"=UpdateSymbolList(na.omit(unique(GattaDOWN1$Genes)),timeout = 10000)

"Cruz et al. (2017) a"=c(`Cruz et al. (2017) UP a`,`Cruz et al. (2017) DOWN a`)
"Cruz et al. (2017) b"=UpdateSymbolList(na.omit(unique(Cruz2$Genes)),timeout = 10000)
"Cruz et al. (2017) c"=UpdateSymbolList(na.omit(unique(Cruz3$Genes)),timeout = 10000)

"Gatta et al. (2013) a"=c(`Gatta et al. (2013) UP a`,`Gatta et al. (2013) DOWN a`)
"Gatta et al. (2013) b"=UpdateSymbolList(na.omit(unique(Gatta2$Genes)),timeout = 10000)

"Barberi et al. (2012)"=UpdateSymbolList(na.omit(unique(Barbieri$Genes)),timeout = 10000)

"Gurgan et al. (2014) a"=UpdateSymbolList(na.omit(unique(Gurgan1$Genes)),timeout = 10000)
"Gurgan et al. (2014) b"=UpdateSymbolList(na.omit(unique(Gurgan2$Genes)),timeout = 10000)
"Gurgan et al. (2014) c"=UpdateSymbolList(na.omit(unique(Gurgan3$Genes)),timeout = 10000)

## HP-hMG & rFSH CC UP

HPhMGrFSHUP=list("Adriaenssens et al. (2010)"=`Adriaenssens et al. (2010) UP`,
               "Assou et al. (2013)"=`Assou et al. (2013) UP`,
               "Cruz et al. (2017)"=`Cruz et al. (2017) UP a`,
               "Gatta et al. (2013)"=`Gatta et al. (2013) UP a`)
HPhMGrFSHUP=rev(HPhMGrFSHUP)
HPhMGrFSHUPintersectionlist=overlapGroups(HPhMGrFSHUP)
HPhMGrFSHUPintersectionlist <- purrr::map(HPhMGrFSHUPintersectionlist, ~ attr(HPhMGrFSHUPintersectionlist, "elements")[.x] )
HPhMGrFSHUP=fromList(HPhMGrFSHUP)
ComplexUpset::upset(HPhMGrFSHUP,
                    colnames(HPhMGrFSHUP),sort_sets=FALSE)

## HP-hMG & rFSH CC DOWN

HPhMGrFSHDOWN=list("Adriaenssens et al. (2010)"=`Adriaenssens et al. (2010) DOWN`,
                 "Assou et al. (2013)"=`Assou et al. (2013) DOWN`,
                 "Cruz et al. (2017)"=`Cruz et al. (2017) DOWN a`,
                 "Gatta et al. (2013)"=`Gatta et al. (2013) DOWN a`)
HPhMGrFSHDOWN=rev(HPhMGrFSHDOWN)
HPhMGrFSHDOWNintersectionlist=overlapGroups(HPhMGrFSHDOWN)
HPhMGrFSHDOWNintersectionlist <- purrr::map(HPhMGrFSHDOWNintersectionlist, ~ attr(HPhMGrFSHDOWNintersectionlist, "elements")[.x] )
HPhMGrFSHDOWN=fromList(HPhMGrFSHDOWN)
ComplexUpset::upset(HPhMGrFSHDOWN,
                    colnames(HPhMGrFSHDOWN),sort_sets=FALSE)

## Cruz comparison

Cruz=list("Cruz et al. (2017) a"=`Cruz et al. (2017) a`,
          "Cruz et al. (2017) b"=`Cruz et al. (2017) b`,
          "Cruz et al. (2017) c"=`Cruz et al. (2017) c`)
Cruz=rev(Cruz)
Cruzintersectionlist=overlapGroups(Cruz)
Cruzintersectionlist <- purrr::map(Cruzintersectionlist, ~ attr(Cruzintersectionlist, "elements")[.x] )
Cruz=fromList(Cruz)
ComplexUpset::upset(Cruz,
                    colnames(Cruz),sort_sets=FALSE)

ggVennDiagram(list(`Cruz et al. (2017) a`,
                   `Cruz et al. (2017) b`,
                   `Cruz et al. (2017) c`),label="count",
              category.names = c("A","B","C"))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

## Gatta comparison

Gatta=list("Gatta et al. (2013) a"=`Gatta et al. (2013) a`,
           "Gatta et al. (2013) b"=`Gatta et al. (2013) b`)
Gatta=rev(Gatta)
Gattaintersectionlist=overlapGroups(Gatta)
Gattaintersectionlist <- purrr::map(Gattaintersectionlist, ~ attr(Gattaintersectionlist, "elements")[.x] )
Gatta=fromList(Gatta)
ComplexUpset::upset(Gatta,
                    colnames(Gatta),sort_sets=FALSE)

ggVennDiagram(list(`Gatta et al. (2013) a`,
                   `Gatta et al. (2013) b`),label="count",
              category.names = c("A","B"))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

## Gurgan comparison

Gurgan=list("Gurgan et al. (2014) a"=`Gurgan et al. (2014) a`,
            "Gurgan et al. (2014) b"=`Gurgan et al. (2014) b`,
            "Gurgan et al. (2014) c"=`Gurgan et al. (2014) c`)
Gurgan=rev(Gurgan)
Gurganintersectionlist=overlapGroups(Gurgan)
Gurganintersectionlist <- purrr::map(Gurganintersectionlist, ~ attr(Gurganintersectionlist, "elements")[.x] )
Gurgan=fromList(Gurgan)
ComplexUpset::upset(Gurgan,
                    colnames(Gurgan),sort_sets=FALSE)

ggVennDiagram(list(`Gurgan et al. (2014) a`,
                   `Gurgan et al. (2014) b`,
                   `Gurgan et al. (2014) c`),label="count",
              category.names = c("A","B","C"))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

## rLH+rFSH & rFSH CC

rLHrFSH=list("Barberi et al. (2012)"=`Barberi et al. (2012)`,
             "Gatta et al. (2013)"=`Gatta et al. (2013) b`)
rLHrFSH=rev(rLHrFSH)
rLHrFSHintersectionlist=overlapGroups(rLHrFSH)
rLHrFSHintersectionlist <- purrr::map(rLHrFSHintersectionlist, ~ attr(rLHrFSHintersectionlist, "elements")[.x] )
rLHrFSH=fromList(rLHrFSH)
ComplexUpset::upset(rLHrFSH,
                    colnames(rLHrFSH),sort_sets=FALSE)


## HP-hMG MGC

`Brannian et al. (2010) UP`=UpdateSymbolList(na.omit(unique(BrannianUP$Genes)),timeout = 10000)
`Brannian et al. (2010) DOWN`=UpdateSymbolList(na.omit(unique(BrannianDOWN$Genes)),timeout = 10000)

`Grondahl et al. (2009) UP`=UpdateSymbolList(na.omit(unique(GrondahlUP$Genes)),timeout = 10000)
`Grondahl et al. (2009) DOWN`=UpdateSymbolList(na.omit(unique(GrondahlDOWN$Genes)),timeout = 10000)

MGCUP=list("Brannian et al. (2010)"=`Brannian et al. (2010) UP`,
         "Grondahl et al. (2009)"=`Grondahl et al. (2009) UP`)

MGCUP=rev(MGCUP)
MGCUPintersectionlist=overlapGroups(MGCUP)
MGCUPintersectionlist <- purrr::map(MGCUPintersectionlist, ~ attr(MGCUPintersectionlist, "elements")[.x] )
MGCUP=fromList(MGCUP)
ComplexUpset::upset(MGCUP,
                    colnames(MGCUP),sort_sets=FALSE)

ggVennDiagram(list(`Brannian et al. (2010) UP`,
                   `Grondahl et al. (2009) UP`),label="count",
              category.names = c("B","G"))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

MGCDOWN=list("Brannian et al. (2010)"=`Brannian et al. (2010) DOWN`,
           "Grondahl et al. (2009)"=`Grondahl et al. (2009) DOWN`)

MGCDOWN=rev(MGCDOWN)
MGCDOWNintersectionlist=overlapGroups(MGCDOWN)
MGCDOWNintersectionlist <- purrr::map(MGCDOWNintersectionlist, ~ attr(MGCDOWNintersectionlist, "elements")[.x] )
MGCDOWN=fromList(MGCDOWN)
ComplexUpset::upset(MGCDOWN,
                    colnames(MGCDOWN),sort_sets=FALSE)

ggVennDiagram(list(`Brannian et al. (2010) DOWN`,
                   `Grondahl et al. (2009) DOWN`),label="count",
              category.names = c("B","G"))+ 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")


### Trigger

BorgboCC <- read_excel("Data COH.xlsx", 
                       sheet = "Borgbo 2013 CC")
BorgboCCUP = BorgboCC[which(BorgboCC$`Fold-Change(agonist vs. hCG) (Description)`=="Up hCG"),]
BorgboCCDOWN = BorgboCC[which(BorgboCC$`Fold-Change(agonist vs. hCG) (Description)`=="Down hCG"),]


BorgboMGC <- read_excel("Data COH.xlsx", 
                        sheet = "Borgbo 2013 MGC")
BorgboMGCUP = BorgboMGC[which(BorgboMGC$`Fold-Change(agonist vs. hCG)...7`=="Up hCG"),]
BorgboMGCDOWN = BorgboMGC[which(BorgboMGC$`Fold-Change(agonist vs. hCG)...7`=="Down hCG"),]


Haas1 <- read_excel("Data COH.xlsx", 
                    sheet = "Haas 2014")
HaasUP1=Haas1[which(Haas1$Direction=="Up hCG"),]
HaasDOWN1=Haas1[which(Haas1$Direction=="Down hCG"),]

Haas2 <- read_excel("Data COH.xlsx", 
                    sheet = "Haas 2016")

HaasUP2=Haas1[which(Haas2$Direction=="Up double"),]
HaasDOWN2=Haas1[which(Haas2$Direction=="Down double"),]

Vuong <- read_excel("Data COH.xlsx", 
                    sheet = "Vuong 2017")

Weizman1 <- read_excel("Data COH.xlsx", 
                       sheet = "Weizman 2019 High responder")
WeizmanHighUP = Weizman1[which(Weizman1$Direction=="Up dual"),]
WeizmanHighDOWN = Weizman1[which(Weizman1$Direction=="Down dual"),]

Weizman2 <- read_excel("Data COH.xlsx", 
                       sheet = "Weizman 2019 Normal responder")
WeizmanNormalUP = Weizman2[which(Weizman2$Direction=="Up dual"),]
WeizmanNormalDOWN = Weizman2[which(Weizman2$Direction=="Down dual"),]

Weizman3 <- read_excel("Data COH.xlsx", 
                       sheet = "Weizman 2019 Poor responder")
WeizmanPoorUP = Weizman3[which(Weizman3$Direction=="Up dual"),]
WeizmanPoorDOWN = Weizman3[which(Weizman3$Direction=="Down dual"),]

`Borgbo et al. (2013) CC UP`=UpdateSymbolList(na.omit(unique(BorgboCCUP$Genes)),timeout = 10000)
`Borgbo et al. (2013) CC DOWN`=UpdateSymbolList(na.omit(unique(BorgboCCDOWN$Genes)),timeout = 10000)

`Borgbo et al. (2013) MGC UP`=UpdateSymbolList(na.omit(unique(BorgboMGCUP$Genes)),timeout = 10000)
`Borgbo et al. (2013) MGC DOWN`=UpdateSymbolList(na.omit(unique(BorgboMGCDOWN$Genes)),timeout = 10000)

`Haas et al. (2014) UP`=UpdateSymbolList(na.omit(unique(HaasUP1$Genes)),timeout = 10000)
`Haas et al. (2014) DOWN`=UpdateSymbolList(na.omit(unique(HaasDOWN1$Genes)),timeout = 10000)

`Haas et al. (2016) UP`=UpdateSymbolList(na.omit(unique(HaasUP2$Genes)),timeout = 10000)
`Haas et al. (2016) DOWN`=UpdateSymbolList(na.omit(unique(HaasDOWN2$Genes)),timeout = 10000)

`Vuong et al. (2017)`=na.omit(unique(Vuong$Genes))
`Weizman et al. (2019) High UP`=UpdateSymbolList(na.omit(unique(WeizmanHighUP$Genes)),timeout = 10000)
`Weizman et al. (2019) High DOWN`=UpdateSymbolList(na.omit(unique(WeizmanHighDOWN$Genes)),timeout = 10000)

`Weizman et al. (2019) Normal UP`=UpdateSymbolList(na.omit(unique(WeizmanNormalUP$Genes)),timeout = 10000)
`Weizman et al. (2019) Normal DOWN`=UpdateSymbolList(na.omit(unique(WeizmanNormalDOWN$Genes)),timeout = 10000)

`Weizman et al. (2019) Poor UP`=UpdateSymbolList(na.omit(unique(WeizmanPoorUP$Genes)),timeout = 10000)
`Weizman et al. (2019) Poor DOWN`=UpdateSymbolList(na.omit(unique(WeizmanPoorDOWN$Genes)),timeout = 10000)

#`Weizman et al. (2019)`=UpdateSymbolList(na.omit(unique(Weizman)),timeout = 10000)

hCGUP=list(
  "Borgbo et al. (2013) CC"=`Borgbo et al. (2013) CC UP`,
  "Borgbo et al. (2013) MGC"=`Borgbo et al. (2013) MGC UP`,
  "Haas et al. (2014)"=`Haas et al. (2014) UP`)

hCGUP=rev(hCGUP)
hCGUPintersectionlist=overlapGroups(hCGUP)
hCGUPintersectionlist <- purrr::map(hCGUPintersectionlist, ~ attr(hCGUPintersectionlist, "elements")[.x] )
hCGUP=fromList(hCGUP)
ComplexUpset::upset(hCGUP,
                    colnames(hCGUP),sort_sets=FALSE)

hCGDOWN=list(
  "Borgbo et al. (2013) CC"=`Borgbo et al. (2013) CC DOWN`,
  "Borgbo et al. (2013) MGC"=`Borgbo et al. (2013) MGC DOWN`,
  "Haas et al. (2014)"=`Haas et al. (2014) DOWN`)

hCGDOWN=rev(hCGDOWN)
hCGDOWNintersectionlist=overlapGroups(hCGDOWN)
hCGDOWNintersectionlist <- purrr::map(hCGDOWNintersectionlist, ~ attr(hCGDOWNintersectionlist, "elements")[.x] )
hCGDOWN=fromList(hCGDOWN)
ComplexUpset::upset(hCGDOWN,
                    colnames(hCGDOWN),sort_sets=FALSE)

DualGnRHUP=list(
  "Haas et al. (2016)"=`Haas et al. (2016) UP`,
  "Weizman et al. (2019) normal"=`Weizman et al. (2019) Normal UP`,
  "Weizman et al. (2019) poor"=`Weizman et al. (2019) Poor UP`)

DualGnRHUP=rev(DualGnRHUP)
DualGnRHUPintersectionlist=overlapGroups(DualGnRHUP)
DualGnRHUPintersectionlist <- purrr::map(DualGnRHUPintersectionlist, ~ attr(DualGnRHUPintersectionlist, "elements")[.x] )
DualGnRHUP=fromList(DualGnRHUP)
ComplexUpset::upset(DualGnRHUP,
                    colnames(DualGnRHUP),sort_sets=FALSE)

DualGnRHDOWN=list(
  "Haas et al. (2016)"=`Haas et al. (2016) DOWN`,
  "Weizman et al. (2019) normal"=`Weizman et al. (2019) Normal DOWN`,
  "Weizman et al. (2019) poor"=`Weizman et al. (2019) Poor DOWN`)

DualGnRHDOWN=rev(DualGnRHDOWN)
DualGnRHDOWNintersectionlist=overlapGroups(DualGnRHDOWN)
DualGnRHDOWNintersectionlist <- purrr::map(DualGnRHDOWNintersectionlist, ~ attr(DualGnRHDOWNintersectionlist, "elements")[.x] )
DualGnRHDOWN=fromList(DualGnRHDOWN)
ComplexUpset::upset(DualGnRHDOWN,
                    colnames(DualGnRHDOWN),sort_sets=FALSE)


### COH itself


`de los Santos` <- read_excel("Data COH.xlsx", 
                              sheet = "de los Santos 2012")
`de los Santos UP` = `de los Santos`[which(`de los Santos`$Direction=="Up COS"),]$Genes
`de los Santos DOWN` = `de los Santos`[which(`de los Santos`$Direction=="Down COS"),]$Genes


Lu <- read_excel("Data COH.xlsx", 
                 sheet = "Lu 2019")
LuUP = Lu[which(Lu$Direction=="Up COS"),]$Genes
LuDOWN = Lu[which(Lu$Direction=="Down COS"),]$Genes

Papler <- read_excel("Data COH.xlsx", 
                     sheet = "Papler 2013")
PaplerUP = Papler[which(Papler$Direction=="Up COS"),]$Genes
PaplerDOWN = Papler[which(Papler$Direction=="Down COS"),]$Genes

`de los Santos et al. (2012)`=UpdateSymbolList(na.omit(unique(`de los Santos`$Genes)),timeout = 10000)
`de los Santos et al. (2012) UP`=UpdateSymbolList(na.omit(unique(`de los Santos UP`)),timeout = 10000)
`de los Santos et al. (2012) DOWN`=UpdateSymbolList(na.omit(unique(`de los Santos DOWN`)),timeout = 10000)

`Lu et al. (2019)`=UpdateSymbolList(na.omit(unique(Lu$Genes)),timeout = 10000)
`Lu et al. (2019) UP`=UpdateSymbolList(na.omit(unique(LuUP)),timeout = 10000)
`Lu et al. (2019) DOWN`=UpdateSymbolList(na.omit(unique(LuDOWN)),timeout = 10000)

`Papler et al. (2013)`=UpdateSymbolList(na.omit(unique(Papler$Genes)),timeout = 10000)
`Papler et al. (2013) UP`=UpdateSymbolList(na.omit(unique(PaplerUP)),timeout = 10000)
`Papler et al. (2013) DOWN`=UpdateSymbolList(na.omit(unique(PaplerDOWN)),timeout = 10000)


COHUP=list("de los Santos et al. (2012)"=`de los Santos et al. (2012) UP`,
         "Lu et al. (2019)"=`Lu et al. (2019) UP`,
         "Papler et al. (2013)"=`Papler et al. (2013) UP`)

COHUP=rev(COHUP)
COHUPintersectionlist=overlapGroups(COHUP)
COHUPintersectionlist <- purrr::map(COHUPintersectionlist, ~ attr(COHUPintersectionlist, "elements")[.x] )
COHUP=fromList(COHUP)
ComplexUpset::upset(COHUP,
                    colnames(COHUP),sort_sets=FALSE)

COHDOWN=list("de los Santos et al. (2012)"=`de los Santos et al. (2012) DOWN`,
           "Lu et al. (2019)"=`Lu et al. (2019) DOWN`,
           "Papler et al. (2013)"=`Papler et al. (2013) DOWN`)

COHDOWN=rev(COHDOWN)
COHDOWNintersectionlist=overlapGroups(COHDOWN)
COHDOWNintersectionlist <- purrr::map(COHDOWNintersectionlist, ~ attr(COHDOWNintersectionlist, "elements")[.x] )
COHDOWN=fromList(COHDOWN)
ComplexUpset::upset(COHDOWN,
                    colnames(COHDOWN),sort_sets=FALSE)


########## In Vitro Maturation ###########


## rescue IVM

Jones <- read_excel("Data IVM.xlsx", 
                    sheet = "Jones 2008")
JonesUP = Jones[which(Jones$Direction=="Up IVM"),]
JonesDOWN = Jones[which(Jones$Direction=="Down IVM"),]

Lee <- read_excel("Data IVM.xlsx", 
                  sheet = "Lee 2021 in vivo vs IVM")
LeeUP = Lee[which(Lee$Direction=="Up IVM"),]
LeeDOWN = Lee[which(Lee$Direction=="Down IVM"),]

`Virant-Klun` <- read_excel("Data IVM.xlsx", 
                            sheet = "Virant-Klun 2018 top25")
`Virant-Klun UP` = `Virant-Klun`[which(`Virant-Klun`$Direction=="Up IVM"),]

Zhao <- read_excel("Data IVM.xlsx", 
                   sheet = "Zhao 2019")

#Li <- read_excel("Data IVM.xlsx", 
#                   sheet = "Li 2019")

Ouandaogo <- read_excel("Data IVM.xlsx", 
                        sheet = "Ouandaogo 2012")
OuandaogoUP = Ouandaogo[which(Ouandaogo$Direction=="Up IVM"),]
OuandaogoDOWN = Ouandaogo[which(Ouandaogo$Direction=="Down IVM"),]

## Medical IVM

Yang <- read_excel("Data IVM.xlsx", 
                   sheet = "Yang 2022")

Ye <- read_excel("Data IVM.xlsx", 
                 sheet = "Ye 2020")

Guzman <- read_excel("Data IVM.xlsx", 
                     sheet = "Guzman 2013 CC")
GuzmanDOWN = Guzman[which(Guzman$Direction=="Down IVM"),]

`Jones et al. (2008) UP`=UpdateSymbolList(na.omit(unique(JonesUP$Genes)),timeout = 10000)
`Jones et al. (2008) DOWN`=UpdateSymbolList(na.omit(unique(JonesDOWN$Genes)),timeout = 10000)
`Lee et al. (2021) UP`=UpdateSymbolList(na.omit(unique(LeeUP$Genes)),timeout = 10000)
`Lee et al. (2021) DOWN`=UpdateSymbolList(na.omit(unique(LeeDOWN$Genes)),timeout = 10000)
`Virant-Klun et al. (2018) top 25`=UpdateSymbolList(na.omit(unique(`Virant-Klun UP`$Genes)),timeout = 10000)
`Zhao et al. (2019)`=na.omit(unique(Zhao$Genes))
`Li et al. (2019)`=UpdateSymbolList(na.omit(unique(Li$Genes)),timeout = 10000)
`Ouandaogo et al. (2012) UP`=UpdateSymbolList(na.omit(unique(OuandaogoUP$Genes)),timeout = 10000)
`Ouandaogo et al. (2012) DOWN`=UpdateSymbolList(na.omit(unique(OuandaogoDOWN$Genes)),timeout = 10000)

`Guzman et al. (2013) DOWN`=UpdateSymbolList(na.omit(unique(GuzmanDOWN$Genes)),timeout = 10000)


rescueIVMUP=list(`Ouandaogo et al. (2012)`=`Ouandaogo et al. (2012) UP`,
                 `Lee et al. (2021)`=`Lee et al. (2021) UP`,
                 `Jones et al. (2008)`=`Jones et al. (2008) UP`,
                 `Virant-Klun et al. (2018) top 25`=`Virant-Klun et al. (2018) top 25`)
rescueIVMUP=rev(rescueIVMUP)
rescueIVMUPintersectionlist=overlapGroups(rescueIVMUP)
rescueIVMUPintersectionlist <- purrr::map(rescueIVMUPintersectionlist, ~ attr(rescueIVMUPintersectionlist, "elements")[.x] )
rescueIVMUP=fromList(rescueIVMUP)
ComplexUpset::upset(rescueIVMUP,
                    colnames(rescueIVMUP),sort_sets=FALSE)

rescueIVMDOWN=list(`Ouandaogo et al. (2012)`=`Ouandaogo et al. (2012) DOWN`,
                   `Lee et al. (2021)`=`Lee et al. (2021) DOWN`,
                   `Jones et al. (2008)`=`Jones et al. (2008) DOWN`)
rescueIVMDOWN=rev(rescueIVMDOWN)
rescueIVMDOWNintersectionlist=overlapGroups(rescueIVMDOWN)
rescueIVMDOWNintersectionlist <- purrr::map(rescueIVMDOWNintersectionlist, ~ attr(rescueIVMDOWNintersectionlist, "elements")[.x] )
rescueIVMDOWN=fromList(rescueIVMDOWN)
ComplexUpset::upset(rescueIVMDOWN,
                    colnames(rescueIVMDOWN),sort_sets=FALSE)


########## Cryoconservation ##########

Monzo <- read_excel("Data cryo.xlsx", 
                    sheet = "Monzo top10 down cryo")

Chamayou <- read_excel("Data cryo.xlsx", 
                       sheet = "Chamayou 2011")

Huo <- read_excel("Data cryo.xlsx", 
                  sheet = "Huo 2020")

Barberet <- read_excel("Data cryo.xlsx", 
                       sheet = "Barberet 2022")

HuoUP = Huo[which(Huo$Direction=="Up cryo"),]
HuoDOWN = Huo[which(Huo$Direction=="Down cryo"),]

BarberetUP = Barberet[which(Barberet$Direction=="Up cryo"),]
BarberetDOWN = Barberet[which(Barberet$Direction=="Down cryo"),]

Stigliani <- read_excel("Data cryo.xlsx", 
                        sheet = "Stigliani 2015")
StiglianiUP = Stigliani[which(Stigliani$Direction=="Up cryo"),]
StiglianiDOWN = Stigliani[which(Stigliani$Direction=="Down cryo"),]


`Monzo et al. (2012)`=UpdateSymbolList(na.omit(unique(Monzo$Genes)),timeout = 10000)
`Chamayou et al. (2011)`=UpdateSymbolList(na.omit(unique(Chamayou$Genes)),timeout = 10000)
`Huo et al. (2021) UP`=UpdateSymbolList(na.omit(unique(HuoUP$Genes)),timeout = 10000)
`Huo et al. (2021) DOWN`=UpdateSymbolList(na.omit(unique(HuoDOWN$Genes)),timeout = 10000)

`Barberet et al. (2022) UP`=UpdateSymbolList(na.omit(unique(BarberetUP$Genes)),timeout = 10000)
`Barberet et al. (2022) DOWN`=UpdateSymbolList(na.omit(unique(BarberetDOWN$Genes)),timeout = 10000)

`Stigliani et al. (2015) UP`=UpdateSymbolList(na.omit(unique(StiglianiUP$Genes)),timeout = 10000)
`Stigliani et al. (2015) DOWN`=UpdateSymbolList(na.omit(unique(StiglianiDOWN$Genes)),timeout = 10000)

CryoUP=list("Huo et al. (2021)"=`Huo et al. (2021) UP`,
            "Barberet et al. (2022)"=`Barberet et al. (2022) UP`,
            "Stigliani et al. (2015)"=`Stigliani et al. (2015) UP`)

CryoUP=rev(CryoUP)
CryoUPintersectionlist=overlapGroups(CryoUP)
CryoUPintersectionlist <- purrr::map(CryoUPintersectionlist, ~ attr(CryoUPintersectionlist, "elements")[.x] )
CryoUP=fromList(CryoUP)
ComplexUpset::upset(CryoUP,
                    colnames(CryoUP),sort_sets=FALSE)

CryoDOWN=list("Huo et al. (2021)"=`Huo et al. (2021) DOWN`,
              "Barberet et al. (2022)"=`Barberet et al. (2022) DOWN`,
              "Stigliani et al. (2015)"=`Stigliani et al. (2015) DOWN`,
              "Chamayou et al. (2011)"=`Chamayou et al. (2011)`,
              "Monzo et al. (2012)"=`Monzo et al. (2012)`)

CryoDOWN=rev(CryoDOWN)
CryoDOWNintersectionlist=overlapGroups(CryoDOWN)
CryoDOWNintersectionlist <- purrr::map(CryoDOWNintersectionlist, ~ attr(CryoDOWNintersectionlist, "elements")[.x] )
CryoDOWN=fromList(CryoDOWN)
ComplexUpset::upset(CryoDOWN,
                    colnames(CryoDOWN),sort_sets=FALSE)
