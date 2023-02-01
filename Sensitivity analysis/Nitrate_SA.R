

library(RxODE)
library(multisensi)
library(sensitivity)

rm(list = ls())


#建立抽样空间

library(dplyr)
library(sensitivity)


para_value<-read.csv("参数及其典型值.csv",header=T,sep=",")

para_list<-read.csv("参数列表.csv",header = T,sep=",")

para_value<-na.omit(para_value)

col_need<-seq(1,183,by=2)

para_value<-para_value[col_need,]

para_value<-na.omit(para_value)


para<-para_value[which(para_value$parameter%in%para_list$para),]

dim(para)


para$typical_value<-as.numeric(para$typical_value)

m <- 1000

Xb<-data.frame(matrix(nrow = m,ncol = 0))

for (i in 1:length(para$parameter)) {
  
  minvalue=0.1*para$typical_value[i]
  maxvalue=10*para$typical_value[i]  
  Xb=cbind(Xb,runif(m,min=minvalue,max = maxvalue))
}

#Xb<-data.frame(matrix(Xb,nrow = m,ncol = length(para$parameter)))

#Xb<-as.data.frame.matrix(Xb)



Xb<-apply(Xb, 2, as.numeric)

Xb<-as.data.frame(Xb)

names(Xb)<-para$parameter

#Xb<-data.frame(matrix(Xb))

X1 =Xb[1:(m/2), ]


X2 = Xb[(1 + m/2):m, ]

#输入模型
#%溶解氧(mg/L)#
# SO=c(6.6,7,6.8,6.7,7.5) 
SO=6.6
# #基质浓度(mg/L)#
# SS=c(9.69 ,7.93,7.33,6.92,6.69)
SS=9.69
# #进水UAP浓度##
#   SUAP=c(0.98,1.65,2.30,3.06,4.11)
SUAP=0.98
# #进水BAP浓度(mg/L)#
#   SBAP=c(0.89,0.91,0.93,0.95,0.95)
SBAP=0.89
# #氨氮浓度#
#   SNH=c(8.41,8.14,7.40,6.24,5.02)
SNH=8.41
# #亚硝氮硝氮浓度#
#   SNO=c(0.81,1.18,2.14,3.45,4.86)
SNO=0.81
# #碱度#
#   SALK=c(0.01,0.03,0.05,0.05,0.1)
SALK=0.01
# #异养菌(gX/m3)#
#   XH=c(5938.5832,4375.3694,1345.9254,266.19069,56.494469)
XH=5938.5832
# #自养菌((mgX/L))#
# XA=c(4937.1518,5695.37,4402.0802,3992.0263,2792.156)
XA=4937.1518
# #胞内贮存物浓度#%
# XSTO=c(7.82,5.36,4.11,2.90,1.62)
XSTO=7.82
# #XSTO=(5.5,5.36,4.11,3,2.5)%
X=5.5
#   #溶解性有机氮浓度#
#   SND=c(0.79,0.96,1.12,1.28,1.35)
SND=0.79


COD_1<-function(KO,KS,KUAP,KBAP,KSTO,uH,KNH,KALK,KNO,
                kH,KX,kSTO,kUSTO,kBSTO,nNO,bHO,bHNO,bSTOO,
                bSTONO,uA,KANH,KAO,KANO,KAALK,bAO,bANO,ka,
                YSTOO,YSTONO,YHO,YHNO,YA,iNSS,iNXI,iNBM,iTSBM,
                iTSSTO,iNUAP,iNBAP,kUAPO,kUAPNO,kBAPO,kBAPNO,
                kUAPAO,kBAPAO,kBAPANO,t){
  
  
  x2=YSTOO-1
  x21=YSTOO-1
  x22=YSTOO-1
  x3=(YSTONO-1)/2.86
  x31=(YSTONO-1)/2.86
  x32=(YSTONO-1)/2.86
  x4=1-1/YHO+kUAPO
  x5=(kUAPNO+1-1/YHNO)/2.86
  x6=kBAPO-1
  x7=(kBAPNO-1)/2.86
  x8=-1
  x9=-1/2.86
  x10=kUAPAO+1-4.57/YA
  x11=kBAPAO-1
  x12=(kBAPANO-1)/2.86
  
  y2=iNSS
  y21=iNUAP
  y22=iNBAP
  y3=iNSS
  y31=iNUAP
  y32=iNBAP
  y4=-iNBM-iNUAP*kUAPO
  y5=-iNBM-iNUAP*kUAPNO
  y6=iNBM-iNBAP*kBAPO
  y7=iNBM-iNBAP*kBAPNO
  y8=0
  y9=0
  y10=-1/YA-iNBM-iNUAP*kUAPAO
  y11=iNBM-iNBAP*kBAPAO
  y12=iNBM-iNBAP*kBAPANO
  
  z2=y2/14
  z21=y21/14
  z22=y22/14
  z3=(y3-x3)/14
  z31=(y31-x31)/14
  z32=(y32-x32)/14
  z4=y4/14
  z5=(y5-x5)/14
  z6=y6/14
  z7=(y7-x7)/14
  z8=0
  z9=-x9/14
  z10=(y10-1/YA)/14
  z11=y11/14
  z12=(y12-x12)/14
  
  #%各组分的动力学模型%%%%%
  p2=kSTO*SO*SS*XH/(KO+SO)/(KS+SS)
  p21=kUSTO*SO*SUAP*XH/(KO+SO)/(KUAP+SUAP)
  p22=kBSTO*SO*SBAP*XH/(KO+SO)/(KBAP+SBAP)
  p3=kSTO*nNO*KO*SNO*SS*XH/(KO+SO)/(KNO+SNO)/(KS+SS)
  p31=kUSTO*nNO*KO*SNO*SUAP*XH/(KO+SO)/(KNO+SNO)/(KUAP+SUAP)
  p32=kBSTO*nNO*KO*SNO*SBAP*XH/(KO+SO)/(KNO+SNO)/(KBAP+SBAP)
  p4=uH*SO*SNH*SALK*XSTO/(KO+SO)/(KNH+SNH)/(KALK+SALK)/(KSTO+XSTO/XH)
  p5=uH*nNO*SO*SNO*SNH*SALK*XSTO/(KO+SO)/(KNO+SNO)/(KNH+SNH)/(KALK+SALK)/(KSTO+XSTO/XH)
  p6=bHO*SO*XH/(KO+SO)
  p7=bHNO*KO*SNO*XH/(KO+SO)/(KNO+SNO)
  p8=bSTOO*SO*XSTO/(KO+SO)
  p9=bSTONO*KO*SNO*XSTO/(KO+SO)/(KNO+SNO)
  p10=uA*SO*SNH*SALK*XA/(KAO+SO)/(KANH+SNH)/(KAALK+SALK)
  p11=bAO*SO*XA/(KAO+SO)
  p12=bANO*KO*SNO*XA/(KAO+SO)/(KANO+SNO)
  p13=ka*SND
  
  r1=(YSTOO-1)*p2+(YSTOO-1)*p21+(YSTOO-1)*p22
  +(YSTONO-1)*p3+(YSTONO-1)*p31+(YSTONO-1)*p32
  +(kUAPO+1-1/YHO)*p4+(kUAPNO+1-1/YHNO)*p5+(kBAPO-1)*p6
  +(kBAPNO-1)*p7-p8-p9+(kUAPAO+1)*p10+(kBAPAO-1)*p11+(kBAPANO-1)*p12
  
  R1=r1*t
  SOC=11.3
  SC=SOC+R1
  return(SC)
  
}



T<-seq(from = 0.01, to = 1, by = 0.01)

COD_2 <- function(X, t =T) {
  
  out <- matrix(nrow = nrow(X), ncol = length(t))
  for (i in 1:nrow(X)) {
    out[i, ] <- COD_1(X$KO[i],X$KS[i],X$KUAP[i],X$KBAP[i],
                      X$KSTO[i],X$uH[i],X$KNH[i],X$KALK[i],
                      X$KNO[i],X$kH[i],X$KX[i],X$kSTO[i],X$kUSTO[i],
                      X$kBSTO[i],X$nNO[i],X$bHO[i],X$bHNO[i],X$bSTOO[i],
                      X$bSTONO[i],X$uA[i],X$KANH[i],X$KAO[i],X$KANO[i],
                      X$KAALK[i],X$bAO[i],X$bANO[i],X$ka[i],X$YSTOO[i],
                      X$YSTONO[i],X$YHO[i],X$YHNO[i],X$YA[i],X$iNSS[i],X$iNXI[i],
                      X$iNBM[i],X$iTSBM[i],X$iTSSTO[i],X$iNUAP[i],X$iNBAP[i],
                      X$kUAPO[i],X$kUAPNO[i],X$kBAPO[i],X$kBAPNO[i],X$kUAPAO[i],
                      X$kBAPAO[i],X$kBAPANO[i], t)
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}



#library(boot)
library(sensitivity)
library(multisensi)

COD.seq.sobol <- multisensi(design = sobol, model = COD_2,
                            reduction = NULL, analysis = analysis.sensitivity, 
                            center = F,scale=TRUE,
                            design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                            analysis.args = list(keep.outputs = FALSE))




COD.seq.morris<-multisensi(design = morris,model=COD_2,reduction = NULL,analysis = analysis.sensitivity,
                           #center = TRUE,scale = TRUE,
                           design.args = list(factors=para$parameter,r=4,
                                              center=F,
                                              design=list(type="oat",levels=100,grid.jump=50),
                                              binf=0.01*para$typical_value,
                                              bsup=10*para$typical_value,scale=TRUE))





library(ggplot2)
library(ggrepel)
data_plot_sobol<-as.data.frame(cbind(rownames(COD.seq.sobol$mSI),COD.seq.sobol$mSI$GSI))
names(data_plot_sobol)<-c("parameter","sensitivity")

data_plot_sobol$sensitivity<-as.numeric(data_plot_sobol$sensitivity)

p1<-ggplot(data=data_plot_sobol,aes(x=parameter,y=sensitivity))+
  geom_text_repel(data = data_plot_sobol[which(data_plot_sobol$sensitivity>mean(data_plot_sobol$sensitivity)),],
                  aes(x=parameter,y=sensitivity,label=parameter))+
  geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Sensitivity of COD ")+
  xlab("Based on sobol")+
  theme(text=element_text(size=16,  family="serif"))

p1



data_plot_morris<-as.data.frame(cbind(rownames(COD.seq.morris$mSI),COD.seq.morris$mSI$t1))
names(data_plot_morris)<-c("parameter","sensitivity")

data_plot_morris$sensitivity<-as.numeric(data_plot_morris$sensitivity)

p2<-ggplot(data=data_plot_morris,aes(x=parameter,y=sensitivity))+
  geom_text_repel(data = data_plot_morris[which(data_plot_morris$sensitivity>mean(data_plot_morris$sensitivity)),],
                  aes(x=parameter,y=sensitivity,label=parameter))+
  geom_point()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("Sensitivity of COD ")+
  xlab("Based on morris")+
  theme(text=element_text(size=16,  family="serif"))

p2

library(cowplot)
plot_grid(p1,p2)


save.image(file = "OCD_model.RData")
