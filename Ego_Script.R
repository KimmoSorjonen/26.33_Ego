

############################# BIENVENUE #############################

############# I SOLEMNLY SWEAR THAT I AM UP TO NO GOOD ##############

## Loading package

library(MASS)
library(metafor)
library(lavaan)

#################
## Data
## Correlations reported in Meng et al. (2026)
## order: pts1-3, ed1-3, dc1-3

rm <- matrix(c(
  
   1.00,  0.40,  0.32, -0.27, -0.17, -0.20, -0.56, -0.28, -0.29,
   0.40,  1.00,  0.47, -0.17, -0.28, -0.26, -0.31, -0.51, -0.37,
   0.32,  0.47,  1.00, -0.15, -0.20, -0.22, -0.24, -0.31, -0.60,
  -0.27, -0.17, -0.15,  1.00,  0.34,  0.31,  0.56,  0.28,  0.22,
  -0.17, -0.28, -0.20,  0.34,  1.00,  0.48,  0.26,  0.60,  0.34,
  -0.20, -0.26, -0.22,  0.31,  0.48,  1.00,  0.26,  0.39,  0.58,
  -0.56, -0.31, -0.24,  0.56,  0.26,  0.26,  1.00,  0.37,  0.33,
  -0.28, -0.51, -0.31,  0.28,  0.60,  0.39,  0.37,  1.00,  0.44,
  -0.29, -0.37, -0.60,  0.22,  0.34,  0.58,  0.33,  0.44,  1.00), nrow=9)

n <- 2471  ## sample size

dfall <- data.frame(mvrnorm(n=n, ## generating data frame with required size and corr.
            mu=rep(0,9), Sigma=rm, empirical=T))

df.pe <- dfall[,1:6] ## pts and ed
df.pd <- dfall[,c(1:3,7:9)] ## pts and dc
df.ed <- dfall[,4:9] ## ed and dc

dflist <- list(df.pe,df.pd,df.ed) ## list with all three combinations

## The six models, with dy1 = Y2-Y1 and dy2 = Y3-Y2

m1 <- "dy1 ~ x1 + y1"
m2 <- "dy2 ~ x2 + y2"
m3 <- "dy1 ~ x1 + y2"
m4 <- "dy2 ~ x2 + y3"
m5 <- "dy1 ~ x1"
m6 <- "dy2 ~ x2"

modlist <- list(m1,m2,m3,m4,m5,m6) ## list with all six models

#################
## Figure

panlab <- c("PTS → ED", "ED → PTS",
            "PTS → DC", "DC → PTS",
            "ED → DC", "DC → ED")

slab <- c("7.RMA","",
          "6.y3-y2.no.adj", "5.y2-y1.no.adj", ## y-labels
          "4.y3-y2.adj.y3", "3.y2-y1.adj.y2", 
          "2.y3-y2.adj.y2", "1.y2-y1.adj.y1")

cx <- 0.8 ## sizing factor
r.low <- -0.6 ## range, lower
r.upp <- 0.6 ## range, upper
f.upp <- r.upp+0.8*(r.upp-r.low) ## room for text
tic1 <- 0.3 ## distance, tics
tic2 <- 0.6 ## distance, labels

if(dev.cur()==2) dev.off() ## removing earlier plots
par(mar=c(1,1,1.2,0), oma=c(1.5,4.5,0.5,1), mfrow=c(3,2)) ## setting margins and layout

pl <- 1 ## keeping track of panels

for(k in 1:3){ ## three rows on panels
  
  df <- dflist[[k]] ## pickin the data
  
  for(i in 1:2){ ## two orders
    
    if(i==1) names(df) <- c("x1","x2","x3","y1","y2","y3") ## order 1
    if(i==2) names(df) <- c("y1","y2","y3","x1","x2","x3") ## order 2
    
    dy1 <- df$y2 - df$y1 ## difference score 1
    dy2 <- df$y3 - df$y2 ## difference score 2
    
    numeff <- 8 ## number of rows per panel (one is empty)
    
    plot(c(r.low,f.upp),c(0.5,8.5), type="n",xaxt="n",yaxt="n",xlab="",ylab="") ## empty plot
    
    mtext(panlab[pl],3, line=0.2, cex=cx) ## panel label
    
    axis(1,at=seq(r.low,r.upp,tic1),labels=F, cex.axis=cx) ## x-axis
    if(k==3) axis(1,at=seq(r.low,r.upp,tic2),labels=seq(r.low,r.upp,tic2), cex.axis=cx) ## x-labels
    axis(2,at=1:8,labels=F, las=1, cex.axis=cx) ## y-tics
    if(i==1) axis(2,at=1:8,labels=slab, las=1, cex.axis=cx) ## y-labels
    
    lines(c(0,0),c(-1,15),col="gray") ## vertical gray line at x=0
    
    ball <- vector() ## to be filled with effects below
    seall <- vector() ## to be filled with standard errors below
    
    for(j in 1:6){ ## for the six models
      
      fit <- lm(modlist[[j]], data=df) ## fitting the model
      b <- fit$coefficients[2] ## reg. coefficient
      low <- confint(fit)[2,1] ## CI, low
      upp <- confint(fit)[2,2] ## CI, upp
      se <- summary(fit)$coefficients[2,2] ## standard error
      tx <- paste(round(b,2)," [", round(low,2),"; ", 
                  round(upp,2),"]", sep="") ## string with values
      
      ball <- c(ball,b) ## adding effect to object
      seall <- c(seall,se) ## adding se to object
      
      lines(c(-2,r.upp),c(numeff,numeff), col="gray") ## vertical gray line
      arrows(low,numeff,upp,numeff,angle=90,code=3,length=0.05,lwd=2) ## CI
      points(b,numeff,pch=21,col="black",bg="black") ## point for effect
      text(r.upp,numeff,tx,pos=4,cex=cx) ## adding string with values
      
      numeff <- numeff-1 ## next row, please
    }
    
    ## Meta-analysis of the six effects
    
    b.fish <- 0.5*log((1+ball)/(1-ball)) ## Fisher's transformation of effects
    se.fish <- 0.5*log((1+seall)/(1-seall)) ## Fisher's transformation of SE
    
    ###### Following Bartos et al.
    
    we <- rep(1/6,6) ## equal weight to all effects
    
    random1 <- rma(yi=b.fish, vi=se.fish^2)
    res.ma <- rma(yi=b.fish, vi=(se.fish^2)/we, tau=random1$tau2)
    
    ######
    
    pred <- predict(res.ma, transf=transf.ztor) ## transforms back from Fisher's
    metb <- pred$pred ## estimate
    metlow <- pred$ci.lb ## lower CI
    metupp <- pred$ci.ub ## upper CI
    
    ##
    
    lines(c(-2,4),c(2,2),lty=2) ## dashed line
    lines(c(-2,r.upp),c(1,1),col="gray") ## gray vertical line
    
    polygon(c(metlow,metb,metupp,metb), c(1,1.5,1,0.5), ## diamond for RMA 
            col = "black", border = "black", lwd = 1)
    
    tx.ma <- paste(round(metb,2)," [", round(metlow,2),"; ", 
                   round(metupp,2),"]", sep="") ## string with values
    text(r.upp,1,tx.ma,pos=4,cex=cx) ## adding string of values
    legend("topleft", title=LETTERS[pl], legend="", bty="n", inset=0, cex=1.7*cx) ## A-B legend
    
    pl <- pl+1
  }
}

####################
## Alternative model

mosla <- "

## Trait

tx =~ 1*x1+1*x2+1*x3
ty =~ 1*y1+1*y2+1*y3

tx ~~ ty

## State

st1 =~ sx*x1+sy*y1
st2 =~ sx*x2+sy*y2
st3 =~ sx*x3+sy*y3

st2 ~ as*st1
st3 ~ as*st2

## (Error) variances

x1 ~~ x1
x2 ~~ x2
x3 ~~ x3

y1 ~~ y1
y2 ~~ y2
y3 ~~ y3

tx ~~ tx
ty ~~ ty

st1 ~~ 1*st1
st2 ~~ 1*st2
st3 ~~ 1*st3

## Intercepts

x1 ~ 0*1
x2 ~ 0*1
x3 ~ 0*1

y1 ~ 0*1
y2 ~ 0*1
y3 ~ 0*1

tx ~ 0*1
ty ~ 0*1

st1 ~ 0*1
st2 ~ 0*1
st3 ~ 0*1

"

df <- df.pd
names(df) <- c("x1","x2","x3","y1","y2","y3") ## pts as X and dc as Y
fit <- lavaan(mosla, data=df) ## fitting mosla
summary(fit, fit.measures=T, standardized=T, ci=T) ## let's have a look


########################## MISCHIEF MANAGED #########################

############################# AU REVOIR #############################

