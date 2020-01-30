library(viridis)
palette(viridis(256))
# a function we will work with
F1<-function(x){
		return(c(x^2,2*x))
}

F2<-function(x){
		return(c(sin(x),cos(x)))
}

F3<-function(x){
	y<-(x-2)^3-6*x
	deriv<-3*(x-2)^2-6
	return(c(y,deriv))
}

F4<-function(x){
	return(c(cos(x)-x,-1*sin(x)-1))
}

F5<-function(x){
  return(c(x^3+3*x^2-5*x-10,3*x^2+6*x-5))
}

F6<-function(z){ # this has three roots 1,-1/2 + sqrt(3)/2i and -1/2 - sqrt(3)/2i
  # a complex function z^3-1
  zcubed<-ComplexPower(z,3)
  zsquared<-ComplexPower(z,2)

  return(c(zcubed[1]-1,zcubed[2],3*zsquared))
}


F7<-function(z){
  # a complex function z^4-1
  y<-ComplexPower(z,4)-c(1,0)
  deriv<-4*ComplexPower(z,3)

  return(c(y[1],y[2],deriv))
}

F8<-function(z){
   # a complex function 35z^9 - 180z^7 + 378z^5 - 420z^3 + 315z    (taken from http://www.chiark.greenend.org.uk/~sgtatham/newton/)
   y<-(35*ComplexPower(z,9)-180*ComplexPower(z,7)+378*ComplexPower(z,5)-420*ComplexPower(z,3)+315*z)/315
   deriv<-(35*9*ComplexPower(z,8)-180*7*ComplexPower(z,6)+378*5*ComplexPower(z,4)-420*3*ComplexPower(z,2)+c(315,0))/315

   return(c(y[1],y[2],deriv))
}


F9<-function(z){
  # a complex function z^3+z^2-z-1   [(z^2-1)(z+1)]
  z<-ComplexPower(z,3)+ComplexPower(z,2)-z+c(1,0)
  deriv<-3*ComplexPower(z,2)+2*z-c(1,0)
  return(c(z[1],z[2],deriv))
}

F10<-function(z){
  # a complex function z^3-z-1 -- look at it from [-1,3] on the x-axis and [-3,3] on the y-axis
  z<-ComplexPower(z,3)-z-c(1,0)
  deriv<-3*ComplexPower(z,2)-c(1,0)
  return(c(z[1],z[2],deriv))
}

F11<-function(z){  # Fatou's function
  # a complex function (z+z^2)/2 --
  z<-(ComplexPower(z,2)+z)/2.0
  deriv<-z+0.5
  return(c(z[1],z[2],deriv))
}


# This is a global variable to control whether we record roots or the number of iterations taken to find them
# This is an inelegant way of doing it. It should really be part of the argument to the relevent functions.
bRootOrIterations<-0   # set <-1 to record roots, or <- 0 to record number of iterations needed


ComplexPower<-function(z,n){
   NewZ<-z
   i<-1
   while (i<n)
   {
     NewZ<-ComplexMult(NewZ,z)
     i <- i+1
   }
   return(NewZ)
 }

ComplexMult<-function(z1,z2){
  # multiplies Z1 by z2 (both complex)
  x<-z1[1]*z2[1]-z1[2]*z2[2]
  y<-z1[1]*z2[2]+z1[2]*z2[1]
  return(c(x,y))
}

ComplexDiv<-function(z1,z2){
  # divides z1 by z2
  Conjugate<-c(z2[1],-1*z2[2])  #form conjugate of z2
  Divisor<-ComplexMult(z2,Conjugate)
  Numerator<-ComplexMult(z1,Conjugate)
  Numerator<-Numerator/Divisor[1]
  return(Numerator)
}


# function declaration
TwoDNewtonRaphson<-function(func,StartingValue,Tolerance,MaxNumberOfIterations){
  i<-0
  a<-1    # in general, a should be =1
  X1 <- StartingValue[1]
  X2 <- StartingValue[2]
  Y1 <- func(StartingValue)[1]
  Y2 <- func(StartingValue)[2]
  NewZ<- StartingValue
  Deviation=Y1*Y1+Y2*Y2
  #segments(X,0,X,Y,lty=2,lwd=2,col="blue")
  #text(X-0.2,-1.2,"X1",cex=0.75)
  #Set up a while loop until we hit the required target accuracy or the max. number of steps
  while ((i<MaxNumberOfIterations)&&(Deviation>Tolerance))
  {
    # Find the next (X1,X2)-value using Newton-Raphson's formula
    Z<-func(c(X1,X2))
    #Z(n+1) <- Z(n)-f(Z(n))/f'(z(n))
    Temp<-ComplexDiv(c(Z[1],Z[2]),c(Z[3],Z[4]))
    NewZ<-c(X1-Temp[1],X2-Temp[2]*c(a,0))
    #cat("Z",i,Deviation,NewZ)

    if ((NewZ[1]=="NaN")||(NewZ[2]=="NaN")){
      cat("Function or derivative not defined error (or deriv is zero.\n")
      cat("\n",NewZ,Z)
      break
    }
    # annotate
    #segments(X,Y,NewX,0,lty=2,lwd=2,col="blue")
    #segments(NewX,0,NewX,func(NewX)[1],lty=2,lwd=2,col="blue")
    #text(NewX-0.2,-1.2,paste("X",i+2,sep=""),cex=0.75)

    # calculate accuracy<- |f(x)-0|
    NewVal <- func(NewZ)
    Deviation <- NewVal[1]*NewVal[1]+NewVal[2]*NewVal[2]
    X1<-NewZ[1]
    X2<-NewZ[2]

    i<-i+1
    #cat(paste("\nIteration ",i,":   X=",NewZ[1],"  Y=",NewZ[2],"  Z=",sqrt(Deviation)))
  }

  # output the result
  if (Deviation<Tolerance){
 #   cat(paste("\nFound the root point: (",NewZ[1],",",NewZ[2],") after ", i, "iterations"))
  }else{
    cat(paste("\nConvergence failure. Deviation:",Deviation, "after ", i, 	"iterations"))}
  if (bRootOrIterations==1){
   return(NewZ)
  }else{
    return(c(i,i))
  }
}


CloseTo<-function(x,y){
  # returns 1 if x is close to y
  if ((x-y)*(x-y)<0.01) {
    return(1)
  }else{
    return (0)
  }
}

RootPlotter<-function(Funcn,xmin,xmax,xsteps,ymin,ymax,ysteps,PtSize)
{
  #plot(0,0,col='white',type="p",xlim=c(-5,5),ylim=c(-5,5))
  #par(new=T)
  #x<-seq(-5.01,4.99,by=1)
  #y<-seq(-5.01,4.99,1)
  x<-seq(xmin,xmax,length.out=xsteps)
  y<-seq(ymin,ymax,length.out=ysteps)
  ThingsToPlot<-c(-9,-9,'white',-9,-9)
  colors<-mat.or.vec(xsteps,ysteps)
  for (i in 1:xsteps){
    for (j in 1:ysteps){
      ThisX<-x[i]
      ThisY<-y[j]
      colors[i,j]<-'black'
      #cat(i,j,sep=" ")
      #cat("\n")
      # The following looks for the roots for z^3-1. You will need to modifyu it to find roots for other functions.
      # the three roots are 1,  -1/2 + sqrt(3)/2i and -1/2 - sqrt(3)/2i
      # find the root for this start point
      #cat("\nFinding root at :",ThisX,ThisY)
      #Root<-TwoDNewtonRaphson(Funcn,c(ThisX,ThisY),1e-6,100)
      Root<-TwoDNewtonRaphson(Funcn,c(ThisX,ThisY),1e-1,200)
      if (bRootOrIterations==0){ # set <-1 to record roots, or <- 0 to record number of iterations needed
        ThisColor<-261+5*Root[1]
        #ThisColor<-261+Root[1]
        colors[i,j]<-ThisColor}else{
          if (CloseTo(Root[1],1)){
            colors[i,j]<-'blue'
          }
          if ((CloseTo(Root[1],-0.5))&&(CloseTo(Root[2],0.8660254))){
            colors[i,j]<-'red'
          }
          if ((CloseTo(Root[1],-0.5))&&(CloseTo(Root[2],-0.8660254))){
            colors[i,j]<-'yellow'
          }
        }
       # if (CloseTo(Root[1],-0.5)){
      #  if (CloseTo(Root[2],-1.366025)){  # this one is -1/2 + sqrt(3)/2i
       #   colors[i,j]<-'red'
      #  }else{   # so this one must be -1/2 - sqrt(3)/2i
      #    if (CloseTo(Root[2],0.366025))  # just to double-check
      #      colors[i,j]<-'yellow'
      #  }
      #}
      ThingsToPlot<-rbind(ThingsToPlot,c(ThisX,ThisY,colors[i,j],Root[1],Root[2]))
      #plot(x,y,col=Mycolor,axes=F)
      #par(new=F)
    }
  }
  plot(ThingsToPlot[,1:2],col=ThingsToPlot[,3],,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16,cex=PtSize)
  return(ThingsToPlot)
}



# The following are a bunch of examples I used to create some of the figures shown in class

#A<-RootPlotter(F8,-1,1,300,-1,1,300,0.2)
#A<-RootPlotter(F8,0.5,1.5,300,-0.5,0.5,300,0.2)
A<-RootPlotter(F8,0.9,1.5,200,-0.2,0.2,200,0.4)
#A<-RootPlotter(F8,1.15,1.35,300,-0.08,0.08,300,0.4)
#A<-RootPlotter(F8,1.27,1.3,300,-0.02,0.02,300,0.4)
#A<-RootPlotter(F8,1.283,1.287,300,-0.003,0.003,300,0.4)
#A<-RootPlotter(F8,1.285,1.2866,300,0,0.001,300,0.4)

#RootPlotter(F11,-5,5,300,-5,5,300,0.2)
#RootPlotter(F11,-0.1,0.25,300,-1.25,-0.25,300,0.2)
# RootPlotter(F11,-0.2,0.1,300,-1.01,-0.4,300,0.2)
#RootPlotter(F11,-0.02,0.05,300,-1.0,-0.8,300,0.2)
#RootPlotter(F11,-0.01,0,300,-1.0,-0.95,300,0.2)
#RootPlotter(F11,-0.006,-0.00575,300,-0.968,-0.967,300,0.2)

