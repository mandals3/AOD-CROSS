### Title: Non-parametric estimation of the age-at-onset distribution from a cross-sectional sample
### Authors: Soutrik Mandal, Jing Qin, Ruth M. Pfeiffer
### Journal: Biometrics
### Date: July 25, 2022

### This code can be used to implement the proposed EM algorithm.
### We provide an example on how to use this code at the end.
### Further details are given in the READ_ME.pdf file.

###################################################
### negative observed data log likelihood function
### tlik is used in parameter estimation using the EM method (unconstrained)
###################################################

   tlik=function(o,pc,na01,nt01,na00,na11,nt11,na10,
                 ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
   {

     K1= length(pc$t1cuts)+1
     K2= pc$K2

     p0=o[1:(K1*K2)]
     p1=o[(K1*K2+1):(2*K1*K2)]
     theta=o[2*K1*K2+1]
     p0=matrix(p0,K1,K2,byrow=T)
     p1=matrix(p1,K1,K2,byrow=T)

### g=0 group calculation

     lik0=0
     for(i1 in 1:na01)
     {
       for(j1 in 1:nt01)
       {
         w1=matrix(0,K1,K2)
         w2=w1
         w3=w1

         Delta0=matrix(0,K1,K2)
         Delta1=Delta0

         aa=ua01[i1]
         bb=ut01[j1]
         cc=tab01[j1,i1]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j==bb & bb<=aa & j+k>aa) 
             w2[j,k]=1

             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) w3[j,k]=1
           }
         }

         q2=sum(p0*w2*theta)
         r2=0

         if(q2>0) r2=cc*log(q2)
         lik0=lik0+r2

         Del0= Delta0*p0
         Del1= Delta1*p1

         Del=theta*Del0+(1-theta)*Del1
         Del=sum(Del)

         te=0
         if(Del>0) te=cc*log(Del)
         lik0=lik0-te

       }
     }

     for(i in 1:na00)
     {
       w1=matrix(0,K1,K2)
       Delta0=matrix(0,K1,K2)
       Delta1=Delta0
       ww3=w1

       aa=ua00[i]
       cc=tab00[i]

       for(j in 1:K1)
       {
         for(k in 1:K2)
         {
           if(j>aa ) w1[j,k]=1
           if(j+k>aa) Delta0[j,k]=1
           if(j+k>aa) Delta1[j,k]=1
           if(j+k<=aa) ww3[j,k]=1
         }
       }

       q1=sum(p0*w1*theta)
       if(q1>0) q1=cc*log(q1)
       lik0=lik0+q1

       r1=theta*p0*Delta0+(1-theta)*p1*Delta1
       if(sum(r1)>0) lik0=lik0-cc*log(sum(r1))
     }

### g==1 group calculation

     for(i1 in 1:na11)
     {
       for(j1 in 1:nt11)
       {
         w1=matrix(0,K1,K2)
         w2=w1
         w3=w1

         Delta0=matrix(0,K1,K2)
         Delta1=Delta0

         aa=ua11[i1]
         bb=ut11[j1]
         cc=tab11[j1,i1]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j==bb & bb<=aa & j+k>aa) 
             w2[j,k]=1

             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) w3[j,k]=1
           }
         }

         q2=sum(p1*w2*(1-theta))
         if(q2>0) q2=cc*log(q2)
      
         lik0=lik0+q2
         Del0= Delta0*p0
         Del1= Delta1*p1
         Del=theta*Del0+(1-theta)*Del1
         Del=sum(Del)

         if(Del>0) lik0=lik0-cc*log(Del)
       }
     }

     for(i in 1:na10)
     {
       w1=matrix(0,K1,K2)
       Delta0=matrix(0,K1,K2)
       Delta1=Delta0
       ww3=w1

       aa=ua10[i]
       cc=tab10[i]

       for(j in 1:K1)
       {
         for(k in 1:K2)
         {
           if(j>aa ) w1[j,k]=1
           if(j+k>aa) Delta0[j,k]=1
           if(j+k>aa) Delta1[j,k]=1
           if(j+k<=aa) ww3[j,k]=1
         }
       }

       q1=sum(p1*w1*(1-theta))
       if(q1>0) lik0=lik0+cc*log(q1)

       r1=theta*p0*Delta0+(1-theta)*p1*Delta1
       if(sum(r1)>0) lik0=lik0-cc*log(sum(r1))
     }

     return(-lik0)
   }


##############################################################
### function used for confidence interval of theta (genotype prevalence):
### calculates negative observed data log likelihood at
### p0 and p1 estimated via EM and theta fixed at \tilde theta
##############################################################

cpthetafn= function(thetatilde,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
{
   K1= length(pc$t1cuts)+1
   K2= pc$K2

   np0= matrix(0,K1,K2)
   np1= matrix(0,K1,K2)

### EM starts

   tt=1000000000

   for(ii in 1:20)
   {
     tol=1
     count=0

     p0=matrix(runif(K1*K2),K1,K2)
     p0=p0/sum(p0)
     p1=matrix(runif(K1*K2),K1,K2)
     p1=p1/sum(p1)
     theta=thetatilde

     tm1=0

     while(tol>0.0001 & count<=500)
     {
       count=count+1

       f01=matrix(0,K1,K2)
       f02=f01
       f03=f01

       f11=f01
       f12=f01
       f13=f01

### g=0 group calculation

       for(i1 in 1:na01)
       {
         for(j1 in 1:nt01)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua01[i1]
           bb=ut01[j1]
           cc=tab01[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1
               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p0*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)
  
           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p0*w3
           if(sum(Del)>0) q3=theta*cc*q3/sum(Del)

           f02=f02+q2
           f03=f03+q3

           q31=p1*w3
           if(sum(Del)>0) q31=(1-theta)*cc*q31/sum(Del)

           f13=f13+q31
         }
       }

       for(i in 1:na00)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua00[i]
         cc=tab00[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p0*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p0
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= theta*cc*q3/sum(r1)

         f01=f01+q1
         f03=f03+q3

         q31=ww3*p1
         if(sum(r1)>0) q31= (1-theta)*cc*q31/sum(r1)

         f13=f13+q31
       }

### g==1 group calculation

       for(i1 in 1:na11)
       {
         for(j1 in 1:nt11)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua11[i1]
           bb=ut11[j1]
           cc=tab11[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1

               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p1*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)

           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p1*w3
           if(sum(Del)>0) q3=(1-theta)*cc*q3/sum(Del)
 
           f12=f12+q2
           f13=f13+q3

           q30=p0*w3
           if(sum(Del)>0) q30=theta*cc*q30/sum(Del)
           f03=f03+q30
         }
       }

       for(i in 1:na10)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua10[i]
         cc=tab10[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p1*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p1
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= (1- theta)*cc*q3/sum(r1)

         f11=f11+q1
         f13=f13+q3

         q30=ww3*p0
         if(sum(r1)>0) q30= theta*cc*q30/sum(r1)

         f03=f03+q30
       }

       f0=f01+f02+f03
       f1=f11+f12+f13

       np0=f0/sum(f0)
       np1=f1/sum(f1)
       ntheta= theta

       tol=max(abs(p0-np0))+max(abs(p1-np1))+abs(theta-ntheta)
       p1=np1
       p0=np0
       theta=ntheta

      o1= c(c(t(np0)),c(t(np1)),ntheta)
 
       lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
       tm1=c(tm1,lik)
     }

      o1= c(c(t(np0)),c(t(np1)),ntheta)

     lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

     if(lik<tt)
     {
       so0=np0
       so1=np1
       tt=lik
     }

   }

   np0=t(so0)
   np1=t(so1)

   p0= colSums(np0)
   p1= colSums(np1)
   outg1= ntheta

  so1= c(c(t(so0)),c(t(so1)),outg1)
  alike0= tlik(so1,pc,na01,nt01,na00,na11,nt11,na10,
              ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

  return(alike0)
}



################################################
### function used for confidence interval of p0:
################################################

cpp0fn= function(pfix0,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
{
   K1= length(pc$t1cuts)+1
   K2= pc$K2

   np0= matrix(0,K1,K2)
   np1= matrix(0,K1,K2)

### EM starts

   tt=1000000000

   for(ii in 1:20)
   {
     tol=1
     count=0

     p0=matrix(runif(K1*K2),K1,K2)
     p0=p0/sum(p0)
     p1=matrix(runif(K1*K2),K1,K2)
     p1=p1/sum(p1)

     theta=runif(1)
     tm1=0

     while(tol>0.0001 & count<=500)
     {
       count=count+1

       f01=matrix(0,K1,K2)
       f02=f01
       f03=f01

       f11=f01
       f12=f01
       f13=f01

### g=0 group calculation

       for(i1 in 1:na01)
       {
         for(j1 in 1:nt01)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua01[i1]
           bb=ut01[j1]
           cc=tab01[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1
               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p0*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)
  
           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p0*w3
           if(sum(Del)>0) q3=theta*cc*q3/sum(Del)

           f02=f02+q2
           f03=f03+q3

           q31=p1*w3
           if(sum(Del)>0) q31=(1-theta)*cc*q31/sum(Del)

           f13=f13+q31
         }
       }

       for(i in 1:na00)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua00[i]
         cc=tab00[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p0*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p0
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= theta*cc*q3/sum(r1)

         f01=f01+q1
         f03=f03+q3

         q31=ww3*p1
         if(sum(r1)>0) q31= (1-theta)*cc*q31/sum(r1)

         f13=f13+q31
       }

### g==1 group calculation

       for(i1 in 1:na11)
       {
         for(j1 in 1:nt11)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua11[i1]
           bb=ut11[j1]
           cc=tab11[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1

               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p1*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)

           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p1*w3
           if(sum(Del)>0) q3=(1-theta)*cc*q3/sum(Del)
 
           f12=f12+q2
           f13=f13+q3

           q30=p0*w3
           if(sum(Del)>0) q30=theta*cc*q30/sum(Del)
           f03=f03+q30
         }
       }

       for(i in 1:na10)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua10[i]
         cc=tab10[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p1*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p1
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= (1- theta)*cc*q3/sum(r1)

         f11=f11+q1
         f13=f13+q3

         q30=ww3*p0
         if(sum(r1)>0) q30= theta*cc*q30/sum(r1)

         f03=f03+q30
       }

       f0=f01+f02+f03
       f1=f11+f12+f13

### p0 estimate constrained

       nu0= (sum(f0)-sum(f0[jj,]))/(1-pfix0)
       lam0= sum(f0[jj,])/pfix0 - nu0

       np0[jj,]= f0[jj,]/(lam0+nu0)
       np0[-jj,]= f0[-jj,]/nu0

### p1 estimate unconstrained

       np1=f1/sum(f1)

### theta estimate

       ntheta=sum(f0)/(sum(f0)+sum(f1))

       tol=max(abs(p0-np0))+max(abs(p1-np1))+abs(theta-ntheta)
       p1=np1
       p0=np0
       theta=ntheta

      o1= c(c(t(np0)),c(t(np1)),ntheta)
 
       lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
       tm1=c(tm1,lik)
     }

      o1= c(c(t(np0)),c(t(np1)),ntheta)

     lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

     if(lik<tt)
     {
       so0=np0
       so1=np1
       tt=lik
     }

   }

   np0=t(so0)
   np1=t(so1)

   out= list(np0, np1, ntheta)
   return(out)

}



#################################################
### function used for confidence interval of p1:
#################################################

cpp1fn= function(pfix1,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
{
   K1= length(pc$t1cuts)+1
   K2= pc$K2

   np0= matrix(0,K1,K2)
   np1= matrix(0,K1,K2)

### EM starts

   tt=1000000000

   for(ii in 1:20)
   {
     tol=1
     count=0

     p0=matrix(runif(K1*K2),K1,K2)
     p0=p0/sum(p0)
     p1=matrix(runif(K1*K2),K1,K2)
     p1=p1/sum(p1)

     theta=runif(1)
     tm1=0

     while(tol>0.0001 & count<=500)
     {
       count=count+1

       f01=matrix(0,K1,K2)
       f02=f01
       f03=f01

       f11=f01
       f12=f01
       f13=f01

### g=0 group calculation

       for(i1 in 1:na01)
       {
         for(j1 in 1:nt01)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua01[i1]
           bb=ut01[j1]
           cc=tab01[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1
               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p0*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)
  
           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p0*w3
           if(sum(Del)>0) q3=theta*cc*q3/sum(Del)

           f02=f02+q2
           f03=f03+q3

           q31=p1*w3
           if(sum(Del)>0) q31=(1-theta)*cc*q31/sum(Del)

           f13=f13+q31
         }
       }

       for(i in 1:na00)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua00[i]
         cc=tab00[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p0*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p0
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= theta*cc*q3/sum(r1)

         f01=f01+q1
         f03=f03+q3

         q31=ww3*p1
         if(sum(r1)>0) q31= (1-theta)*cc*q31/sum(r1)

         f13=f13+q31
       }

### g==1 group calculation

       for(i1 in 1:na11)
       {
         for(j1 in 1:nt11)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua11[i1]
           bb=ut11[j1]
           cc=tab11[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1

               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p1*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)

           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p1*w3
           if(sum(Del)>0) q3=(1-theta)*cc*q3/sum(Del)
 
           f12=f12+q2
           f13=f13+q3

           q30=p0*w3
           if(sum(Del)>0) q30=theta*cc*q30/sum(Del)
           f03=f03+q30
         }
       }

       for(i in 1:na10)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua10[i]
         cc=tab10[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p1*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p1
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= (1- theta)*cc*q3/sum(r1)

         f11=f11+q1
         f13=f13+q3

         q30=ww3*p0
         if(sum(r1)>0) q30= theta*cc*q30/sum(r1)

         f03=f03+q30
       }

       f0=f01+f02+f03
       f1=f11+f12+f13

### p0 estimate unconstrained

       np0=f0/sum(f0)

### p1 estimate constrained

       nu1= (sum(f1)-sum(f1[jj,]))/(1-pfix1)
       lam1= sum(f1[jj,])/pfix1 - nu1

       np1[jj,]= f1[jj,]/(lam1+nu1)
       np1[-jj,]= f1[-jj,]/nu1

### theta estimate

       ntheta=sum(f0)/(sum(f0)+sum(f1))

       tol=max(abs(p0-np0))+max(abs(p1-np1))+abs(theta-ntheta)
       p1=np1
       p0=np0
       theta=ntheta

       o1= c(c(t(np0)),c(t(np1)),ntheta)
 
       lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
       tm1=c(tm1,lik)
     }

     o1= c(c(t(np0)),c(t(np1)),ntheta)

     lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

     if(lik<tt)
     {
       so0=np0
       so1=np1
       tt=lik
     }

   }

   np0=t(so0)
   np1=t(so1)

   out= list(np0, np1, ntheta)
   return(out)

}

###########################################################
### function used for confidence interval of p0 cumulative:
###########################################################

cpp0fn_cumsum= function(pfix0,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
{

   K1= length(pc$t1cuts)+1
   K2= pc$K2

   np0= matrix(0,K1,K2)
   np1= matrix(0,K1,K2)

### EM starts

   tt=1000000000

   for(ii in 1:20)
   {
     tol=1
     count=0

     p0=matrix(runif(K1*K2),K1,K2)
     p0=p0/sum(p0)
     p1=matrix(runif(K1*K2),K1,K2)
     p1=p1/sum(p1)

     theta=runif(1)
     tm1=0

     while(tol>0.0001 & count<=500)
     {
       count=count+1

       f01=matrix(0,K1,K2)
       f02=f01
       f03=f01

       f11=f01
       f12=f01
       f13=f01

### g=0 group calculation

       for(i1 in 1:na01)
       {
         for(j1 in 1:nt01)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua01[i1]
           bb=ut01[j1]
           cc=tab01[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1
               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p0*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)
  
           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p0*w3
           if(sum(Del)>0) q3=theta*cc*q3/sum(Del)

           f02=f02+q2
           f03=f03+q3

           q31=p1*w3
           if(sum(Del)>0) q31=(1-theta)*cc*q31/sum(Del)

           f13=f13+q31
         }
       }

       for(i in 1:na00)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua00[i]
         cc=tab00[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p0*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p0
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= theta*cc*q3/sum(r1)

         f01=f01+q1
         f03=f03+q3

         q31=ww3*p1
         if(sum(r1)>0) q31= (1-theta)*cc*q31/sum(r1)

         f13=f13+q31
       }

### g==1 group calculation

       for(i1 in 1:na11)
       {
         for(j1 in 1:nt11)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua11[i1]
           bb=ut11[j1]
           cc=tab11[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1

               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p1*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)

           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p1*w3
           if(sum(Del)>0) q3=(1-theta)*cc*q3/sum(Del)
 
           f12=f12+q2
           f13=f13+q3

           q30=p0*w3
           if(sum(Del)>0) q30=theta*cc*q30/sum(Del)
           f03=f03+q30
         }
       }

       for(i in 1:na10)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua10[i]
         cc=tab10[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p1*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p1
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= (1- theta)*cc*q3/sum(r1)

         f11=f11+q1
         f13=f13+q3

         q30=ww3*p0
         if(sum(r1)>0) q30= theta*cc*q30/sum(r1)

         f03=f03+q30
       }

       f0=f01+f02+f03
       f1=f11+f12+f13

### p0 estimate constrained

       nu0= (sum(f0)-sum(f0[1:jj,]))/(1-pfix0)
       lam0= sum(f0[1:jj,])/pfix0 - nu0

       np0[1:jj,]= f0[1:jj,]/(lam0+nu0)
       np0[(jj+1):K1,]= f0[(jj+1):K1,]/nu0

### p1 estimate unconstrained

       np1=f1/sum(f1)

### theta estimate

       ntheta=sum(f0)/(sum(f0)+sum(f1))


       tol=max(abs(p0-np0))+max(abs(p1-np1))+abs(theta-ntheta)
       p1=np1
       p0=np0
       theta=ntheta

       o1= c(as.vector(t(np0)),as.vector(t(np1)),ntheta)
 
       lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
       tm1=c(tm1,lik)

     }

     o1= c(as.vector(t(np0)),as.vector(t(np1)),ntheta)

     lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

     if(lik<tt)
     {
       so0=np0
       so1=np1
       tt=lik
     }

   }

   np0=t(so0)
   np1=t(so1)

   out= list(np0, np1, ntheta)
   return(out)

}


###########################################################
### function used for confidence interval of p1 cumulative:
###########################################################

cpp1fn_cumsum= function(pfix1,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
{
   K1= length(pc$t1cuts)+1
   K2= pc$K2

   np0= matrix(0,K1,K2)
   np1= matrix(0,K1,K2)

### EM starts

   tt=1000000000

   for(ii in 1:20)
   {
     tol=1
     count=0

     p0=matrix(runif(K1*K2),K1,K2)
     p0=p0/sum(p0)
     p1=matrix(runif(K1*K2),K1,K2)
     p1=p1/sum(p1)

     theta=runif(1)
     tm1=0

     while(tol>0.0001 & count<=500)
     {
       count=count+1

       f01=matrix(0,K1,K2)
       f02=f01
       f03=f01

       f11=f01
       f12=f01
       f13=f01

### g=0 group calculation

       for(i1 in 1:na01)
       {
         for(j1 in 1:nt01)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua01[i1]
           bb=ut01[j1]
           cc=tab01[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1
               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p0*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)
  
           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p0*w3
           if(sum(Del)>0) q3=theta*cc*q3/sum(Del)

           f02=f02+q2
           f03=f03+q3

           q31=p1*w3
           if(sum(Del)>0) q31=(1-theta)*cc*q31/sum(Del)

           f13=f13+q31
         }
       }

       for(i in 1:na00)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua00[i]
         cc=tab00[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p0*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p0
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= theta*cc*q3/sum(r1)

         f01=f01+q1
         f03=f03+q3

         q31=ww3*p1
         if(sum(r1)>0) q31= (1-theta)*cc*q31/sum(r1)

         f13=f13+q31
       }

### g==1 group calculation

       for(i1 in 1:na11)
       {
         for(j1 in 1:nt11)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua11[i1]
           bb=ut11[j1]
           cc=tab11[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1

               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p1*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)

           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p1*w3
           if(sum(Del)>0) q3=(1-theta)*cc*q3/sum(Del)
 
           f12=f12+q2
           f13=f13+q3

           q30=p0*w3
           if(sum(Del)>0) q30=theta*cc*q30/sum(Del)
           f03=f03+q30
         }
       }

       for(i in 1:na10)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua10[i]
         cc=tab10[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p1*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p1
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= (1- theta)*cc*q3/sum(r1)

         f11=f11+q1
         f13=f13+q3

         q30=ww3*p0
         if(sum(r1)>0) q30= theta*cc*q30/sum(r1)

         f03=f03+q30
       }

       f0=f01+f02+f03
       f1=f11+f12+f13

### p0 estimate unconstrained

       np0=f0/sum(f0)

### p1 estimate constrained

       nu1= (sum(f1)-sum(f1[1:jj,]))/(1-pfix1)
       lam1= sum(f1[1:jj,])/pfix1 - nu1

       np1[1:jj,]= f1[1:jj,]/(lam1+nu1)
       np1[(jj+1):K1,]= f1[(jj+1):K1,]/nu1

### theta estimate

       ntheta=sum(f0)/(sum(f0)+sum(f1))

       tol=max(abs(p0-np0))+max(abs(p1-np1))+abs(theta-ntheta)
       p1=np1
       p0=np0
       theta=ntheta

       o1= c(as.vector(t(np0)),as.vector(t(np1)),ntheta)
 
       lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
       tm1=c(tm1,lik)

     }

     o1= c(as.vector(t(np0)),as.vector(t(np1)),ntheta)

     lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

     if(lik<tt)
     {
       so0=np0
       so1=np1
       tt=lik
     }

   }

   np0=t(so0)
   np1=t(so1)

   out= list(np0, np1, ntheta)
   return(out)

}


#######################
###    EM code
#######################

func.em= function(pc,input.data)
{

### original data

   tt1=input.data[,pc$ageatdx]
   A1=input.data[,pc$ageatstudy]
   g=input.data[,pc$mutation]

   n= length(g)
   K1= length(pc$t1cuts)+1
   K2= pc$K2

   tt1[is.na(tt1)]=A1[is.na(tt1)]+1000

   A1= findInterval(A1,pc$acuts)
   tt1= findInterval(tt1,pc$t1cuts)

   A1= A1+1
   tt1= tt1+1

   t= tt1

### divide into g=0 or 1 groups
   A=A1

   t0=t[g==0]
   A0=A[g==0]

   t1=t[g==1]
   A1=A[g==1]


### tabulate data for censored and uncensored T1 separately
### in g=0 group

   m=length(A0)

   d=rep(0,m)
   d[t0<=A0]=1

   a00=A0[d==0]

   tab00=table(a00)
   ua00=sort(unique(a00))
   na00=length(ua00)

   a01=A0[d==1]
   t01=t0[d==1]

   tab01=table(t01,a01)
   ua01=sort(unique(a01))
   na01=length(ua01)

   ut01=sort(unique(t01))
   nt01=length(ut01)


### tabulate data for censored and uncensored T1 separately
### in g=1 group

   m=length(A1)

   d=rep(0,m)
   d[t1<=A1]=1

   a10=A1[d==0]

   tab10=table(a10)
   ua10=sort(unique(a10))
   na10=length(ua10)

   a11=A1[d==1]
   t11=t1[d==1]

   tab11=table(t11,a11)
   ua11=sort(unique(a11))
   na11=length(ua11)

   ut11=sort(unique(t11))
   nt11=length(ut11)

### EM starts

   set.seed(1)

   tt=1000000000

   for(ii in 1:20)
   {
     tol=1
     count=0

     p0=matrix(runif(K1*K2),K1,K2)
     p0=p0/sum(p0)
     p1=matrix(runif(K1*K2),K1,K2)
     p1=p1/sum(p1)

     theta=runif(1)
     tm1=0

     while(tol>0.0001 & count<=500)
     {
       count=count+1

       f01=matrix(0,K1,K2)
       f02=f01
       f03=f01

       f11=f01
       f12=f01
       f13=f01

### g=0 group calculation

       for(i1 in 1:na01)
       {
         for(j1 in 1:nt01)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua01[i1]
           bb=ut01[j1]
           cc=tab01[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1
               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p0*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)
  
           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p0*w3
           if(sum(Del)>0) q3=theta*cc*q3/sum(Del)

           f02=f02+q2
           f03=f03+q3

           q31=p1*w3
           if(sum(Del)>0) q31=(1-theta)*cc*q31/sum(Del)

           f13=f13+q31
         }
       }

       for(i in 1:na00)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua00[i]
         cc=tab00[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p0*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p0
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= theta*cc*q3/sum(r1)

         f01=f01+q1
         f03=f03+q3

         q31=ww3*p1
         if(sum(r1)>0) q31= (1-theta)*cc*q31/sum(r1)

         f13=f13+q31
       }

### g==1 group calculation

       for(i1 in 1:na11)
       {
         for(j1 in 1:nt11)
         {
           w1=matrix(0,K1,K2)
           w2=w1
           w3=w1

           Delta0=matrix(0,K1,K2)
           Delta1=Delta0

           aa=ua11[i1]
           bb=ut11[j1]
           cc=tab11[j1,i1]

           for(j in 1:K1)
           {
             for(k in 1:K2)
             {
               if(j==bb & bb<=aa & j+k>aa) 
               w2[j,k]=1

               if(j+k>aa) Delta0[j,k]=1
               if(j+k>aa) Delta1[j,k]=1

               if(j+k<=aa) w3[j,k]=1
             }
           }

           q2=p1*w2
           if(sum(q2)>0) q2=cc*q2/sum(q2)

           Del0= Delta0*p0
           Del1= Delta1*p1
           Del=theta*Del0+(1-theta)*Del1

           q3=p1*w3
           if(sum(Del)>0) q3=(1-theta)*cc*q3/sum(Del)
 
           f12=f12+q2
           f13=f13+q3

           q30=p0*w3
           if(sum(Del)>0) q30=theta*cc*q30/sum(Del)
           f03=f03+q30
         }
       }

       for(i in 1:na10)
       {
         w1=matrix(0,K1,K2)
         Delta0=matrix(0,K1,K2)
         Delta1=Delta0
         ww3=w1

         aa=ua10[i]
         cc=tab10[i]

         for(j in 1:K1)
         {
           for(k in 1:K2)
           {
             if(j>aa) w1[j,k]=1
             if(j+k>aa) Delta0[j,k]=1
             if(j+k>aa) Delta1[j,k]=1
             if(j+k<=aa) ww3[j,k]=1
           }
         }

         q1=p1*w1
         if(sum(q1)>0) q1=cc*q1/sum(q1)

         q3=ww3*p1
         r1=theta*p0*Delta0+(1-theta)*p1*Delta1
         if(sum(r1)>0) q3= (1- theta)*cc*q3/sum(r1)

         f11=f11+q1
         f13=f13+q3

         q30=ww3*p0
         if(sum(r1)>0) q30= theta*cc*q30/sum(r1)

         f03=f03+q30
       }

       f0=f01+f02+f03
       f1=f11+f12+f13

       np0=f0/sum(f0)
       np1=f1/sum(f1)
       ntheta=sum(f0)/(sum(f0)+sum(f1))

       tol=max(abs(p0-np0))+max(abs(p1-np1))+abs(theta-ntheta)
       p1=np1
       p0=np0
       theta=ntheta

       o1= c(c(t(np0)),c(t(np1)),ntheta)
 
       lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
       tm1=c(tm1,lik)
     }

      o1= c(c(t(np0)),c(t(np1)),ntheta)

     lik=tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

     if(lik<tt)
     {
       so0=np0
       so1=np1
       tt=lik
     }

   }

   np0=t(so0)
   np1=t(so1)

   p0= colSums(np0)
   p1= colSums(np1)
   outg1= ntheta

   o1= c(c(t(so0)),c(t(so1)),outg1)
   alike= tlik(o1,pc,na01,nt01,na00,na11,nt11,na10,
               ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

###-------------------------------------------
### confidence interval calculation for theta
###-------------------------------------------

### upper bound

  chicalc= 0
  chicount= 0
  g1fix= 0

  while(chicalc < qchisq(0.95,df=1))
  { 
    chicount= chicount+1
    g1fix= outg1+chicount*0.001

    if(g1fix>1){g1fix= 1; break}

    alike0= cpthetafn(thetatilde= g1fix,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

    chicalc=2*(alike0-alike)

  }

  g1ub= g1fix

### lower bound

  chicalc= 0
  chicount= 0
  g1fix= 1

  while(chicalc < qchisq(0.95,df=1))
  { 
    chicount= chicount+1
    g1fix= outg1-chicount*0.001

    if(g1fix<0){g1fix= 0; break}

    alike0= cpthetafn(thetatilde= g1fix,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

    chicalc= 2*(alike0-alike)
  }

  g1lb= g1fix

###----------------------------------------
### confidence interval calculation for p0
###----------------------------------------

  ubvec0= rep(0,K1)
  lbvec0= rep(0,K1)
  ubvec1= rep(0,K1)
  lbvec1= rep(0,K1)

  ubvec_cumsum= rep(0,K1)
  lbvec_cumsum= rep(0,K1)


### upper bound for p0

  for(jj in 1:K1)
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix0= sum(np0[,jj])+chicount*0.001

      newpout= cpp0fn(pfix0,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p0ub= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p0ub-alike)

      if(pfix0 >= 1)
      {
        pfix0= 1
        break
      }
    }

    ubvec0[jj]= pfix0
  }

### lower bound for p0

  for(jj in 1:K1)
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix0= sum(np0[,jj])-chicount*0.001

      newpout= cpp0fn(pfix0,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p0lb= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p0lb-alike)

      if(pfix0 <= 0)
      {
        pfix0= 0
        break
      }
    }

    lbvec0[jj]= pfix0
  }


###----------------------------------------
### confidence interval calculation for p1
###----------------------------------------

### upper bound for p1

  for(jj in 1:K1)
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix1= sum(np1[,jj])+chicount*0.001

      newpout= cpp1fn(pfix1,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p1ub= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p1ub-alike)

      if(pfix1 >= 1)
      {
        pfix1= 1
        break
      }
    }

    ubvec1[jj]= pfix1
  }


### lower bound for p1

  for(jj in 1:K1)
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix1= sum(np1[,jj])-chicount*0.001

      newpout= cpp1fn(pfix1,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p1lb= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p1lb-alike)

      if(pfix1 <= 0)
      {
        pfix1= 0
        break
      }
    }

    lbvec1[jj]= pfix1
  }

###---------------------------------------------------
### confidence interval calculation for p0 cumulative
###---------------------------------------------------

  ubvec_cumsum0= rep(0,K1)
  lbvec_cumsum0= rep(0,K1)

  ubvec_cumsum1= rep(0,K1)
  lbvec_cumsum1= rep(0,K1)

### upper bound for p0

  for(jj in 2:(K1-1))
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix0= sum(np0[,1:jj])+chicount*0.001

      newpout= cpp0fn_cumsum(pfix0,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p0ub= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p0ub-alike)

      if(pfix0 >= 1)
      {
        pfix0= 1
        break
      }
    }

    ubvec_cumsum0[jj]= pfix0
  }


### lower bound for p0

  for(jj in 2:(K1-1))
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix0= sum(np0[,1:jj])-chicount*0.001

      newpout= cpp0fn_cumsum(pfix0,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p0lb= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p0lb-alike)

      if(pfix0 <= 0)
      {
        pfix0= 0
        break
      }
    }

    lbvec_cumsum0[jj]= pfix0
  }


###---------------------------------------------------
### confidence interval calculation for p1 cumulative
###---------------------------------------------------

### upper bound for p1

  for(jj in 2:(K1-1))
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix1= sum(np1[,1:jj])+chicount*0.001

      newpout= cpp1fn_cumsum(pfix1,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p1ub= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p1ub-alike)

      if(pfix1 >= 1)
      {
        pfix1= 1
        break
      }
    }

    ubvec_cumsum1[jj]= pfix1
  }


### lower bound for p1

  for(jj in 2:(K1-1))
  {
    chicalc= 0
    chicount= 0

    while(chicalc < qchisq(0.95,df=1) & chicount<=200)
    {
      chicount= chicount+1

      pfix1= sum(np1[,1:jj])-chicount*0.001

      newpout= cpp1fn_cumsum(pfix1,jj,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)
      newp0= newpout[[1]]
      newp1= newpout[[2]]
      newtheta= newpout[[3]]

      newso0= t(newp0)
      newso1= t(newp1)

      oci= c(c(t(newso0)),c(t(newso1)),newtheta)
      glikfn_p1lb= tlik(oci,pc,na01,nt01,na00,na11,nt11,na10,
                ua01,ut01,tab01,ua00,tab00,ua11,ut11,tab11,ua10,tab10)

      chicalc= 2*(glikfn_p1lb-alike)

      if(pfix1 <= 0)
      {
        pfix1= 0
        break
      }
    }

    lbvec_cumsum1[jj]= pfix1
  }

  out_cumsum0_lb= c(lbvec0[1],lbvec_cumsum0[2:(K1-1)])
  out_cumsum0_ub= c(ubvec0[1],ubvec_cumsum0[2:(K1-1)])
  out_cumsum1_lb= c(lbvec1[1],lbvec_cumsum1[2:(K1-1)])
  out_cumsum1_ub= c(ubvec1[1],ubvec_cumsum1[2:(K1-1)])

  outtheta= 1-outg1
  outthetalb= 1-g1ub
  outthetaub= 1-g1lb

  myout= list(p0,lbvec0,ubvec0,out_cumsum0_lb,out_cumsum0_ub,
              p1,lbvec1,ubvec1,out_cumsum1_lb,out_cumsum1_ub,
              outtheta,outthetalb,outthetaub)
  names(myout)= c("p0.est","p0.lower","p0.upper","p0.cumulative.lower","p0.cumulative.upper",
                  "p1.est","p1.lower","p1.upper","p1.cumulative.lower","p1.cumulative.upper",
                  "theta.est","theta.lower","theta.upper")
  return(myout)

}


### A data example for the proposed EM algorithm
### input list

pc= list(K2= 4,                   ## number of categories for T2 (time from diagnosis to death)
         t1cuts= c(0.3,0.6,1.0),  ## cut-points for T1 (age at diagnosis)
         acuts= c(0.05,0.15,0.35),## cut-points for A (age at study)
         ageatdx= 'AGEATDX',      ## column name in data file for age at diagnosis
                                  ## (note: must be NA for non-cases i.e. healthy subjects)
         ageatstudy= 'AGEATSTUDY',## column name in data file for age at study
         mutation= 'MUTATION'     ## column name in data file for mutation status 
                                  ## (coded as: 1 for mutation present, 0 otherwise)
         )

input.data= read.csv("data-sim.csv", header=T)

output.em= func.em(pc,input.data)
output.em$'p1.est'                ## to extract the estimates of p1 
