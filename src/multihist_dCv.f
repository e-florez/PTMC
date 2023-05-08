        implicit real*8(a-h,o-z)
	parameter(OUT=1e35)
	parameter(NSIZEMAX=1000,NTMAX=100,NENRMAX=1000)
	parameter(NVMAX=10000,nmarge=2,nini=2)
	parameter(pi=3.141592654)
	real*8 kB
	real*8 emin,emax
	double precision v,vmin,vmax,dv,s,xp,pivot,u,r2,r3,T(NTMAX)
	dimension wn(NTMAX,NVMAX),wnpa(NVMAX),nmax(NTMAX)
	dimension jpmax(NTMAX)
	double precision b(NTMAX,NVMAX),Etot(NTMAX)
	double precision AM(NTMAX,NTMAX),BM(NTMAX),x(NTMAX),y(NVMAX)
	double precision xnpa(NVMAX),tmin,tmax,tl,dt
	double precision Z,Cv,dCv,vminnew,vmaxnew,kadm

	double precision kdib(NVMAX)
	double precision etotmin,etotmax,detot,eetot

	character*29 filename,word,filerec,filedipole

	character*2 nb(60)
	data nb/'1','2','3','4','5','6','7','8','9','10',
     #       '11','12','13','14','15','16','17','18','19','20',
     #       '21','22','23','24','25','26','27','28','29','30',
     #       '31','32','33','34','35','36','37','38','39','40',
     #       '41','42','43','44','45','46','47','48','49','50',
     #       '51','52','53','54','55','56','57','58','59','60'/

	open(19,file='histE.data',status='unknown')
        read(19,*)
        read(19,*) kb
        read(19,*) NT
        read(19,*) (T(i),i=1,NT)
        read(19,*) emin,emax,nhe
	close(19)

	de=(emax-emin)/(NHE-1)
	NV=NHE

	eminf=emin
	emaxf=emax

	do it=1,NT
	   do j=1,NV
	      wn(it,j)=0.
	      b(it,j)=0.
	   enddo

	   open(10,file='histE.'//nb(it))
           do i=1,NV
              read(10,*) edum,nh
	      wn(it,i)=wn(it,i)+nh*1.D0
	      wnsum=wnsum+1.D0
           enddo
	   close(10)
	   do j=1,NV
              ener=(j-1)*de+emin
	      wn(it,j)=wn(it,j)/wnsum
	      b(it,j)=log(wn(it,j))+ener/T(it)/kB
	   enddo

	enddo

	print *,'n(i,j) found.'

1985	do j=1,NV
	   wnpa(j)=0
	   do i=1,NT
	      wnpa(j)=wnpa(j)+wn(i,j)
	   enddo
	   if (wnpa(j).eq.0.) wnpa(j)=1.
	   xnpa(j)=wnpa(j)
	enddo

	do i=1,NT
	   BM(i)=0.
	   do j=1,NV
	      if (wn(i,j).ne.0.) BM(i)=BM(i)+b(i,j)*wn(i,j)
	      xp=0.
	      do k=1,NT
		 if (wn(k,j).ne.0.) xp=xp+wn(k,j)*b(k,j)
	      enddo
	      BM(i)=BM(i)-xp*wn(i,j)/xnpa(j)
	   enddo
	enddo

	print *,'B vector calculated.'
	s=0.
	do i=1,NT
	   s=s+BM(i)
	enddo
	print *,'   sum=',s
	
	do i=1,NT
	   do ip=1,NT
	      AM(i,ip)=0.
	      do k=1,NV
		 AM(i,ip)=AM(i,ip)-wn(i,k)*wn(ip,k)/wnpa(k)
	      enddo
	      
	      if (i.eq.ip) then
		 do k=1,NV
		    AM(i,ip)=AM(i,ip)+wn(i,k)
		 enddo
	      endif
	      
	   enddo
	
	enddo
	
	print *,'A matrix calculated'

	do j=1,NT
	   s=0.
	   do i=1,NT
	      s=s+AM(i,j)
	   enddo
	   print *,'Column ',j,':',s
	enddo

	do i=1,NT-1
	   
	   jpivot=0
	   pivot=0.
	   do j=i,NT
	      u=AM(j,i)
	      if (pivot.lt.abs(u)) then
		 pivot=abs(u)
		 jpivot=j
	      endif
	   enddo
	   j=jpivot
	   u=AM(j,i)

	   do k=1,NT
	      s=AM(i,k)
	      AM(i,k)=AM(j,k)
	      AM(j,k)=s
	   enddo
	   s=BM(i)
	   BM(i)=BM(j)
	   BM(j)=s
	   
	   do k=i+1,NT
	      s=AM(k,i)	
	      AM(k,i)=0.
	      do nk=i+1,NT
		 AM(k,nk)=AM(k,nk)-s*AM(i,nk)/u
	      enddo
	      BM(k)=BM(k)-s*BM(i)/u
	   enddo
	   
	   do k=1,i-1
	      s=AM(k,i)
	      AM(k,i)=0.
	      do nk=i+1,NT
		 AM(k,nk)=AM(k,nk)-s*AM(i,nk)/u
	      enddo
	      BM(k)=BM(k)-s*BM(i)/u
	   enddo

	enddo

	AM(NT,NT)=1.
	BM(NT)=0.
	
	x(NT)=BM(NT)/AM(NT,NT)

	do i=NT-1,1,-1
	   s=BM(i)
	   x(i)=s/AM(i,i)
	enddo

	print *,'System solved'
	
	open(17,file='solx',status='unknown')
	do j=1,NT
	   write(17,*) j,x(j)
	enddo
	close(17)

	do j=1,NV
	   y(j)=0.
	   do i=1,NT
	      if (wn(i,j).ne.0.) y(j)=y(j)+wn(i,j)*(b(i,j)-x(i))
	   enddo
	   y(j)=y(j)/wnpa(j)
	enddo

	open(20,file='S',status='unknown')
	do i=1,NV
	   if (y(i).ne.0) write(20,*) emin+(i-1)*de,y(i)
	enddo
	close(20)

	print *,'Entropy saved'

	Tmin=T(1)
	Tmax=T(NT)
	nexp=-1000000000
	npoints=1000
	dT=(Tmax-Tmin)/npoints

	open(11,file='U.NVT',status='unknown')
	open(12,file='Cv.NVT',status='unknown')
	open(13,file='dCv.NVT',status='unknown')

	do i=1,NV
	   if (y(i).ne.0) then
	      npopo=int(y(i)-(emin+(i-1)*de)/Tmin/kb)
	      if (npopo.gt.nexp) nexp=npopo
	   endif
	enddo
	
	do it=1,npoints
	   Tl=tmin+(it-1)*dT
	   
 2001	   Z=0.
	   do j=1,NV
	      if (y(j).ne.0) then
		 Z=Z+exp(y(j)-(emin+(j-1)*de)/Tl/kb-nexp*1.)
	      endif
	   enddo
	   if (Z.lt.1) then
	      nexp=nexp-2
	      goto 2001
	   endif
	   if (Z.gt.100) then
	      nexp=nexp+2
	      goto 2001
	   endif
	   
	   xp=0
	   Z=0.
	   u=0.
	   s=0.
	   do j=1,NV
	      v=emin+(j-1)*de
	      if (y(j).ne.0) then
		 xp=exp(y(j)-v/Tl/kb-nexp*1.)
	      else
		 xp=0.0
	      endif
	      Z=Z+xp
	      u=u+v*xp   !E
	      s=s+v*v*xp !E^2
	   enddo
	   u=u/Z
	   s=s/Z
	   Cv=(s-u*u)/kb/Tl**2

	   write(11,*) Tl,u
	   write(12,*) Tl,Cv
	   
           !dCv/dT
           r2=0.
           r3=0.
           do j=1,NV
              v=emin+(j-1)*de
              if (y(j).ne.0) then
                 xp=exp(y(j)-v/Tl/kb-nexp*1.)
              else
                 xp=0.0
              endif
              Z=Z+xp
              r2=r2+(v-u)**2*xp
              r3=r3+(v-u)**3*xp
           enddo
           r2=r2/Z      
           r3=r3/Z      
           dCv=r3/kb**2/Tl**4 - 2*r2/kb/Tl**3

           write(13,*) Tl,dCv
   
        enddo
	
	close(11)
	close(12)
	close(13)
	
	print *,'Canonical data saved.'

	stop

	end

