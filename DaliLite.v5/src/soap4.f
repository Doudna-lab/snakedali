c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c------------------------------------------------------------------------------
c
	subroutine initequal_fast_domains(table0,nres1,nres2,ali,
     $		domnres1,domres1,domnres2,domres2)
	implicit none
	integer nres1,nres2,ali(nres1)
	integer domnres1,domres1(nres1),domnres2,domres2(nres2)
	integer*2 table0(nres1,nres2)
c
	integer a,b,i,j
	real r,fa,fb

c	write(*,*),'init',nres1,nres2,domnres1,domnres2
	do i=1,nres1
		ali(i)=0
	end do
	if(domnres1.lt.domnres2) then
		r=float(domnres1)/float(domnres2)
	else
		r=float(domnres2)/float(domnres1)
	end if
	fa=1.0
	fb=1.0
	a=0
	b=0
	i=0
	do while(a.lt.domnres1)
		if(domnres1.lt.domnres2) then
			fa=fa+r
			fb=fb+1.0
		else
			fa=fa+1.0
			fb=fb+r
		end if
		a=int(fa)
		b=int(fb)
		if(a.gt.0.and.b.gt.0) then
			if(ali(domres1(a)).eq.0) ali(domres1(a))=domres2(b)
		end if
	end do
	do i=1,nres1
		do j=1,nres2
			table0(i,j)=0	! penalty
		end do
	end do

	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine initequal_fast(table0,nres1,nres2,ali,icyc)
	implicit none
	integer nres1,nres2,ali(nres1),icyc
	integer*2 table0(nres1,nres2)
c
	integer a,b,i,j
	real r,fa,fb
c
	write(*,*),'init',nres1,nres2,icyc
	do i=1,nres1
		ali(i)=0
	end do
	if(icyc.eq.1) then
		if(nres1.lt.nres2) then
			r=float(nres1)/float(nres2)
		else
			r=float(nres2)/float(nres1)
		end if
		fa=1.0
		fb=1.0
		a=0
		b=0
		i=0
		do while(a.lt.nres1)
			if(nres1.lt.nres2) then
				fa=fa+r
				fb=fb+1.0
			else
				fa=fa+1.0
				fb=fb+r
			end if
			a=int(fa)
			b=int(fb)
			if(a.gt.0.and.b.gt.0) then
				if(ali(a).eq.0) ali(a)=b
			end if
		end do
	else if(icyc.eq.2) then
		do i=1,min(nres1,nres2)
			ali(i)=i
		end do
	else if(icyc.eq.3) then
		if(nres1.gt.nres2) then
			j=nres1-nres2
			do i=1,nres2
				ali(i+j)=i
			end do
		else
			j=nres2-nres1
			do i=1,nres1
				ali(i)=i+j
			end do
		end if
	end if
	do i=1,nres1
		do j=1,nres2
			table0(i,j)=0
		end do
	end do

	return
	end
c
c------------------------------------------------------------------------------
c
c
c------------------------------------------------------------------------------
c
	subroutine initequal(table0,nres1,nres2,ali,icyc)
	implicit none
	integer nres1,nres2,ali(nres1+nres2),icyc
	integer*2 table0(nres1,nres2)
c
	integer a,b,olda,oldb,i,j,k
	real r,fa,fb
c
	if(icyc.eq.1) then
		if(nres1.lt.nres2) then
			r=float(nres1)/float(nres2)
		else
			r=float(nres2)/float(nres1)
		end if
		fa=1.0
		fb=1.0
		olda=0
		oldb=0
		i=0
		do while(i.lt.nres1+nres2)
			if(nres1.lt.nres2) then
				fa=fa+r
				fb=fb+1.0
			else
				fa=fa+1.0
				fb=fb+r
			end if
			a=int(fa)
			b=int(fb)
			if(a.gt.olda) then
				olda=a
				i=i+1
				ali(i)=1
			end if
			if(b.gt.oldb) then
				oldb=b
				i=i+1
				ali(i)=2
			end if
		end do
	else if(icyc.eq.2) then
		i=0
		do j=1,min(nres1,nres2)
			i=i+1
			ali(i)=1
			i=i+1
			ali(i)=2
		end do
		if(nres1.gt.nres2) then
			k=1
		else
			k=2
		end if
		do j=i+1,nres1+nres2
			ali(j)=k
		end do
	else if(icyc.eq.3) then
		if(nres1.gt.nres2) then
			do i=1,nres1-nres2
				ali(i)=1
			end do
			i=nres1-nres2
		else
			do i=1,nres2-nres1
				ali(i)=2
			end do
			i=nres2-nres1
		end if
		do while (i.lt.nres1+nres2)
			i=i+1
			ali(i)=1
			i=i+1
			ali(i)=2
		end do
	end if
	do i=1,nres1
		do j=1,nres2
			table0(i,j)=100
		end do
	end do

	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine fill1(nres1,nres2,ca1,ca2,table0,rcutoff)
	implicit none
	include 'parsizes.for'
	integer nres1,nres2,i,j
	real ca1(3,nres1),ca2(3,nres2),distance,r,rcutoff
	integer*2 table0(nres1,nres2)

	do i=1,nres1
		do j=1,nres2
			r=distance(ca1(1,i),ca1(2,i),ca1(3,i),
     $				ca2(1,j),ca2(2,j),ca2(3,j))
			! minimize penalty
			table0(i,j)=-nint(100.0*exp(-r/rcutoff))
		end do
	end do

c	do i=1,20
c		write(*,520) (table0(i,j),j=1,20)
c	end do

520	format(20i4)

	return
	end
c
c------------------------------------------------------------------------------
c
 	subroutine fill1_fast(nres1,nres2,ca1,ca2,table0,rcutoff)
	implicit none
	include 'parsizes.for'
	integer nres1,nres2,i,j
	real ca1(3,nres1),ca2(3,nres2),squaredistance,r,rcutoff,r2
	integer*2 table0(nres1,nres2)

	r2=rcutoff*rcutoff
	do i=1,nres1
		do j=1,nres2
			r=squaredistance(ca1(1,i),ca1(2,i),ca1(3,i),
     $				ca2(1,j),ca2(2,j),ca2(3,j))
			! minimize penalty
			if(r.ge.r2) then
				table0(i,j)=100
			else
				table0(i,j)=nint(100.0*r/r2)
			end if
		end do
	end do

c	do i=1,20
c		write(*,520) (table0(i,j),j=1,20)
c	end do

520	format(20i4)

	return
	end
c
c------------------------------------------------------------------------------
c
       function distance(a1,a2,a3,b1,b2,b3)
        implicit none
        real distance,a1,a2,a3,b1,b2,b3
c
        distance=sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3))
c
        return
        end
c
c------------------------------------------------------------------------------
c
       function squaredistance(a1,a2,a3,b1,b2,b3)
        implicit none
        real squaredistance,a1,a2,a3,b1,b2,b3
c
        squaredistance=(a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3)
c
        return
        end
c
c----------------------------------------------------------------------
c
	subroutine soap_getut(nx,ali,x,y,u,t,lali,rms,table0,nres1,
     $          nres2,totw)
        implicit none
	include 'parsizes.for'
        integer nx,ali(nx),lali,nres1,nres2
        real u(3,3),t(3),rms,x(3,*),y(3,*)
	integer*2 table0(nres1,nres2)
c
        integer i,j,ier,a,b,k,s
        real ux(3,maxres+maxres),uy(3,maxres+maxres),w(maxres+maxres)
	real ssq,totw
c
        lali=0
	a=1
	b=1
	totw=0.0
        do i=2,nx-2
		k=ali(i)
		if(k.eq.0) goto 19
 		if(k.eq.1) then
			a=a+1
		else if(k.eq.2) then
			b=b+1
		else if(k.eq.3) then
			a=a+1
			b=b+1
		end if
		s=table0(a,b)
		if(s.eq.100) goto 19
		lali=lali+1
c		write(*,*),'pair:',lali,k,a,b,s
		do j=1,3
			if(k.eq.1) then
 				ux(j,lali)=(x(j,a)+x(j,a-1))/2
			else
 				ux(j,lali)=x(j,a)
			end if
			if(k.eq.2) then
				uy(j,lali)=(y(j,b)+y(j,b-1))/2 ! testing
			else
				uy(j,lali)=y(j,b)
			end if
			w(lali)=1.0-float(s)/100.0 ! min.penalty
			totw=totw+w(lali)
		end do
19 	end do
        if(lali.gt.2) then
		call u3b(w,ux,uy,lali,1,ssq,u,t,ier)
	else
		ssq=0.0
	end if
        rms=sqrt(ssq/max(1.0,totw))

        return
        end
c
c-------------------------------------------------------------------------------
c
	subroutine getut_w(nx,ali,x,y,u,t,lali,rms,table0,nres1,nres2)
        implicit none
	include 'parsizes.for'
        integer nx,ali(maxres),lali,nres1,nres2
        real u(3,3),t(3),rms,x(3,*),y(3,*)
	integer*2 table0(nres1,nres2)
c
        integer i,j,ier,a,b
        real ux(3,maxres),uy(3,maxres),ssq,w(maxres),totw
c
        lali=0
	totw=0.0
	do a=1,nx
		b=ali(a)
		if(b.eq.0) goto 19
		if(table0(a,b).ge.99) goto 19		! high penalty
		lali=lali+1
c		write(*,*),'pair:',lali,ali(i),a,b
		do j=1,3
			ux(j,lali)=x(j,a)
			uy(j,lali)=y(j,b)
		end do
		w(lali)=abs(float(100-table0(a,b))/100.0) ! minimizing penalty
		totw=totw+w(lali)
19 	end do
        call u3b(w,ux,uy,lali,1,ssq,u,t,ier)
        rms=sqrt(ssq/max(1.0,totw))

        return
        end
c
c-------------------------------------------------------------------------------
c
	subroutine soap_nw(m,n,table0,ali,trace,table)
	implicit none
	integer m,n
	integer*2 trace(m,n),table0(m,n)
	integer table(m,n),ali(m+n),score,a,b,infinit,v
	parameter(infinit=1e9)
c
	integer i,j,x,y
c
	do i=1,m+n
		ali(i)=0
	end do
	score=0
c
c	fill trace
c
	trace(1,1)=0
	table(1,1)=table0(1,1)
	do i=1,m
		do j=1,n
		  if(i.gt.1.or.j.gt.1) then
			if(i.gt.1) then
				a=table(i-1,j)
			else
				a=infinit
			end if
			if(j.gt.1) then
				b=table(i,j-1)
			else
				b=infinit
			end if
			if(a.lt.b) then
				trace(i,j)=1
				table(i,j)=table0(i,j)+a
			else
				trace(i,j)=2
				table(i,j)=table0(i,j)+b
			end if
		  end if
		end do
	end do
c
c
c
c	do i=1,20
c		write(*,500) (table(i,j),j=1,15)
c	end do
c
c	backtrack alignment
c
	x=m
	y=n
	v=m+n
	do while(x.gt.1.or.y.gt.1)
		i=trace(x,y)
		ali(v)=i
		v=v-1
		if(i.eq.1) then
			x=x-1
		else if(i.eq.2) then
			y=y-1
		else
			!stop 'unknown trace'
                        return
		end if
	end do
c
500	format(15i5)
c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine soap_nw_fast(m,n,table0,ali,trace)
	implicit none
	include 'parsizes.for'
	integer m,n
	integer*2 trace(m,n),table0(m,n),s
	integer shift(0:2,0:maxres)
	integer ali(m+n)
c
	integer i,j,k,x,y,y0,a,b,c,k0,k1,diag0,diag1,diag2,v
c
c	fill trace
c
	diag0=2
	diag1=1
	diag2=0
	trace(1,1)=0
	trace(2,1)=1
	trace(1,2)=2
	shift(diag2,0)=table0(1,1)
	shift(diag1,0)=table0(2,1)+shift(diag2,0)
	shift(diag1,1)=table0(1,2)+shift(diag2,0)
	do y0=3,m+n-1
		k0=max(0,y0-n)		! x==1
		x=1+k0
		y=y0-k0
		trace(x,y)=2
		shift(diag0,k0)=shift(diag1,k0)+table0(x,y)
		k1=min(m-1,y0-1)	! y==1
		x=1+k1
		y=y0-k1
		trace(x,y)=1
		shift(diag0,k1)=shift(diag1,k1-1)+table0(x,y)
		do k=k0+1,k1-1		! x>1,y>1
			x=1+k
			y=y0-k
			s=table0(x,y)
			a=shift(diag1,k-1)
			b=shift(diag1,k)
			c=shift(diag2,k-1)
			if(a.le.b.and.a.le.c) then
				trace(x,y)=1
				shift(diag0,k)=a+s
			else if(b.le.a.and.b.le.c) then
				trace(x,y)=2
				shift(diag0,k)=b+s
			else
				trace(x,y)=3
				shift(diag0,k)=c+s
			end if
		end do
		diag0=mod(diag0+1,3)			! rolling indices
		diag1=mod(diag1+1,3)
		diag2=mod(diag2+1,3)
	end do
c
c
c
!	do i=1,20
!		write(*,500) (table0(i,j),j=1,15)
!	end do
! 	do i=1,20
!		write(*,500) (trace(i,j),j=1,15)
!	end do
c
c	backtrack alignment
c
!	do i=1,m
!		ali(i)=0
!	end do
!	x=m
!	y=n
!	do while(x.gt.1.and.y.gt.1)
!		i=trace(x,y)
!		if(i.eq.1) then
!			x=x-1
!		else if(i.eq.2) then
!			y=y-1
!		else if(i.eq.3) then
!			ali(x)=y
!			x=x-1
!			y=y-1
!		else
!			type 500,(ali(i),i=1,m)
!			stop 'unknown trace'
!		end if
!	end do
c
c	backtrack alignment
c
	do i=1,m+n
		ali(i)=0
	end do
	x=m
	y=n
	v=m+n
	do while(x.gt.1.or.y.gt.1)
		i=trace(x,y)
		ali(v)=i
		v=v-1
		if(i.eq.1) then
			x=x-1
		else if(i.eq.2) then
			y=y-1
		else if(i.eq.3) then
			x=x-1
			y=y-1
			ali(v+1)=1		! patch 3->(1,2)
			ali(v)=2
			v=v-1
		else
			!stop 'unknown trace'
                        return
		end if
	end do
c
500	format(15i5)
c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine getlist(filnam,iunit,list,nprot)
	implicit none
	include 'parsizes.for'
	character*5 list(maxprot)
	integer nprot,iunit
	character*(*) filnam
c
	character*5 cd
c
	nprot=0
	open(iunit,file=filnam,status='old')
10	read(iunit,500,end=19) cd
	if(nprot.eq.maxprot) then
c		write(*,*) 'WARNING: skip reading list after maxprot',maxprot
		goto 19
	end if
	nprot=nprot+1
	list(nprot)=cd
	goto 10
19	close(iunit)
c	write(*,*) nprot,' proteins in list from ',filnam
c
500	format(a5)
c

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine setup(cd1,ndom,domns,domseglist,nres,ca,iunit,
     $		dalidatpath1,dalidatpath2)
	implicit none
	include 'parsizes.for'
	integer ndom,domns(maxdom),domseglist(maxseg,2,maxdom),iunit,nres
	real ca(3,maxres)
	character*5 cd1
c
	character*80 filnam,line,constructfilnam
        character*80 dalidatpath1,dalidatpath2
	character node_type(maxdom)
	integer i,j,k,idom,nseg
c
!	write(*,*),'setup ',cd1
	filnam=constructfilnam(cd1,dalidatpath1,'.dat')
	open(iunit,file=filnam,status='old',err=11)
        goto 12
11      write(*,*) 'file not found, try other path ',dalidatpath2
        filnam=constructfilnam(cd1,dalidatpath2,'.dat')
        open(iunit,file=filnam,status='old',err=19)
12      continue
		i=0
200		read(iunit,530) line
		if(line(1:4).ne.'>>>>') goto 200
		read(line,550) ndom
		i=i+1
		if(i.eq.1) then
			read(line(11:20),*) nres,nseg
			do j=1,nseg
				read(iunit,530) line
			end do
			read(iunit,*) ((ca(k,j),k=1,3),j=1,nres)
		end if
		if(i.lt.3) goto 200
		do idom=1,ndom
			read(iunit,710,end=219) j,node_type(j),domns(j),
     $				((domseglist(i,k,j),k=1,2),i=1,domns(j))
		end do
219	close(iunit)
c
c	keep only (+,*)-units
c
	ndom=0
	do i=1,j
		if(i.eq.1.or.node_type(i).eq.'+'.or.node_type(i).eq.'*') then
			ndom=ndom+1
			domns(ndom)=domns(i)
			do k=1,domns(i)
				domseglist(k,1,ndom)=domseglist(k,1,i)
				domseglist(k,2,ndom)=domseglist(k,2,i)
			end do
		end if
	end do
c
c	normal exit
c
	return
c
c	error exit
c
19	write(*,*) 'ERROR in fetchdomseglist: could not open file',filnam
	ndom=0
c
530	format(a80)
550	format(10x,i5)
710	format(i4,1x,a1,13x,400i4)

	return
	end
c
c------------------------------------------------------------------------------
c
