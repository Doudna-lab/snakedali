c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c------------------------------------------------------------------------------
c
        subroutine getut(nx,ali,x,y,u,t,lali,rms)
        implicit none
	include 'parsizes.for'
        integer nx,ali(nx),lali
        real u(3,3),t(3),rms,x(3,*),y(3,*)
c
        integer i,j,ier
        real ux(3,maxres),uy(3,maxres),ssq,w(maxres)
c
        lali=0
        do i=1,min(maxres,nx)
                if(ali(i).ne.0) then
                        lali=lali+1
                        do j=1,3
                                ux(j,lali)=x(j,i)
                                uy(j,lali)=y(j,ali(i))
                                w(lali)=1.0
                        end do
                end if
        end do
        call u3b(w,ux,uy,lali,1,ssq,u,t,ier)
        rms=sqrt(ssq/max(1,lali))

        return
        end
c
c-------------------------------------------------------------------------------
c
        subroutine transrotate(x,nx,u,t)
c
c       overwrites x with (ux+t)
c
        implicit none
        integer nx
        real x(3,nx),u(3,3),t(3),tmp(3)
c
        integer i,j,k
c
        do i=1,nx
                do j=1,3
                        tmp(j)=t(j)
                        do k=1,3
                                tmp(j)=tmp(j)+u(j,k)*x(k,i)
                        end do
                end do
                do j=1,3
                        x(j,i)=tmp(j)
                end do
        end do

        return
        end
c
c-------------------------------------------------------------------------------
c
	subroutine fitz(x,nx,y,ny,ali,rcut,maxiter,rms,lali,iter)
c
c	find trace of closest neighbours within rcut:
c		score(i,j)=max(rcut-dist(i,j),0.0)
c	repeat u3b until ali freezes or maxiter
c
c	returns last ali, superimposed x,y coordinates
c
	implicit none
	include 'parsizes.for'	! only for maxres
	integer nx,ny,ali(nx),maxiter
	real rcut,x(3,nx),y(3,ny)
	integer s,iter,i
	logical lfrozen
	integer aliold(maxres)
	integer*2 trace(maxres*maxres),table0(maxres*maxres)
	real rms,u(3,3),t(3)
	integer lali
c
c	check dimensions
c
        if(nx.lt.1.or.ny.lt.1) return ! HACK: error exit
	if(nx.gt.maxres.or.ny.gt.maxres) then
                write(*,*),' WARNING: maxres overflow:',nx,ny,maxres
		nx=min(nx,maxres)
		ny=min(ny,maxres)
	end if
!	type *,'fitz:',nx,ny,rcut,maxiter
c
	iter=0
	lfrozen=.false.
	do i=1,nx
		ali(i)=0
	end do
	do while(iter.lt.maxiter.and..not.lfrozen)
		iter=iter+1
		call filltable_maxsim(nx,ny,table0,x,y,rcut)
		do i=1,nx
			aliold(i)=ali(i)
		end do
		call nw_maxsim(nx,ny,table0,s,ali,trace)
!		type 500,s
!		write(*,*) 'nw ali ',nx,(ali(i),i=1,nx)
		lfrozen=.true.
		i=0
		do while(lfrozen.and.i.lt.nx)
			i=i+1
			lfrozen=(ali(i).eq.aliold(i))
		end do
c
c		do new superimposition
c
		call getut(nx,ali,x,y,u,t,lali,rms)
!		write(*,510) rms,lali,iter,s,nx,(ali(i),i=1,nx)
		call transrotate(x,nx,u,t)
	end do

500	format('score:',i10)
505	format(20i4)
510	format('rmsd:',f10.1,' lali:',i5,' iteration:',i5,' score:',i20)

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine filltable_maxsim(nx,ny,table0,x,y,rcut)
	implicit none
	integer nx,ny
	integer*2 table0(nx,ny)
	real x(3,nx),y(3,ny),rcut
c
	real distance_real
c
	integer i,j
	real r
c
	do i=1,nx
		do j=1,ny
		r=distance_real(x(1,i),x(2,i),x(3,i),y(1,j),y(2,j),y(3,j))
			if(r.ge.rcut) then
				table0(i,j)=0
			else
				table0(i,j)=nint((rcut-r)*10.0)
			end if
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
        subroutine filltable_maxsim_cumulative(nx,ny,table0,x,y,rcut)
        implicit none
        integer nx,ny
        integer*2 table0(nx,ny)
        real x(3,nx),y(3,ny),rcut
c
        real distance_real
c
        integer i,j
        real r
c
        do i=1,nx
          do j=1,ny
            r=distance_real(x(1,i),x(2,i),x(3,i),y(1,j),y(2,j),y(3,j))
            if(r.lt.rcut) then
              table0(i,j)=max(table0(i,j),nint((rcut-r)*10.0))
            end if
          end do
        end do

        return
        end

c
c-------------------------------------------------------------------------------
c
	subroutine nw_maxsim(m,n,table0,score,ali,trace)
	implicit none
	include 'parsizes.for'
	integer m,n
	integer*2 trace(m,n),table0(m,n)
	integer ali(m),score
c
	integer shift(0:2,0:maxres)
	integer i,j,k,x,y,y0,a,b,c,k0,k1,diag0,diag1,diag2,s
c
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
	shift(diag1,0)=max(table0(2,1),table0(1,1))
	shift(diag1,1)=max(table0(1,2),table0(1,1))
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
			if(c.ge.a.and.c.ge.b) then	! maximize similarity
				trace(x,y)=3
				shift(diag0,k)=c+s
			else if(b.ge.a.and.b.ge.c) then
				trace(x,y)=2
				shift(diag0,k)=b
			else
				trace(x,y)=1
				shift(diag0,k)=a
			end if
		end do
		diag0=mod(diag0+1,3)			! rolling indices
		diag1=mod(diag1+1,3)
		diag2=mod(diag2+1,3)
	end do
c
c
c
c	do i=1,20
c		write(*,500) (table0(i,j),j=1,15)
c	end do
c 	do i=1,20
c		write(*,500) (trace(i,j),j=1,15)
c	end do
c
c	backtrack alignment
c
	do i=1,m
		ali(i)=0
	end do
	x=m
	y=n
	do while(x.gt.1.or.y.gt.1)
		i=trace(x,y)
		if(i.eq.1) then
			x=x-1
		else if(i.eq.2) then
			y=y-1
		else if(i.eq.3) then
			if(table0(x,y).gt.0) then
				ali(x)=y	! clean-up hack
			end if
			x=x-1
			y=y-1
		else
			y=0
			x=0 !stop 'unknown trace'
		end if
	end do
	score=0
	do i=1,m
		if(ali(i).gt.0) score=score+table0(i,ali(i))
	end do
c
500	format(25i3)
c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine filltable_minpen(nx,ny,table0,x,y,rcut,t0)
	implicit none
	integer nx,ny
	integer*2 table0(nx,ny),t0
	real x(3,nx),y(3,ny),rcut
c
	real distance
c
	integer i,j
	real r
c
	do i=1,nx
		do j=1,ny
			r=distance(x(1,i),x(2,i),x(3,i),y(1,j),y(2,j),y(3,j))
			if(r.ge.rcut) then
				table0(i,j)=t0
			else
				table0(i,j)=nint(r*10.0)
			end if
		end do
	end do

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine nw_minpen(m,n,table0,score,ali,trace,t0)
	implicit none
	include 'parsizes.for'
	integer m,n
	integer*2 trace(m,n),table0(m,n),t0
	integer ali(m),score
c
	integer shift(0:2,0:maxres)
	integer i,j,k,x,y,y0,a,b,c,k0,k1,diag0,diag1,diag2,s
c
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
			if(c.le.a.and.c.le.b) then	! minimize penalty
				trace(x,y)=3
				shift(diag0,k)=c+s
			else if(b.le.a.and.b.le.c) then
				trace(x,y)=2
				shift(diag0,k)=b
			else
				trace(x,y)=1
				shift(diag0,k)=a
			end if
		end do
		diag0=mod(diag0+1,3)			! rolling indices
		diag1=mod(diag1+1,3)
		diag2=mod(diag2+1,3)
	end do
c
c
c
c	do i=1,20
c		write(*,500) (table0(i,j),j=1,15)
c	end do
c 	do i=1,20
c		write(*,500) (trace(i,j),j=1,15)
c	end do
c
c	backtrack alignment
c
	do i=1,m
		ali(i)=0
	end do
	score=0
	x=m
	y=n
	do while(x.gt.1.or.y.gt.1)
		i=trace(x,y)
		if(i.eq.1) then
			x=x-1
		else if(i.eq.2) then
			y=y-1
		else if(i.eq.3) then
			if(table0(x,y).lt.t0) then
				ali(x)=y	! clean-up hack
				score=score+table0(x,y)
			end if
			x=x-1
			y=y-1
		else
			y=0
			x=0 !stop 'unknown trace'
		end if
	end do
c
500	format(25i3)
c
	return
	end
c
c------------------------------------------------------------------------------
c
        function distance_real(a1,a2,a3,b1,b2,b3)
        implicit none
        real distance_real,a1,a2,a3,b1,b2,b3
c
        distance_real=
     $		sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3))
c
        return
        end
c
c----------------------------------------------------------------------
c

