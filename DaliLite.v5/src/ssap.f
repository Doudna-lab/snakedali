c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c------------------------------------------------------------------------------
c
c	iteration:
c		trim net-negative assignments from prealignment
c		fill table0 with scores(i,j) relative to prealignment
c		needle
c		use best alignment as new prealignment
c	until bestscore no longer increases
c
c------------------------------------------------------------------------------
c
	subroutine trimming(nres1,nres2,preali1,d1,d2)
	implicit none
	include 'gagasizes.for'
	integer nres1,nres2
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres)
c
	real scorefun,s(maxres)
	integer i,j,k,l,n,a(maxres)
c
c	write(*,*) 'this is trimming',nres1,nres2
	n=0
	do i=1,nres1
		if(preali1(i).gt.0) then
			n=n+1
			a(n)=i
		end if
	end do
c	write(*,*) n,' positions aligned',(a(i),i=1,n)
	do i=1,n
		s(i)=0.0
		k=a(i)
		do j=1,n
	l=a(j)
	s(i)=s(i)+scorefun(d1(k,l),d2(preali1(k),preali1(l)))
		end do
c		write(*,*) 'trimming score:',i,k,s(i),preali1(k)
	end do
	do i=1,n
		k=a(i)
		if(s(i).lt.0) then
c			write(*,*) 'trimming',i,k,s(i)
			preali1(k)=0
		end if
	end do
c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine refinepreali(nres1,nres2,preali1,d1,d2,ali1,score)
c
c	ssap routine
c
c	fix (i,j) of preali, fill score table in search space, accumulate traces
c
	implicit none
	include 'gagasizes.for'
	integer nres1,nres2
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres),
     $		ali1(maxres,popsize)
	real score(popsize)
c
	integer left(maxres),rite(maxres),ali(maxres)
	integer table(maxres,maxres),newscore,table0(maxres,maxres)
	integer i,j,k,l,niter
	real scorefun,x,xold,map(maxres,maxres)
	character struc(maxres)
c
	xold=-9.0
	x=0.0
	niter=0
	do while((xold.lt.x).and.(niter.lt.10))
	  niter=niter+1
	  xold=x
	  do i=1,nres1
		do j=1,nres2
			table0(i,j)=0
		end do
	  end do
c
	  call getleftrite(nres1,nres2,preali1,left,rite,10)
	  do k=1,nres1
	    l=preali1(k)
	    if(l.gt.0) then
c
c		fix one pair from prealignment
c
		do i=1,nres1
			do j=1,left(i)
				table(i,j)=0
			end do
			do j=left(i)+1,rite(i)-1
		table(i,j)=nint(100.0*scorefun(d1(k,i),d2(l,j)))
			end do
			do j=rite(i),nres2
				table(i,j)=0
			end do
		end do
		call nw(nres1,nres2,table,newscore,ali)
		write(*,*) 'fixed pair:',k,l,' score:',newscore
		do i=1,nres1
			preali1(i)=ali(i)
		end do
		call trimming(nres1,nres2,preali1,d1,d2)
		call gettotscore(preali1,d1,d2,nres1,x)
		do i=1,nres1
			j=preali1(i)
	if(j.gt.0) table0(i,j)=table0(i,j)+nint(x/100.0)
		end do
	    end if
	  end do
c
c	  consensus alignment
c
	  do i=1,nres1
		do j=1,nres2
			map(i,j)=table0(i,j)
		end do
		struc(i)=' '
	  end do
	  call postplot(map,nres1,nres2,.false.,.true.,struc)
	  call nw(nres1,nres2,table0,newscore,ali)
	  do i=1,nres1
		preali1(i)=ali(i)
	  end do
	  call trimming(nres1,nres2,preali1,d1,d2)
	  do i=1,nres1
		ali1(i,1)=preali1(i)
	  end do
	  call gettotscore(preali1,d1,d2,nres1,x)
	  score(1)=x/100
	  write(*,*) 'iteration:',niter,x,xold
	end do
c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine refinepreali1(nres1,nres2,preali1,d1,d2,ali1,score)
	implicit none
	include 'gagasizes.for'
	integer nres1,nres2
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres),
     $		ali1(maxres,popsize)
	real score(popsize)
c
	integer left(maxres),rite(maxres),ali(maxres)
	integer i,niter,newscore,oldscore,table0(maxres,maxres)
c
	oldscore=-9
	newscore=0
	niter=0
	call trimming(nres1,nres2,preali1,d1,d2)
	do while((newscore.gt.oldscore).and.(niter.lt.10))
		niter=niter+1
		oldscore=newscore
		call getleftrite(nres1,nres2,preali1,left,rite,10)
		call filltable0(table0,left,rite,preali1,nres1,nres2,
     $                  d1,d2)
		call nw(nres1,nres2,table0,newscore,ali)
		do i=1,nres1
			preali1(i)=ali(i)
		end do
		call trimming(nres1,nres2,preali1,d1,d2)
		write(*,*) 'iteration:',niter,newscore,oldscore
		score(1)=float(newscore)/1000
c		write(*,500) (preali1(i),i=1,nres1)
	end do
	do i=1,nres1
		ali1(i,1)=preali1(i)
	end do
c
500	format(20i4)
c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine filltable0(table0,left,rite,preali1,nres1,nres2,
     $          d1,d2)
	implicit none
	include 'gagasizes.for'
	integer nres1,nres2
	integer left(maxres),rite(maxres),table0(maxres,maxres)
	integer*2 d1(maxres,maxres),d2(maxres,maxres),preali1(maxres)
c
	integer i,j,k,l,n,a(maxres)
	real s,scorefun
c
	n=0
	do i=1,nres1
		if(preali1(i).gt.0) then
			n=n+1
			a(n)=i
		end if
	end do
	do i=1,nres1
	  doj=1,left(i)
		table0(i,j)=0
	  end do
	  do j=left(i)+1,rite(i)-1
		s=0
		do k=1,n
			l=a(k)
			s=s+scorefun(d1(i,l),d2(j,preali1(l)))
		end do
		table0(i,j)=max(0,nint(10.0*s))
	  end do
	  do j=rite(i),nres2
		table0(i,j)=0
	  end do
	end do

c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine getleftrite(nres1,nres2,preali1,left,rite,width)
	implicit none
	include 'gagasizes.for'
	integer nres1,nres2,left(maxres),rite(maxres),width
	integer*2 preali1(maxres)
c
	integer i
c
c	search space
c
	do i=1,nres1
		left(i)=0
		rite(i)=nres2
	end do
	do i=2,nres1
		if(preali1(i).gt.0) then
			left(i)=max(left(i-1),preali1(i)-width)
		else
			left(i)=left(i-1)
		end if
	end do
	do i=nres1-1,1,-1
		if(preali1(i).gt.0) then
			rite(i)=min(preali1(i)+width,rite(i+1))
		else
			rite(i)=rite(i+1)
		end if
	end do
c
	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine nw(m,n,table0,score,ali)
	implicit none
	include 'gagasizes.for'
	integer m,n,table0(maxres,maxres),trace(maxres,maxres,2)
	integer table(maxres,maxres),ali(maxres),score
c
	integer i,j,k,x,y
c
	do i=1,m
		ali(i)=0
	end do
c
c	copy input table to work table
c
	score=0
	do i=1,m
		do j=1,n
			trace(i,j,1)=0
			trace(i,j,2)=0
			table(i,j)=table0(i,j)
		end do
		ali(i)=0
	end do
c	write(*,*) ' table input'
c	do i=1,m
c		write(*,500) (table(i,j),j=1,n)
c	end do
c
c	fill trace
c
	do k=2,m
		j=k
		do i=k,m
			call best(m,n,table,i,j,x,y)
			trace(i,j,1)=x
			trace(i,j,2)=y
			table(i,j)=table(i,j)+table(x,y)
		end do
		i=k
		do j=k+1,n
			call best(m,n,table,i,j,x,y)
			trace(i,j,1)=x
			trace(i,j,2)=y
			table(i,j)=table(i,j)+table(x,y)
		end do
	end do
c	write(*,*) ' table after fill'
c	do i=1,m
c		write(*,500) (table(i,j),j=1,n)
c	end do
c
c	find best end-sum
c
	score=-9
	x=0
	y=0
	i=m
	do j=1,n
		if(table(i,j).gt.score) then
			score=table(i,j)
			x=i
			y=j
		end if
	end do
	j=n
	do i=1,m
		if(table(i,j).gt.score) then
			score=table(i,j)
			x=i
			y=j
		end if
	end do
c
c	backtrack alignment
c
c	write(*,*) ' backtrack from ',x,y
	do k=x,1,-1
		ali(x)=y
		i=trace(x,y,1)
		j=trace(x,y,2)
c		write(*,*) i,j,x,y
		x=i
		y=j
		if(x.eq.0) goto 19
		if(y.eq.0) goto 19
	end do
19	continue
c	write(*,*) ' backtrack done'
c
c	a clean-up hack
c
	do i=1,m
		j=ali(i)
		if(j.gt.0) then
			if(table0(i,j).le.0) ali(i)=0
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
	subroutine best(m,n,table,i0,j0,x,y)
	implicit none
	include 'gagasizes.for'
	integer m,n,table(maxres,maxres),i0,j0,x,y
c
	integer s,i,j
c
	x=0
	y=0
	s=-99
	i=i0-1
	do j=1,j0-1
		if(table(i,j).gt.s) then
			s=table(i,j)
			x=i
			y=j
		end if
	end do
	j=j0-1
	do i=1,i0-1
		if(table(i,j).gt.s) then
			s=table(i,j)
			x=i
			y=j
		end if
	end do
c	write(*,*) 'best: ',i0,j0,s,x,y
c
	return
	end
c
c------------------------------------------------------------------------------
c
