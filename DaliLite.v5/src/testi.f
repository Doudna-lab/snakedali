c This module/program is part of DaliLite (c) L. Holm 1999
c
	subroutine testi(nres1,nres2,preali1,d1,d2,ntetra,tetrapool)
	implicit none
	include 'gagasizes.for'
	integer*2 preali1(maxres0),d1(maxres,maxres),d2(maxres,maxres)
	integer nres1,nres2,ntetra
	integer*2 tetrapool(2,maxpair)
c
	integer left(maxres),rite(maxres),a(maxres)
	integer i,j,k,l,n,i0,j0
	real x,x0,x1,x2,scorefun,map(maxres,maxres)
c	character struc(maxres)
c
	write(*,*) 'testi:',nres1,nres2
	n=0
	do i=1,nres1
		if(preali1(i).ne.0) then
			n=n+1
			a(n)=i
		end if
c		struc(i)=' '
	end do
	call getleftrite(nres1,nres2,preali1,left,rite,10)
	do i=1,nres1
		do j=1,nres2
			map(i,j)=0.0
		end do
		do j=left(i)+1,rite(i)-1
			x=0.0
			do k=1,n
				i0=a(k)
				j0=preali1(i0)
				if(i.ne.i0.and.j.ne.j0) x=x+
     $					scorefun(d1(i0,i),d2(j0,j))
			end do
			map(i,j)=x
		end do
	end do
c
c	accept all bands +-3 residues around segments !
c
	do i=1,nres1
		k=preali1(i)
		if(k.gt.0) then
			do j=max(1,k-3),min(nres1-3,k+3)
				map(i,j)=1.0
			end do
		end if
	end do
c
c!	call postplot(map,nres1,nres2,.false.,.true.,struc)
c
c	count tetrapeptides
c
	n=0
	do i=1,nres1-3
		do j=left(i)+1,rite(i)-1
			x=0.0
			do l=0,3
				x=x+map(i+l,j+l)
			end do
			if(x.gt.0.and.j.le.nres2-3.and.n.lt.maxpair) then
				n=n+1
				tetrapool(1,n)=i
				tetrapool(2,n)=j
			end if
		end do
	end do
	ntetra=n
	write(*,*) n,' positive tetrapeptides'

	return
	end


	subroutine postplot(map,nx,ny,lzero,lcolor,struc)
c
c	read a real matrix and plot it a la Conan
c	dimensions given by maxres
c	lzero=.true. => do not plot cells.eq.0.0
c	lcolor=.true. => plot +red/-green 
c
	implicit none
	include 'gagasizes.for'
	integer nx,ny
	real map(maxres,maxres),m,p
	logical lzero,lcolor
	character struc(maxres)
c	
	character*80 line,filnam
	integer i,j
	real x
c
	write(*,*) ' This is PostPlot '
	write(*,*) maxres,nx,ny,lzero
	write(*,*) ' enter name of output file ?'
	read(*,500) filnam
	open(91,file=filnam,status='new')
	open(90,file='/home/holm/conan/prolog.ps',status='old')
10	read(90, 500,end=19) line
	write(91, 500) line
	goto 10
19	close(90)
	write(91,510) filnam
	write(91,*) '/HasLegend true def'
	write(91,*) '/HasAutoScale true def'
	write(91,*) '/HasInfo false def'
	write(91,*) '/XGrid [()] def'
	write(91,*) '/XSequence [(',('X',i=1,nx),')] def'
	write(91,*) '/XDSSP [(',(struc(i),i=1,nx),')]def'
	write(91,*) '/XT0 [()] def'
	write(91,*) '/XT1 [()] def'
	write(91,*) '/XT2 [()] def'
	write(91,*) '/YGrid [()] def'
	write(91,*) '/YSequence [(',('Y',j=1,ny),')] def'
	write(91,*) '/YDSSP [(',(struc(j),j=1,ny),')]def'
	write(91,*) '/YT0 [()] def'
	write(91,*) '/YT1 [()] def'
	write(91,*) '/YT2 [()] def'
	write(91,*) '/FirstResX 1 def'
	write(91,*) '/LastResX ',nx,' def'
	write(91,*) '/FirstResY 1 def'
	write(91,*) '/LastResY ',ny,' def'
	write(91,*) 'InitPage'
c
c	scale black=maximum value,white=minimum value
c
	if(.not.lcolor) then
c
c		black/white
c
		m=0.0
		p=0.0
		do i=1,nx
			do j=1,ny
				if(map(i,j).gt.m) m=map(i,j)
				if(map(i,j).lt.p) p=map(i,j)
			end do
		end do
		write(91,*) '/Black ',m,' def'
		write(91,*) '/White ',p,' def'
		do i=1,nx
			do j=1,ny
				if(lzero.and.(map(i,j).eq.0.0)) goto 129
				write(91,*) i,j,map(i,j),' cs'
129			end do
		end do
	else
c
c		use red/green
c
		m=0.0
		p=0.0
		do i=1,nx
			do j=1,ny
				x=map(i,j)
				if(x.gt.m) m=x
				if(x.lt.p) p=x
			end do
		end do
		write(91,*) '/Black ',p,' def'
		write(91,*) '/White ',0.0,' def'
		do i=1,nx
		  do j=1,ny
			if(lzero.and.(map(i,j).eq.0.0)) goto 109
			if(map(i,j).lt.0.0) write(91,*) i,j,map(i,j),' csg'
109		  end do
		end do
		write(91,*) '/Black ',m,' def'
		write(91,*) '/White ',0.0,' def'
		do i=1,nx
			do j=1,ny
			 if(lzero.and.(map(i,j).eq.0.0)) goto 119
			 if(map(i,j).gt.0.0) write(91,*) i,j,map(i,j),' csr'
119			end do
		end do
	end if
	write(91,*) 'FinishPage'

	close(91)

500	format(a79) 
510	format('/CMName (',a60,') def')
	
	return
	end
