c This module/program is part of DaliLite (c) L. Holm 1999

        program filter95fitz
c
c       f77 filter95filtz.f subfitz.f u3b-8.f
c
c       reads fort.95 (output by parsi)
c       writes selection to fort.96 (input for pipe|dalicon)
c
        implicit none
        real zcut1,fitzrcut
        integer fitzmaxiter,inunit,outunit
        character*80 line,dalidatpath_1,dalidatpath_2
c
        if(iargc().lt.7) stop "USAGE: filter95fitz " //
     $          "dat_1 dat_2 zcut1 fitzrcut fitzmaxiter inunit outunit"
        call getarg(1,dalidatpath_1)
        call getarg(2,dalidatpath_2)
        call getarg(3,line)
        read(line,*) zcut1
        call getarg(4,line)
        read(line,*) fitzrcut
        call getarg(5,line)
        read(line,*) fitzmaxiter
        call getarg(6,line)
        read(line,*) inunit
        call getarg(7,line)
        read(line,*) outunit
        call filter95(zcut1,fitzrcut,fitzmaxiter,inunit,outunit,9,
     $          dalidatpath_1,dalidatpath_2)

500     format(a80)

        end
c
c----------------------------------------------------------------------------
c
	subroutine filter95(zcut1,fitzrcut,fitzmaxiter,inunit,outunit,
     $		tmpunit,dalidatpath_1,dalidatpath_2)
	implicit none
        include 'parsizes.for'
	real zcut1
	integer inunit,outunit,tmpunit
c
	character*5 oldcd1,oldcd2,cd1,cd2
	integer idom,score,nseg,ali(4*maxres),i,j,k,l,ndom
        integer node_size(maxdom)
	integer minscore(2,maxdom),nres,nres2,ndom2,node_size2(maxdom)
	real ca(3,maxres),ca2(3,maxres),fitzrcut
	integer fitzmaxiter
	logical lkeep(maxdom)
	character node_type(maxdom),node_type2(maxdom)
	character*80 constructfilnam,filnam,dalidatpath_1,dalidatpath_2
	real mean,sigma,x,rms,u(3,3),t(3)
	integer*2 d1(maxres,maxres),d2(maxres,maxres)
	integer lali,niter,fitzali(maxres),xscore,resix(maxres)
        integer xiser(maxres)
	integer node_nseg(maxdom),node_nseg2(maxdom),nx,oldidom
	integer segmentrange(2,maxseg,maxdom)
        integer segmentrange2(2,maxseg,maxdom)
	real xca(3,maxres),zscore
	logical lnew1,ldebug
	parameter(ldebug=.false.)
	integer tmpali(maxres)
c
	oldcd1='?????'
	oldcd2='?????'
	oldidom=0
	close(inunit)
10	read(inunit,500,end=19,err=10) cd1,cd2,idom,score,nseg,
     $		(ali(i),i=1,nseg*4)
	if(ldebug) write(*,*),'input:', oldcd1,cd1,cd2,idom,score,nseg
        if(cd1.eq.oldcd1.and.cd2.eq.oldcd2.and.idom.eq.oldidom) goto 10
	lnew1=.false.
	if(cd1.ne.oldcd1) then
		filnam=constructfilnam(cd1,dalidatpath_1,'.dat')
		open(tmpunit,file=filnam,status='old',err=11)
                goto 12
                ! try second dalidatpath
11              filnam=constructfilnam(cd1,dalidatpath_2,'.dat')
                open(tmpunit,file=filnam,status='old',err=10)
12              continue
		call readproteindata95(ndom,node_type,node_size,nres,
     $			ca,segmentrange,node_nseg,tmpunit)
		close(tmpunit)
		call getdist95(ca,nres,d1)
		do i=1,ndom
		  lkeep(i)=(node_type(i).eq.'*'.or.node_type(i).eq.'+')
		  if(ldebug) write(*,*) 'idom,lkeep ',i,lkeep(i)
		end do
		!
		! minscore(1,2) is Z-score mean,sigma using query domain size
		!
		do i=1,ndom
		  x=min(float(node_size(i)),400.0)
 		  mean=7.9494+0.70852*x+2.5895e-4*x*x-1.9156e-6*x*x*x
 		  sigma=max(1.0,0.50*mean)
		  minscore(1,i)=nint(10000*mean)
		  minscore(2,i)=nint(10000*sigma)
		end do
		oldcd1=cd1
		lnew1=.true.
	end if
        if(.not.lkeep(idom)) goto 10
	! keep directly if Z-score > 2.0
	zscore=float(score-minscore(1,idom))/minscore(2,idom)
	if(zscore.ge.zcut1) then
	  write(outunit,510) zscore,cd1,cd2,idom,score,nseg,
     $		(ali(i),i=1,nseg*4)
	  if(ldebug)write(*,*),'save 1:',idom,score,zscore,' ',cd1,cd2
          oldcd2=cd2
          oldidom=idom
	  goto 10
	end if
	!
	! give lower scoring pairs a chance using fitz
	!
	if(idom.ne.oldidom.or.lnew1) then
		do i=1,nres
			resix(i)=0
			xiser(i)=0
		end do
		nx=0
		do i=1,node_nseg(idom)
                if(ldebug) write(*,*)
     $                  'domain:',idom,i,(segmentrange(j,i,idom),j=1,2)
			do j=segmentrange(1,i,idom),segmentrange(2,i,idom)
				nx=nx+1
				resix(nx)=j
				xiser(j)=nx
				do k=1,3
					xca(k,nx)=ca(k,j)
				end do
			end do
		end do
		oldidom=idom
		if(ldebug.and.nx.ne.node_size(idom)) write(*,*),
     $			'ERROR: nx.ne.node_size !',cd1,cd2,
     $			idom,nx,node_size(idom),node_nseg(idom),
     $			(segmentrange(1,i,idom),segmentrange(2,i,idom),
     $			i=1,node_nseg(idom))
	end if
	if(cd2.ne.oldcd2) then
		filnam=constructfilnam(cd2,dalidatpath_2,'.dat')
		open(tmpunit,file=filnam,status='old')
		call readproteindata95(ndom2,node_type2,node_size2,
     $		  nres2,ca2,segmentrange2,node_nseg2,tmpunit)
		close(tmpunit)
		call getdist95(ca2,nres2,d2)
		oldcd2=cd2
		!write(*,*),cd1,cd2
	end if
	if(ldebug)write(*,*),'ali',(ali(i),i=1,nseg*4)
	do i=1,nx
		fitzali(i)=0
	end do
	do i=1,nseg
		k=(i-1)*2
		if(ali(k+1).gt.0) then
			do j=ali(k+1),ali(k+2)
				l=xiser(j)
				if(l.gt.0.and.l.le.nres) fitzali(l)=
     $					ali(k+nseg*2+1)+j-ali(k+1)
			end do
		end if
	end do
	if(ldebug)write(*,*),'fitzali1',(fitzali(i),i=1,nx)
	call getut(nx,fitzali,xca,ca2,u,t,lali,rms)
	if(ldebug)write(*,*),'getut done',lali,rms
	call transrotate(xca,nx,u,t)
	if(ldebug)write(*,*),'transrotate done'
	call fitz(xca,nx,ca2,nres2,fitzali,fitzrcut,fitzmaxiter,rms,
     $		lali,niter)
	if(ldebug)write(*,*),'fitzali2',(fitzali(i),i=1,nx)
	do i=1,nres
		tmpali(i)=0
	end do
	do i=1,nx
		if(fitzali(i).gt.0) then
			tmpali(resix(i))=fitzali(i)
		end if
	end do
	call gettotscore95(tmpali,d1,nres,d2,nres2,x)
	if(ldebug)write(*,*),'gettotscore95 done',x
	xscore=nint(10000*x)
	if(ldebug)write(*,*),'fitz done',rms,lali,niter,xscore,score
	if(xscore.gt.score) then	! refine fitzed ali
		nseg=0
		do i=1,nx
			if(fitzali(i).gt.0) nseg=nseg+1
		end do
		j=0
		do i=1,nx
			if(fitzali(i).gt.0) then
				j=j+1
				ali((j-1)*2+1)=resix(i)
				ali((j-1)*2+2)=resix(i)
				ali((j-1)*2+nseg*2+1)=fitzali(i)
				ali((j-1)*2+nseg*2+2)=fitzali(i)
			end if
		end do
		score=xscore
	end if
	zscore=float(score-minscore(1,idom))/minscore(2,idom)
	if(ldebug)write(*,*),'save 2:',idom,score,zscore,' ',cd1,cd2
        if(zscore.lt.zcut1) goto 10 ! reject low-scoring even after fitz
	write(outunit,510) zscore,cd1,cd2,idom,score,nseg,
     $          (ali(i),i=1,nseg*4)
	goto 10
19	close(inunit)
	close(outunit)

500	format(6x,2a5,i4,i20,i4,16000i4)
510	format(f5.1,1x,2a5,i4,i20,i4,16000i4)

	return
	end
c
c----------------------------------------------------------------------------
c
	subroutine getdist95(ca,nres,d1)
	implicit none
	include 'parsizes.for'
	integer nres
	real ca(3,maxres)
	integer*2 d1(nres,nres)
c
	integer i,j
	integer*2 intdistance
c
	do i=1,nres
		d1(i,i)=0.0
		do j=i+1,nres
		  d1(i,j)=intdistance(ca(1,i),ca(2,i),ca(3,i),
     $			 ca(1,j),ca(2,j),ca(3,j))
		  d1(j,i)=d1(i,j)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------------
c
	subroutine readproteindata95(ndom,node_type,node_size,nres,ca,
     $		segmentrange,node_nseg,iunit)
        implicit none
        include 'parsizes.for'
        integer ndom,nres,nseg,iunit,segmentrange(2,maxseg,maxdom)
	integer node_size(maxdom) ,node_nseg(maxdom)
        character node_type(maxdom)
        real ca(3,maxres)
c
        integer i,j,k,idom,iseg
c
        read(iunit,500) nres,nseg
        do iseg=1,nseg
		read(iunit,510) i
        end do
        read(iunit,520) ((ca(j,i),j=1,3),i=1,nres)
        read(iunit,500) ndom
        do idom=1,ndom
        	read(iunit,530) i
        end do
        read(iunit,500) ndom
        do idom=1,ndom
          read(iunit,530) i,node_type(i),j,j,node_size(i),node_nseg(i),
     $	  ((segmentrange(j,k,i),j=1,2),k=1,node_nseg(i))
        end do

500     format(10x,4i5,2x,64000a1)
510     format(6i10)
520     format(10f8.1)
530     format(i4,1x,a1,1x,4i4,64000i4)

        return
        end
c
c-------------------------------------------------------------------------------
c
        subroutine gettotscore95(ali1,d1,nres1,d2,nres2,totscore)
        implicit none
        include 'parsizes.for'
        real totscore
        integer nres1,nres2
        integer ali1(maxres)
c
        real scorefun95
        integer*2 d1(nres1,nres1),d2(nres2,nres2)
c
        integer i,j,k,l,q,r,a(maxres),n
        real x
c
        totscore=0.0
        n=0
        do i=1,nres1
                if(ali1(i).ne.0) then
                        n=n+1
                        a(n)=i
                end if
        end do
        do i=1,n
                k=a(i)
                do j=1,n
                        l=a(j)
                        q=abs(ali1(k))
                        r=abs(ali1(l))
                        x=scorefun95(d1(k,l),d2(q,r))
                        totscore=totscore+x
                end do
        end do
c
500     format(20i4)
c
        return
        end
c
c-------------------------------------------------------------------------------
c
	function scorefun95(r1,r2)
	implicit none
	real scorefun95
	integer*2 r1,r2
	real r,s
c
	r=float(r1+r2)/200.0
	s=(r1-r2)/100.0
	if(r.gt.0.01) then
		scorefun95=(0.20-abs(s)/r)*exp(-r*r/400.0)
	else
		scorefun95=0.20
	end if
c
        return
        end
c
c-------------------------------------------------------------------------------
c
	function intdistance(v1,v2,v3,u1,u2,u3)
	implicit none
	integer*2 intdistance
	real v1,v2,v3,u1,u2,u3,x1,x2,x3
c
	x1=v1-u1
	x2=v2-u2
	x3=v3-u3
	intdistance=nint(100.0*sqrt(x1*x1+x2*x2+x3*x3))
c
        return
        end
c
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
c
