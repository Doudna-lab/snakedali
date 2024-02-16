c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
	function lsave(top1,bestest,stack1,stack)
	implicit none
	include 'parsizes.for'
	integer top1,bestest,stack1(maxstack),stack(maxstack)
	logical lsave
c
	if(top1.le.0) then
		lsave=.false.
	else if(bestest.le.0) then
		lsave=.true.
	else
		lsave=(stack1(top1).ge.stack(bestest))
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine singletstack(top1,stack1,est1,ali1,ali,ns,
     $ seglist,chunk)
	implicit none
	include 'parsizes.for'
	integer top1,stack1(maxstack),est1,ali1(maxseg),ns,
     $ seglist(maxseg)
	integer chunk,ali(maxseg)
c
	integer i,j,ibest,best,x,irray(maxseg*2+1)
c
c	est, ali1(1,...,ns), ali(1,...,ns)
c
c	write(*,*) 'singletstack',top1,est1,ns,chunk
	top1=top1+chunk
	irray(1)=est1
	do i=1,ns
		irray(1+i)=ali1(seglist(i))
		irray(1+ns+i)=ali(seglist(i))
	end do
c
c	check degeneracies
c
c	do i=1,top1,chunk
c		if(stack1(i).eq.est1) then
c			do j=i+1,i+ns
c				if(stack1(j).ne.irray(j-i+1)) goto 19
c			end do
c			write(*,*) 'singletstack degenerate'
c			return
c		end if
c19	end do
c
c	save new
c
	call putchunk(top1,irray,stack1,chunk)
	if(lskipflag) return
c
c	move best to top
c
	ibest=1
	best=stack1(ibest)
	do i=1,top1,chunk
		if(stack1(i).gt.best) then
			ibest=i
			best=stack1(i)
		end if
	end do
	if(ibest.ne.top1) then
		j=top1
		do i=ibest,ibest+chunk-1
			x=stack1(i)
			stack1(i)=stack1(j)
			stack1(j)=x
			j=j+1
		end do
	end if
c
c	write(*,*) 'stack1',top1,(stack1(i),i=1,top1+chunk)
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine perresiduescore(ns,seglist,cutoff,lnul,ali1,upper,
     $ lower,segmentrange,dist,nres1,dist2,nres2,dist2sum,ss,
     $ nseg,idom,domns,domseglist,lseqtl,est,string,ali2,
     $ dist1sum,minseglen,segmentrange0,ldb,outputunit)
	implicit none
	include 'parsizes.for'
	integer ns,seglist(maxseg),ali1(maxseg),nseg,cutoff,idom
	integer domns(maxdom)
	integer domseglist(maxseg,maxdom),segmentrange(2,maxseg),
     $ ali2(maxseg)
	logical lnul,lseqtl,ldb
	integer upper(maxseg,maxseg),lower(maxseg,maxseg),
     $ minseglen(maxseg)
	integer nres1,nres2,dist2sum(nres2,0:nres2),ss(maxseg,maxseg)
	integer*2 dist(nres1*nres1),dist2(nres2*nres2)
	character*10 string
	integer dist1sum(nres1,0:nres1),segmentrange0(2,maxseg)
c
	integer*2 s_beg(-29:maxres2,maxseg),s_end(-29:maxres2,maxseg)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim),outputunit
	real firstest
c
	integer top,ni(maxseg),ci(maxres0,maxseg),start(maxseg,maxseg)
	integer trans(maxres0,maxseg),ali(maxseg)
	logical lexdone(maxseg,maxseg),lhom
	integer laststart,nix,bestest,ixstart(2,maxseg*maxseg)
	integer link(2,maxstack2),ex(exdim),stack(maxstack),
     $ ess(maxseg,maxseg)
	integer is,iseg,chunk,setchunk,est,flag,mi(maxseg)
	integer i,j,ires,preali(maxseg)
c
c	write(*,*) 'perresiduescore',ns,cutoff,lnul,lseqtl,nres2,idom,ns,
c     $		(ali1(seglist(is)),is=1,ns),(ali2(seglist(is)),is=1,ns),
c     $          outputunit
	call initlexdonestart(nseg,ss,start,lexdone)
	laststart=0
	nix=0
	do i=1,maxseg
		do j=1,maxres0
			trans(j,i)=0
		end do
	end do
	do i=1,maxstack2
		link(1,i)=0
		link(2,i)=0
	end do
	do is=1,ns
		iseg=seglist(is)
		i=0
		if(ali1(iseg).ne.nul) then
		  do ires=ali1(iseg),min(ali1(iseg)+9,nres2)
			i=i+1
			ci(i,iseg)=i
			trans(i,iseg)=ires
		  end do
		end if
		if(lnul) then
			i=i+1
			ci(i,iseg)=i
			trans(i,iseg)=nul
		end if
		ni(iseg)=i
		mi(iseg)=i
	end do
c
c	overwrite fixed segments with homolog alignment given in ali2
c
	lhom=.false.
	do is=1,ns
		iseg=seglist(is)
		if(ali2(iseg).ne.nul) then
c			write(*,*) 'fix homolog',iseg,ali2(iseg)
			lhom=.true.
			i=1
			ci(i,iseg)=i
			trans(i,iseg)=ali2(iseg)
			if(lnul) then
				i=i+1
				ci(i,iseg)=i
				trans(i,iseg)=nul
			end if
			ni(iseg)=i
			mi(iseg)=i
		end if
	end do
c
	call update_ex(ni,trans,upper,lower,segmentrange,idom,
     $ domns,domseglist,
     $ ex,dist,nres1,dist2,nres2,dist2sum,nix,laststart,ixstart,1,
     $ start,lexdone,lseqtl,s_beg,s_end,dist1sum,minseglen)
	if(lskipflag) return
c	write(*,*) 'ex',(ex(ir),ir=1,laststart)
c	write(*,*) 'trans'
c	do is=1,ns
c		iseg=seglist(is)
c		write(*,*) iseg,ni(iseg),(trans(i,iseg),i=1,ni(iseg))
c	end do
c	write(*,*) 'ci',((ci(i,seglist(is)),i=1,ni(seglist(is))),is=1,nseg)
	chunk=setchunk(ns,seglist,ni)
	do is=1,ns
		ali1(seglist(is))=nul
	end do
	top=1-chunk
	bestest=0
	call init_link_pointer(link_pointer,forward_pointer,
     $ backward_pointer)
	call get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,trans,lseqtl)
	firstest=float(est)
c	write(*,*) 'firstest singlet',firstest
c
c	-> always separate first segment candidates ! <-
c
	iseg=seglist(1)
	do i=1,mi(iseg)
	  ni(iseg)=1
	  ci(1,iseg)=i
	  call get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,
     $ trans,lseqtl)
	  if(est.gt.cutoff) then
		top=top+chunk
		if(top+chunk.gt.maxstack) return !stop 'perres BUG'
		call putrequest(est,ni,ci,ns,seglist,mi,top,stack,
     $			chunk,bestest,link,link_pointer,
     $			forward_pointer,backward_pointer,firstest)
c		write(*,*) 'singlet',top,bestest,est
	  end if
	end do
	do i=1,maxseg
		preali(i)=nul
	end do
10	call dostack(top,stack,chunk,mi,ns,seglist,cutoff,ex,bestest,
     $		link,trans,start,idom,flag,ali,ali1,98,est,string,.true.,
     $		lseqtl,s_beg,s_end,segmentrange0,link_pointer,
     $		forward_pointer,backward_pointer,firstest,ldb,
     $          outputunit)
c	write(*,*) 'perres returned',flag,est,(ali1(seglist(is)),is=1,ns)
c!!!	if(flag.eq.5.and..not.lcheck(lseqtl,ali1,ns,seglist)) then
c		write(*,*) 'got nonsequential singletali'
c		goto 10
c	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine findhomolog(idom1,ali1,ali_nali1,ali_start1,
     $		ali_save1,ali2,domns,domseglist,lfix)
	implicit none
	include 'parsizes.for'
	integer idom1,ali1(maxseg),ali2(maxseg)
	integer ali_start1(maxdom),ali_nali1(maxdom),ali_save1(100000)
	integer domns(maxdom),domseglist(maxseg,maxdom)
	logical lfix(maxseg,maxdom)
c
	integer is,iseg,ix,iali,x
	integer ix1,ix2
c	ix1=ali_start1(idom1)
c	ix2=ix1+ali_nali1(idom1)*domns(idom1)
c	write(*,*)'findhomol',(ali1(domseglist(is,idom1)),is=1,domns(idom1))
c	write(*,*) idom1,ali_start1(idom1),ali_nali1(idom1),ix1,ix2
c	write(*,*) (ali_save1(ix),ix=ix1,ix2)
	do is=1,domns(idom1)
		  iseg=domseglist(is,idom1)
		  ali2(iseg)=nul
	end do
	do iali=1,ali_nali1(idom1)
		ix=ali_start1(idom1)+(iali-1)*domns(idom1)
		do is=1,domns(idom1)
		  iseg=domseglist(is,idom1)
		  if(lfix(iseg,idom1)) then
			x=ali_save1(ix+is)
c			write(*,*) 'matching',x,ali1(iseg),iseg,is,iali,ix
			if(x.ne.nul) then
				if(x.lt.ali1(iseg)) goto 19
				if(x.gt.ali1(iseg)+9) goto 19
			end if
			ali2(iseg)=x
		  end if
		end do
		return
19	end do
c	write(*,*) 'ali2',(ali2(domseglist(is,idom1)),is=1,domns(idom1))

	return
	end
c
c----------------------------------------------------------------------
c
	function ltransversion(ns,seglist,ali,ali1,alix,trans,mi,ci)
	implicit none
	include 'parsizes.for'
	integer ns,seglist(maxseg),ali(maxseg),ali1(maxseg),mi(maxseg)
	integer trans(maxres0,maxseg),ci(maxres0,maxseg),alix(maxseg)
	logical ltransversion
c
	integer is,iseg,ir
c
	ltransversion=.false.
	do is=1,ns
		iseg=seglist(is)
		alix(iseg)=ali(iseg)
		if((ali1(iseg).eq.nul).and.(trans(ali(iseg),iseg).ne.nul))then
			ltransversion=.true.
			do ir=mi(iseg),1,-1
				if(trans(ci(ir,iseg),iseg).eq.nul) then
					alix(iseg)=ir
					goto 12
				end if
			end do
c			no nul residue among candidate set??
			ir=-infinit
12			continue
c			write(*,*) 'transversion',iseg,ir,alix(iseg),mi(iseg)
		end if
	end do

	return
	end
c
c----------------------------------------------------------------------
c
