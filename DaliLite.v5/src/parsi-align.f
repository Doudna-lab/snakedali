c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
c	alignment module
c
c----------------------------------------------------------------------
c
	subroutine align(ndom,node_child,nres2,lnul,nseg,cut0,ldom,
     $ domns,domseglist,upper,lower,segmentrange,dist,nres1,dist2,
     $ dist2sum,lfix,lfix1,start,lexdone,lseqtl,ss,string,
     $ segment2,nseg2,checkrange,secstr,secstr2,checkx,dist1sum,
     $ minseglen,ngap,segmentrange0,ldb,lfirstonly,lpreali,outputunit)
	implicit none
	include 'parsizes.for'
	integer ndom,node_child(2,maxdom),nres2,nseg,cut0(maxdom)
	logical lseqtl,lfix1(maxseg,maxdom),ldb,lfirstonly
	integer ngap(maxseg),outputunit
	logical lnul,ldom(0:maxdom),lfix(maxseg,maxdom),
     $ lexdone(maxseg,maxseg)
	integer domns(maxdom),domseglist(maxseg,maxdom),
     $ upper(maxseg,maxseg)
	integer lower(maxseg,maxseg),segmentrange(2,maxseg),nres1
	integer*2 dist(nres1*nres1),dist2(nres2*nres2)
	integer dist2sum(nres2,0:nres2),start(maxseg,maxseg),
     $ ss(maxseg,maxseg)
	character*10 string
	integer segment2(maxres2),nseg2,
     $ checkrange(2,maxseg),checkx(maxseg)
	character secstr(maxseg),secstr2(maxseg)
	integer dist1sum(nres1,0:nres1),minseglen(maxseg)
	integer segmentrange0(2,maxseg)
	logical lpreali
c
	integer idom1,idom2
	integer idom,mi(maxseg),ci0(maxres0,maxseg),
     $ trans(maxres0,maxseg)
	integer ex(exdim),nix,laststart,ixstart(2,maxseg*maxseg)
	integer top,stack(maxstack),ali_laststart
	integer ali_start(maxdom),ali_nali(maxdom),ali_save(100000)
	integer ni(maxseg),trans1(maxres0,maxseg),i,j,chunk1
	integer ir,iseg,chunk,link(2,maxstack2),bestest,ali(maxseg)
	integer is,seglist(maxseg),flag,setchunk,stack1(maxstack)
	integer slist(maxseg),tlist(maxseg),nali1,nali2,
     $ ali1(maxseg),est1,top1
	logical lsave,ldegenerate,ltransversion
	integer nmem,alimem(maxseg,1000),ali_laststart1,ali2(maxseg)
	integer ali_start1(maxdom),ali_nali1(maxdom),ali_save1(100000)
	integer alix(maxseg)
c
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $ backward_pointer(0:boxdim)
	real firstest
c
	integer*2 s_beg(-29:maxres2,maxseg),s_end(-29:maxres2,maxseg)
	integer est10,checkscore,checkscore1,est,
     $ ess(maxseg,maxseg),icyc
	integer cut(maxdom)
c
c	initial search space mi,ci0: all residues + nul
c	work search space ni,ci: ship up from mi(iseg),ci0(ir,iseg)+all nuls
c
c	security
c
c	write(*,*) 'align',ndom,nres2,lnul,nseg,nres1,string,nseg2,
c     $		(secstr(i),i=1,nseg),(secstr2(i),i=1,nseg2),ldb,
c     $          lfirstonly,outputunit
	do idom=1,ndom
		cut(idom)=cut0(idom)
	end do
	laststart=0
	nix=0
	top=0
	bestest=0
	ali_laststart=0
	ali_laststart1=0
	do idom=1,ndom
		ali_start(idom)=0
		ali_start1(idom)=0
		ali_nali(idom)=0
		ali_nali1(idom)=0
	end do
	do i=1,maxseg
		do j=1,maxres0
			trans(j,i)=0
			trans1(j,i)=0
		end do
	end do
c
	call init_searchspace(nseg,nres2,10,mi,ci0,trans,lnul,ngap,
     $ lpreali,string,segmentrange0,cut)
	do i=1,ndom
c	write(*,*) 'cut:',i,cut(i),cut0(i)
	end do
c	write(*,*) 'search space',nseg,nres2
c	do iseg=1,nseg
c		write(*,*) iseg,(trans(ir,iseg),ir=1,mi(iseg))
c
c	end do
	do i=1,maxstack2
		link(1,i)=0
		link(2,i)=0
	end do
	do idom=ndom,1,-1
	  icyc=0
	  if(ldom(idom)) then
		idom1=node_child(1,idom)
		idom2=node_child(2,idom)
c		write(*,*) 'idom',idom,idom1,idom2,ndom,domns(idom)
c		write(*,*) 'old search space',idom
		do is=1,domns(idom)
			seglist(is)=domseglist(is,idom)
			iseg=seglist(is)
c			write(*,*) iseg,(trans(ir,iseg),ir=1,mi(iseg))
		end do
c		!trans1 goes with ni
		call shipup(ni,trans1,idom,ldom,domns,domseglist,
     $			mi,ci0,trans,idom1,idom2,
     $			ali_start,ali_nali,ali_save,10,nres2,lnul)
c		write(*,*) 'shipup search space'
c		do is=1,domns(idom)
c			iseg=seglist(is)
c			write(*,*) iseg,(trans1(ir,iseg),ir=1,ni(iseg))
c		end do
c
c		compress ex() for subnodes -> work with mi,ci0,trans
c
		if(idom1.gt.0) call compress(idom1,ni,trans1,mi,ci0,trans,ex,
     $			nix,laststart,ixstart,domns,domseglist,start)
		if(lskipflag) return
		if(idom2.gt.0) call compress(idom2,ni,trans1,mi,ci0,trans,ex,
     $			nix,laststart,ixstart,domns,domseglist,start)
		if(lskipflag) return
		call update_ex(mi,trans,upper,lower,segmentrange,idom,domns,
     $			domseglist,ex,dist,nres1,dist2,nres2,dist2sum,
     $			nix,laststart,ixstart,10,start,lexdone,lseqtl,s_beg,
     $			s_end,dist1sum,minseglen)
		if(lskipflag) return
		call get_ess(ess,est,domns(idom),seglist,mi,ci0,ex,start,mi,
     $			trans,lseqtl)
		firstest=float(est)
c		write(*,*) 'firstest',firstest
		call init_link_pointer(link_pointer,forward_pointer,
     $			backward_pointer)
c		write(*,*) 'update_ex done',domns(idom),idom
c		write(*,*) 'ex',(ex(ir),ir=1,laststart)
c		write(*,*) 'trans',((trans(ir,iseg),ir=1,mi(iseg)),iseg=1,nseg)
c		write(*,*) 'ci0',((ci0(ir,iseg),ir=1,mi(iseg)),iseg=1,nseg)
		if(domns(idom).eq.1) then
			call output1(domseglist(1,idom),mi(domseglist(1,idom)),
     $				ex,cut(idom),ali_nali,ali_start,ali_laststart,
     $				ali_save,trans,idom,start,ldb)
		else if(domns(idom).eq.2) then
			call output2(domseglist(1,idom),domseglist(2,idom),ex,
     $			  cut(idom),mi,ci0,idom,ali_nali,ali_start,
     $			  ali_laststart,ali_save,trans,start,mi,lseqtl,ldb)
		else
			chunk=setchunk(domns(idom),seglist,mi)
			call shortlist_lfix(idom1,domns(idom1),domseglist,lfix,
     $				slist,nali1)
			call shortlist_lfix(idom2,domns(idom2),domseglist,lfix,
     $				tlist,nali2)
			call loadstack(top,stack,ali_start,ali_nali,ali_save,
     $				mi,ci0,trans,idom,domns,domseglist,ex,
     $				chunk,link,bestest,start,slist,tlist,nali1,
     $				nali2,idom1,idom2,domns(idom),seglist,
     $				cut(idom),link_pointer,forward_pointer,
     $				backward_pointer,firstest,lseqtl,lnul)
			if(lskipflag) return
c			write(*,*) 'loadstack done'
			chunk1=1+domns(idom)*2
			top1=1-chunk1
			nmem=0
			flag=5
			do while(flag.eq.5)
10				call dostack(top,stack,chunk,mi,
     $				  domns(idom),seglist,cut(idom),ex,bestest,
     $				  link,trans,start,idom,flag,ali,ali1,97,
     $				  est1,string,.false.,lseqtl,s_beg,s_end,
     $				  segmentrange0,link_pointer,forward_pointer,
     $				  backward_pointer,firstest,ldb,outputunit)
				if(lskipflag) return
c				write(*,*) 'dostack done'
				if(flag.ne.5) goto 99
				est10=est1
c
c				declump all 10-block homologs
c
	call declump(cut(idom),stack,top,chunk,link,bestest,domns(idom),
     $		seglist,mi,ali,ex,start,est10,link_pointer,
     $		forward_pointer,backward_pointer,firstest,trans,lseqtl)
c!	write(*,*) checkscore(ali1,10,segmentrange,domns(idom),
c     $		domseglist,dist,nres1,dist2,nres2,dist2sum,
c     $		s_beg,s_end,dist1sum,idom,start),est1,'check10',string
c	do i=1,domns(idom)
c		iseg=domseglist(i,idom)
c		write(*,*) 'ex:',(checkscore1(ali,ex,start,iseg,
c     $			domseglist(j,idom),mi),j=1,domns(idom)),ali(iseg)
c	end do
c
				if(ldegenerate(ali,alimem,
     $			 	    nmem,domns(idom),seglist)) then
					goto 10
				else
	call findhomolog(idom1,ali1,ali_nali1,ali_start1,
     $		ali_save1,ali2,domns,domseglist,lfix1)
	call findhomolog(idom2,ali1,ali_nali1,ali_start1,
     $		ali_save1,ali2,domns,domseglist,lfix1)
c
				  call perresiduescore(domns(idom),seglist,
     $				    cut(idom),lnul,ali1,upper,lower,
     $				    segmentrange,dist,nres1,dist2,nres2,
     $				    dist2sum,ss,nseg,idom,domns,domseglist,
     $				    lseqtl,est1,string,ali2,
     $				    dist1sum,minseglen,segmentrange0,ldb,
     $                                  outputunit)
				  if(lskipflag) return
				  if(est1.le.cut(idom)) goto 10
c				  write(*,*) 'perres done'
c				  ali1 returns real residue numbers!
c?				  call singletstack(top1,stack1,est1,ali1,ali,
c     $				    domns(idom),seglist,chunk1)
c				  write(*,*) 'singletstack done',top1,est1
c	call declump(cut(idom),stack,top,chunk,link,bestest,domns(idom),
c     $		seglist,mi,ali,ex,start,est1,link_pointer,
c     $		forward_pointer,backward_pointer,firstest,trans,lseqtl)
				end if
c?				if(lsave(top1,bestest,stack1,stack).or.
c?     $	(nali1.gt.0).or.(nali2.gt.0)) then
c				if(lsave(top1,bestest,stack1,stack)) then
c!	write(*,*) checkscore(ali1,1,segmentrange,domns(idom),
c     $		domseglist,dist,nres1,dist2,nres2,dist2sum,
c     $		s_beg1,s_end1,dist1sum,idom,start),est1,'check1'
c
c	case: remove top1, save its ali,
c		declump stack [translate to proper ci-indices!], stack1
c
c				  overwrite ali with best singlets'
c
c	 write(*,*) 'bests',top1,bestest,stack1(top1),stack(bestest)
c???				  if(top1.le.0) then???
c				  est1=stack1(top1)
c				  do is=1,domns(idom)
c					iseg=seglist(is)
c					ali(iseg)=stack1(top1+domns(idom)+is)
c					ali1(iseg)=stack1(top1+is)
c				  end do
c				  write(*,*) (ali(seglist(is)),is=1,domns(idom))
c	write(*,*) (trans(ali(seglist(is)),seglist(is)),is=1,domns(idom))
c				  overwrite singlet-nul-transversions into 10-ali
				  do is=1,domns(idom)
				    iseg=seglist(is)
				    alix(iseg)=ali(iseg)
				    if(ali1(iseg).eq.nul.and.
     $	trans(ali(iseg),iseg).ne.nul) then
				if(trans(ci0(mi(iseg),iseg),iseg).eq.nul)then
					alix(iseg)=mi(iseg)
				      else
					write(*,*) 'transversion bug'
				      end if
				    end if
				  end do
				  call saveit(idom,ali_nali,ali_start,
     $				    ali_laststart,ali_save,trans,domns(idom),
     $				    seglist,alix)
				  if(lskipflag) return
				  top1=top1-chunk1
c
c				  call declump1(est1,ali1,ali)
c				print out singlet-ali; declump using 10:s
c
	if(ali_nali1(idom).eq.0) ali_start1(idom)=ali_laststart1
	ali_nali1(idom)=ali_nali1(idom)+1
	do is=1,domns(idom)
		iseg=seglist(is)
		ali_laststart1=ali_laststart1+1
		ali_save1(ali_laststart1)=ali1(iseg)
	end do
c				  do ir=1,top1,chunk1
c				    do is=1,domns(idom)
c				      iseg=seglist(is)
c		     if(ali(iseg).eq.stack1(ir+domns(idom)+iseg)) then
c					stack1(ir)=-infinit
c					goto 19
c				      end if
c				    end do
c19				  end do
c				  i=1-chunk1
c				  do ir=1,top1,chunk1
c					if(stack1(ir).gt.0) then
c						i=i+chunk1
c						do j=i,i+chunk1-1
c							stack1(j)=stack1(ir+j-i)
c						end do
c					end if
c				  end do
c				  top1=i
c
				  if(est1.gt.0) call segmentpairs(ali1,seglist,
     $				  	domns(idom),idom,segment2,nseg2,
     $					checkrange,string,est1,secstr,secstr2,
     $					nres2,checkx)

				  if(lfirstonly) then
				  	if(idom.eq.1) return
					icyc=icyc+1
				  	if(icyc.ge.5) goto 99
				  end if
c
c?				end if
c
c			go'n'get suboptimal alignment !
			end do
c			write(*,*) 'stack',(stack(ir),ir=1,top+chunk-1)
99		end if
	  end if
	end do
c
c	clean up
c
c	do i=1,ali_laststart1
c		ali_save1(i)=0
c	end do
c	do i=1,ali_laststart
c		ali_save(i)=0
c	end do
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine stackdump(stack,top,chunk,link)
	implicit none
	include 'parsizes.for'
	integer stack(maxstack),top,chunk,link(2,maxstack2)
c
	integer i,j
c
c	write(*,*) 'programmed stack dump',top,chunk
c	do i=1,max(13,top),chunk
c		write(*,*) i,(link(j,i),j=1,2),(stack(i+j),j=0,chunk-1)
c	end do

	return
	end
c
c----------------------------------------------------------------------
c
	function lcheck(lseqtl,ali1,ns,seglist)
	implicit none
	include 'parsizes.for'
	logical lcheck,lseqtl
	integer ali1(maxseg),ns,seglist(maxseg)
c
	integer is,top,x
c
	lcheck=.true.
	top=ali1(seglist(1))
	if(lseqtl) then
		do is=2,ns
			x=ali1(seglist(is))
			if(x.ne.nul) then
				if(x.lt.top) then
					lcheck=.false.
					return
				end if
			end if
			if(x.gt.top) top=x
		end do
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	function ldegenerate(ali,alimem,nmem,ns,seglist)
	implicit none
	include 'parsizes.for'
	logical ldegenerate
	integer ali(maxseg),alimem(maxseg,1000),nmem,ns,seglist(maxseg)
c
	integer i,j,iseg
c
c	write(*,*) 'ldegenerate',(ali(seglist(i)),i=1,ns)
c	do j=1,nmem
c		write(*,*) j,(alimem(seglist(i),j),i=1,ns)
c	end do
	ldegenerate=.false.
c
c	check if ali is already in alimem; remember new !
c
	do i=1,nmem
		do j=1,ns
			iseg=seglist(j)
			if(ali(iseg).eq.nul) goto 19
			if(ali(iseg).ne.alimem(iseg,i)) goto 19
		end do
		ldegenerate=.true.
		goto 99
19	end do
	if(nmem.lt.1000) then
		nmem=nmem+1
		do j=1,ns
			iseg=seglist(j)
			alimem(iseg,nmem)=ali(iseg)
		end do
	end if
99	continue
c	write(*,*) ldegenerate,i,nmem

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine loadstack(top,stack,ali_start,ali_nali,ali_save,mi,ci0,
     $		trans,idom,domns,domseglist,ex,chunk,link,bestest,start,
     $		slist,tlist,nali1,nali2,idom1,idom2,ns,seglist,cutoff,
     $		link_pointer,forward_pointer,backward_pointer,firstest,
     $		lseqtl,lnul)
	implicit none
	include 'parsizes.for'
	logical lnul,lseqtl
	integer idom1,idom2
	integer top,stack(maxstack),mi(maxseg),
     $ ci0(maxres0,maxseg),idom,chunk
	integer trans(maxres0,maxseg),domns(maxdom),
     $ domseglist(maxseg,maxdom)
	integer ali_start(maxdom),ali_nali(maxdom),ali_save(100000),
     $ ex(exdim)
	integer link(2,maxstack2),bestest,slist(maxseg),tlist(maxseg),
     $ ns,cutoff
	integer nali1,nali2,start(maxseg,maxseg),seglist(maxseg)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer ess(maxseg,maxseg),est,is,ia1,ia2,iali,jali,iseg,ir,k,l
	integer ni(maxseg),ci(maxres0,maxseg)
c
c	load stack/fixed segments; include nul with all fixed segments
c
	top=1-chunk
	bestest=0
c
	call get_ess(ess,est,ns,seglist,mi,ci0,ex,start,mi,trans,lseqtl)
c	call putrequest(est,mi,ci0,ns,seglist,mi,top,stack,chunk,
c     $		bestest,link)
c
c	copy mi,ci0 to work arrays ni,ci
c
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
		ni(iseg)=mi(iseg)
		do ir=1,mi(iseg)
			ci(ir,iseg)=ci0(ir,iseg)
		end do
	end do
c
c	overwrite fixed segments in ni,ci with iali*jali
c
	ia1=1
	ia2=1
	if(nali1.gt.1) ia1=ali_nali(idom1)
	if(nali2.gt.1) ia2=ali_nali(idom2)
c
	if(chunk*ia1*ia2.lt.maxstack) then
	  do iali=1,ia1
		if(nali1.gt.1) then
		 do k=1,nali1
		  iseg=domseglist(slist(k),idom1)
		  ni(iseg)=0
		  l=ali_save(ali_start(idom1)+(iali-1)*domns(idom1)+slist(k))
		  if(l.ne.nul) then
			call selectcand(l,iseg,trans,ni,ci,mi)
			if(lskipflag) return
		  end if
		  if(lnul) then
			call selectcand(nul,iseg,trans,ni,ci,mi)
			if(lskipflag) return
		  end if
		 end do
		end if
		do jali=1,ia2
		  if(nali2.gt.1) then
		   do k=1,nali2
		    iseg=domseglist(tlist(k),idom2)
		    ni(iseg)=0
		    l=ali_save(ali_start(idom2)+(jali-1)*domns(idom2)+tlist(k))
		    if(l.ne.nul) then
			call selectcand(l,iseg,trans,ni,ci,mi)
			if(lskipflag) return
		    end if
		    if(lnul) then
			call selectcand(nul,iseg,trans,ni,ci,mi)
			if(lskipflag) return
		    end if
		   end do
		  end if
c		  write(*,*) 'loadstack ci',iali,jali,nali1,nali2
c		  do is=1,domns(idom)
c			iseg=domseglist(is,idom)
c			write(*,*) iseg,(trans(ci(ir,iseg),iseg),ir=1,ni(iseg))
c		  end do
		  call get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,
     $ trans,lseqtl)
		  if(est.gt.cutoff) then
			top=top+chunk
			call putrequest(est,ni,ci,ns,seglist,mi,top,stack,
     $				chunk,bestest,link,link_pointer,
     $				forward_pointer,backward_pointer,firstest)
			if(lskipflag) return
		  end if
19		end do
29	  end do
	else if(chunk*ia1.lt.maxstack) then
	  do iali=1,ia1
		if(nali1.gt.1) then
		 do k=1,nali1
		  iseg=domseglist(slist(k),idom1)
		  ni(iseg)=0
		  l=ali_save(ali_start(idom1)+(iali-1)*domns(idom1)+slist(k))
		  if(l.ne.nul) then
			call selectcand(l,iseg,trans,ni,ci,mi)
			if(lskipflag) return
		  end if
		  if(lnul) then
			call selectcand(nul,iseg,trans,ni,ci,mi)
			if(lskipflag) return
		  end if
		 end do
		end if
		call get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,trans,lseqtl)
		if(est.gt.cutoff) then
			top=top+chunk
			call putrequest(est,ni,ci,ns,seglist,mi,top,stack,
     $				chunk,bestest,link,link_pointer,
     $				forward_pointer,backward_pointer,firstest)
			if(lskipflag) return
		end if
	  end do
	else
	  call get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,trans,lseqtl)
	  if(est.gt.cutoff) then
	    top=top+chunk
	    call putrequest(est,mi,ci0,ns,seglist,mi,top,stack,
     $ chunk,bestest,
     $ link,link_pointer,forward_pointer,backward_pointer,firstest)
	    if(lskipflag) return
	  end if
	end if
c
c	write(*,*) 'loadstack',idom,top,chunk,bestest,nali1,nali2,ia1,ia2

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine selectcand(l,iseg,trans,ni,ci,mi)
	implicit none
	include 'parsizes.for'
	integer l,iseg,trans(maxres0,maxseg),ni(maxseg),ci(maxres0,maxseg)
	integer mi(maxseg)
c
	integer ir,i
c
c	write(*,*) 'selectcand',l,iseg,ni(iseg),mi(iseg)
c	write(*,*) (trans(ir,iseg),ir=1,mi(iseg))
	ni(iseg)=ni(iseg)+1
	ir=1
	do while(trans(ir,iseg).ne.l)
		ir=ir+1
		if(ir.gt.mi(iseg)) then
			write(*,*) l,iseg,(trans(i,iseg),i=1,mi(iseg))
	write(*,*)  'FATAL ERROR -- selectcand ir overflow'
			lskipflag=.true.
			return
		end if
	end do
	ci(ni(iseg),iseg)=ir

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine shortlist_lfix(idom,ns,domseglist,lfix,slist,k)
	implicit none
	include 'parsizes.for'
	logical lfix(maxseg,maxdom)
	integer k,idom,domseglist(maxseg,maxdom),ns,slist(maxseg)
c
	integer is
c
	k=0
	do is=1,ns
		if(lfix(domseglist(is,idom),idom)) then
			k=k+1
			slist(k)=is
		end if
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	function setchunk(ns,seglist,mi)
	implicit none
	include 'parsizes.for'
	integer ns,seglist(maxseg),mi(maxseg),setchunk
c
	integer i,is
c
c	chunk=packed size
c
	i=0
	do is=1,ns
		i=i+mi(seglist(is))
	end do
	setchunk=1+(i+30)/31
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine dostack(top,stack,chunk,mi,ns,seglist,scorecutoff,ex,
     $ bestest,link,trans,start,idom,flag,ali,ali1,iunit,est,
     $ string, lprint,lseqtl,s_beg,s_end,segmentrange,link_pointer,
     $ forward_pointer,backward_pointer,firstest,ldb,outputunit)
	implicit none
	include 'parsizes.for'
	integer top,stack(maxstack),chunk,ns,mi(maxseg),ex(exdim)
	integer ali(maxseg),outputunit
	integer seglist(maxseg),trans(maxres0,maxseg),start(maxseg,maxseg)
	integer est,iunit,segmentrange(2,maxseg)
	integer scorecutoff,link(2,maxstack2),bestest,idom,
     $ flag,ali1(maxseg)
	integer*2 s_beg(-29:maxres2,maxseg),s_end(-29:maxres2,maxseg)
	character*10 string
	logical lprint,lseqtl,lcheck,ldb
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
c	loop calling nextbest
c	cases - perresiduescore, output, clearstack, declump, return
c
c	returns: ali=candidate-numbers, ali1=real residue numbers
c
	integer iseg,is,i,j,k
	integer a1(maxseg),a2(maxseg),b1(maxseg),b2(maxseg),ires,ibeg,iend
c
c	write(*,*) 'dostack',scorecutoff,top,chunk,ns,outputunit,iunit
	do is=1,ns
		ali(seglist(is))=0
	end do
10	call getnextbest(scorecutoff,stack,top,chunk,ns,
     $		seglist,mi,bestest,ali,flag,link,ex,start,est,link_pointer,
     $		forward_pointer,backward_pointer,firstest,trans,lseqtl)
	if(lskipflag) return
	if(flag.eq.1) then
c		write(*,*) 'stack empty'
	else if(flag.eq.2) then
c		write(*,*) 'top score below cutoff'
	else if(flag.eq.3) then
c		write(*,*) 'stack overflow: delete bottom 90 %',top
c		write(iunit,*) 'stack overflow: delete bottom 90 %',top
		if(top.eq.1) then
                        write(*,*) 'BUG in dostack'
			flag=4
			return
		end if
c		delete 90 %
		j=bestest
		do i=1,top/chunk/10
			k=link(1,j)
			j=k
		end do
c		write(*,*) 'call clear',bestest,stack(bestest),j,stack(j),top
		k=stack(j)
		call clearstack(k,stack,top,chunk,link,bestest,
     $		  link_pointer,forward_pointer,backward_pointer,firstest)
		if(lskipflag) return
		flag=0
	else if(flag.eq.4) then
		write(*,*) 'no candidates'
	else if(flag.eq.5) then
c		return real residue numbers in ali
		do is=1,ns
			iseg=seglist(is)
			ali1(iseg)=trans(ali(iseg),iseg)
		end do
		if(.not.ldb)write(*,500) string,idom,est,
     $			(seglist(is),is=1,ns),(ali1(seglist(is)),is=1,ns)
c!!		if(lprint.and.lcheck(lseqtl,ali1,ns,seglist)) then
		if(lprint) then
c			  write(iunit,500) string,idom,est,(seglist(is),
c     $				is=1,ns),(ali1(seglist(is)),is=1,ns)
c
c			write out aligned residue ranges
c
			do is=1,ns
				iseg=seglist(is)
				ires=ali1(iseg)
				if(ires.ne.nul) then
					ibeg=s_beg(ires,iseg)
					iend=s_end(ires,iseg)
					a1(is)=segmentrange(1,iseg)+ibeg
					a2(is)=segmentrange(2,iseg)-iend
					b1(is)=ires+ibeg
					b2(is)=b1(is)+a2(is)-a1(is)
				else
					a1(is)=nul
					a2(is)=nul
					b1(is)=nul
					b2(is)=nul
				end if
			end do
!                        write(outputunit,510) string,idom,est,ns,
!     $          (a1(is),a2(is),is=1,ns),(b1(is),b2(is),is=1,ns)
			write(outputunit,*) 'refine'//string,idom,est,
     $		ns,(a1(is),a2(is),is=1,ns),(b1(is),b2(is),is=1,ns)
		end if
	end if
	if(flag.eq.0) goto 10
c
500	format('unique',a10,i4,i20,12i4,8(/8x,18i4))
510	format('refine',a10,i4,i20,i4,800i4)
c The above line has changed from 320i4 to 800i4 CW
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine output2(aseg,bseg,ex,cut,ni,ci,idom,ali_nali,
     $ ali_start,ali_laststart,ali_save,trans,start,mi,lseqtl,ldb)
	implicit none
	include 'parsizes.for'
	logical lseqtl,ldb
	integer aseg,bseg
	integer ex(exdim),cut,ni(maxseg),ci(maxres0,maxseg),mi(maxseg)
	integer idom,ali_nali(maxdom),ali_start(maxdom),ali_laststart
	integer ali_save(100000),trans(maxres0,maxseg),
     $ start(maxseg,maxseg)
c
	integer ires,jres,iseg,jseg
	integer iistart,jjstart,ijstart,ir,jr
	integer iwhere,jw,n,x,ali(maxseg),seglist(maxseg)
c
	iseg=min(aseg,bseg)
	jseg=max(aseg,bseg)
	iistart=start(iseg,iseg)
	jjstart=start(jseg,jseg)
	ijstart=start(iseg,jseg)
	seglist(1)=aseg
	seglist(2)=bseg
	n=0
	jw=mi(jseg)
c	write(*,*) 'output2',aseg,bseg,cut,iseg,jseg,iistart,ijstart,jjstart,jw
	do ir=1,ni(iseg)
	  ires=ci(ir,iseg)
	  iwhere=ijstart+(ires-1)*jw
	  do jr=1,ni(jseg)
		jres=ci(jr,jseg)
		if(lseqtl.and.(trans(jres,jseg).lt.trans(ires,iseg))) goto 19
		x=ex(iistart+ires)+ex(jjstart+jres)
		if(iwhere.ge.0) x=x+ex(iwhere+jres)
		if(x.gt.cut) then
		    ali(iseg)=ires
		    ali(jseg)=jres
		    call saveit(idom,ali_nali,ali_start,ali_laststart,
     $			ali_save,trans,2,seglist,ali)
		    if(lskipflag) return
		    n=n+1
c		    write(*,*) 'saved',aseg,bseg,ires,jres,trans(ires,iseg),
c     $			trans(jres,jseg),x,n
		end if
c		if(trans(ires,iseg).lt.0.or.trans(jres,jseg).lt.0)
c     $			write(*,*) 'negat',aseg,bseg,ires,jres,trans(ires,iseg),
c     $			trans(jres,jseg),x
19	  end do
	end do
	if(.not.ldb)write(*,*) n,' doublets above',cut

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine output1(iseg,nir,ex,cut,ali_nali,ali_start,
     $		ali_laststart,ali_save,trans,idom,start,ldb)
	implicit none
	include 'parsizes.for'
	integer iseg,nir,cut,ex(exdim),idom
	integer ali_nali(maxdom),ali_start(maxdom),ali_laststart
	integer ali_save(100000),trans(maxres0,maxseg),
     $ start(maxseg,maxseg)
	logical ldb
c
	integer n,ir,ali(maxseg),seglist(maxseg)
c
c	write(*,*) 'output1',iseg,nir,cut,idom,ali_laststart
	n=0
	seglist(1)=iseg
	do ir=1,nir
		if (ex(start(iseg,iseg)+ir).gt.cut) then
			n=n+1
			ali(iseg)=ir
			call saveit(idom,ali_nali,ali_start,ali_laststart,
     $				ali_save,trans,1,seglist,ali)
			if(lskipflag) return
c!	write(*,*) 'saved',iseg,ir,ex(start(iseg,iseg)+ir),trans(ir,iseg)
		end if
	end do
	if(.not.ldb)write(*,*) n,' singlets above',cut

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine saveit(idom,ali_nali,ali_start,ali_laststart,ali_save,
     $		trans,ns,seglist,ali)
	implicit none
	include 'parsizes.for'
	integer idom,ali_nali(maxdom),ali_start(maxdom),ali_laststart
	integer ali_save(100000),trans(maxres0,maxseg),
     $ ns,seglist(maxseg)
	integer ali(maxseg)
c
	integer is,iseg,ix,i
c
c	check degeneracy
c
	ix=ali_start(idom)
	do i=1,ali_nali(idom)
		do is=1,ns
			iseg=seglist(is)
			if(trans(ali(iseg),iseg).ne.ali_save(ix+is)) goto 19
		end do
		write(*,*) 'saveit degenerate'
		return
19		ix=ix+ns
	end do
c
c
c
	if(ali_nali(idom).eq.0) ali_start(idom)=ali_laststart
	ali_nali(idom)=ali_nali(idom)+1
	do is=1,ns
		iseg=seglist(is)
		ali_laststart=ali_laststart+1
		if(ali_laststart.gt.100000) then
	write(*,*)  'FATAL ERROR -- ali_laststart overflow'
			lskipflag=.true.
			return
		end if
		ali_save(ali_laststart)=trans(ali(iseg),iseg)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine compress(idom,ni,trans1,mi,ci0,trans,ex,nix,laststart,
     $		ixstart,domns,domseglist,start)
	implicit none
	include 'parsizes.for'
	integer ex(exdim),nix,laststart,ixstart(2,maxseg*maxseg),idom
	integer mi(maxseg),ci0(maxres0,maxseg),trans(maxres0,maxseg)
	integer ni(maxseg),trans1(maxres0,maxseg)
	integer domns(maxdom),domseglist(maxseg,maxdom),
     $ start(maxseg,maxseg)
c
	integer iseg,ires,ir,irold(maxres0,maxseg),is,iwhere,iwhereold
	integer aseg,bseg,jseg,jr,js,ix,oldstart,i
c
c	write(*,*) 'compress',idom,nix,laststart,
c     $		(domseglist(is,idom),is=1,domns(idom))
c
c	make table current ir (ni,ci) <- old ir (mi,ci0)
c
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
c		write(*,*) 'iseg is',iseg
		do ir=1,ni(iseg)
			ires=1
			do while(trans(ires,iseg)
     $			  .ne.trans1(ir,iseg))
				ires=ires+1
				if(ires.gt.mi(iseg)) then
c	write(*,*) ir,ires,iseg,mi(iseg),ni(iseg),trans1(ir,iseg)
c	write(*,*) (trans1(i,iseg),i=1,ni(iseg))
c	write(*,*) (trans(i,iseg),i=1,mi(iseg))
	write(*,*)  'FATAL ERROR compress'
					lskipflag=.true.
					return
				end if
			end do
			irold(ir,iseg)=ires
		end do
c		write(*,*) 'irold',iseg,(irold(ir,iseg),ir=1,ni(iseg))
	end do
c
c	ex(start(iseg,iseg)+ir) <- ex(start(iseg,iseg)+irold(ir))
c	ex(start(iseg,jseg)+(ir-1)*ni(jseg)+jr)
c		<- ex(start(iseg,jseg)+(irold(ir)-1)*mi(jseg)+irold(jr))
c
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
		iwhere=start(iseg,iseg)
c		write(*,*) 'ex',is,iseg,iwhere,ni(iseg),mi(iseg)
		do ir=1,ni(iseg)
c			write(*,*) iseg,ir,irold(ir,iseg),trans1(ir,iseg),
c     $				ex(iwhere+irold(ir,iseg))
			ex(iwhere+ir)=ex(iwhere+irold(ir,iseg))
		end do
		do js=1,is-1
			jseg=domseglist(js,idom)
			aseg=min(iseg,jseg)
			bseg=max(iseg,jseg)
c			write(*,*) 'jsloop',is,js,jseg,aseg,bseg,
c     $	start(aseg,bseg),ni(bseg),mi(bseg),ni(aseg),mi(aseg)
c
c			skip segment-segment pairs that have no contacts, ex()
c
			if(start(aseg,bseg).ge.0) then
			  do ir=1,ni(aseg)
				iwhere=start(aseg,bseg)+(ir-1)*ni(bseg)
				iwhereold=start(aseg,bseg)+
     $					(irold(ir,aseg)-1)*mi(bseg)
c	if(ir.eq.1)
c	 write(*,*) 'isjs',is,js,ir,aseg,bseg,iwhere,iwhereold,
c     $		(ex(iwhereold+jr),jr=1,mi(bseg)),ni(bseg),mi(bseg),
c     $		(irold(jr,bseg),jr=1,ni(bseg))
				do jr=1,ni(bseg)
				  ex(iwhere+jr)=ex(iwhereold+irold(jr,bseg))
				end do
			  end do
			end if
		end do
	end do
c
c	mi,ci0,trans(ir) <- ni,ci,trans1(ir)
c
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
		do ir=1,ni(iseg)
			trans(ir,iseg)=trans1(ir,iseg)
			ci0(ir,iseg)=ir
		end do
		mi(iseg)=ni(iseg)
	end do
c
c	! left-shift ex-chunks ! ixstart,start,laststart !
c
	laststart=0
	do ix=1,nix
		iseg=ixstart(1,ix)
		jseg=ixstart(2,ix)
		oldstart=start(iseg,jseg)
		start(iseg,jseg)=laststart
		start(jseg,iseg)=laststart
c	write(*,*) 'leftshift',iseg,jseg,oldstart,laststart,mi(iseg),mi(jseg)
		if(iseg.eq.jseg) then
			i=0
			do ir=1,mi(iseg)
				i=i+1
				ex(laststart+i)=ex(oldstart+i)
			end do
			laststart=laststart+i
		else
			i=0
			do ir=1,mi(iseg)
				do jr=1,mi(jseg)
					i=i+1
					ex(laststart+i)=ex(oldstart+i)
				end do
			end do
			laststart=laststart+i
		end if
	end do
c	write(*,*) 'laststart compressed to',laststart

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine shipup(ni,trans1,idom,ldom,domns,domseglist,mi,
     $ ci0,trans, idom1,idom2,ali_start,ali_nali,ali_save,bl,
     $ nres2,lnul)
	implicit none
	include 'parsizes.for'
	integer ni(maxseg),idom,domns(maxdom)
	integer domseglist(maxseg,maxdom),mi(maxseg),ci0(maxres0,maxseg)
	integer trans(maxres0,maxseg),idom1,idom2,bl,nres2
	integer trans1(maxres0,maxseg)
	integer ali_start(maxdom),ali_nali(maxdom),ali_save(100000)
	logical ldom(0:maxdom),lnul
c
	integer is,ir,iseg
c
c	ldom(idom1)=T: fill ni,trans1 by candidates given in ali_save[idom1]+nuls
c	ldom(idom1)=F: fill ni,trans1 directly from mi,trans
c
	if(idom1.eq.0) then
	  call loadall(idom,ni,trans1,mi,ci0,trans,domns,domseglist)
	else
	  if(ldom(idom1)) then
		call loadali(idom1,domns,domseglist,0,0,nres2,ni,trans1,bl,
     $			ali_start,ali_nali,ali_save)
	  else
		call loadall(idom1,ni,trans1,mi,ci0,trans,domns,domseglist)
	  end if
	  if(ldom(idom2)) then
		call loadali(idom2,domns,domseglist,0,0,nres2,ni,trans1,bl,
     $			ali_start,ali_nali,ali_save)
	  else
		call loadall(idom2,ni,trans1,mi,ci0,trans,domns,domseglist)
	  end if
	end if
c
c	add nuls if missing
c
	if(lnul) then
		do is=1,domns(idom)
			iseg=domseglist(is,idom)
			do ir=1,ni(iseg)
				if(trans1(ir,iseg).eq.nul) then
c					write(*,*) 'nul already in',iseg,ir
					goto 19
				end if
			end do
			ni(iseg)=ni(iseg)+1
			trans1(ni(iseg),iseg)=nul
19		end do
	end if
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine loadali(idom,domns,domseglist,fuzzN,fuzzC,nres2,ni,
     $		trans,bl,ali_start,ali_nali,ali_save)
	implicit none
	include 'parsizes.for'
	integer idom,domns(maxdom),domseglist(maxseg,maxdom)
	integer ni(maxseg),nres2,bl
	integer ali_start(maxdom),ali_nali(maxdom),ali_save(100000)
	integer trans(maxres0,maxseg),fuzzN,fuzzC
c
	integer ires,is,iseg,i,ix,iali,stres
	parameter (stres=-29)
	logical lcand(stres:maxres2,maxseg)
c
c	load prealignment/all candidates
c
c	write(*,*) 'loadali',idom,domns(idom),nres2,fuzzN,fuzzC,
c     $		ali_start(idom),ali_nali(idom),
c     $		(ali_save(ali_start(idom)+i),i=1,ali_nali(idom)*domns(idom))
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
		do ires=stres,nres2
			lcand(ires,iseg)=.false.
		end do
	end do
	ix=ali_start(idom)
	do iali=1,ali_nali(idom)
		do is=1,domns(idom)
			ix=ix+1
			ires=ali_save(ix)
			iseg=domseglist(is,idom)
			if(ires.ne.nul) then
			  do i=max(stres,ires-fuzzN),min(nres2,ires+fuzzC)
				lcand(i,iseg)=.true.
			  end do
			end if
		end do
	end do
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
		ni(iseg)=0
		do ires=stres,nres2,bl
		  if(lcand(ires,iseg)) then
			ni(iseg)=ni(iseg)+1
			trans(ni(iseg),iseg)=ires
		  end if
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine loadall(idom,ni,trans1,mi,ci0,trans,domns,domseglist)
	implicit none
	include 'parsizes.for'
	integer ni(maxseg),idom,domns(maxdom)
	integer domseglist(maxseg,maxdom),mi(maxseg),ci0(maxres0,maxseg)
	integer trans(maxres0,maxseg),trans1(maxres0,maxseg)
c
	integer is,ir,iseg
c
c	write(*,*) 'loadall',idom,domns(idom),((trans(ir,domseglist(is,idom)),
c     $		ir=1,mi(domseglist(is,idom))),is=1,domns(idom))
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
		ni(iseg)=mi(iseg)
		do ir=1,ni(iseg)
			trans1(ir,iseg)=trans(ir,iseg)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine init_searchspace(nseg,nres2,bl,mi,ci0,trans,lnul,ngap,
     $		lpreali,string,segmentrange0,cut)
c
c	initialize search space: all residues + nul + N-terminal gap residues
c
	implicit none
	include 'parsizes.for'
	integer nseg,nres2,bl,mi(maxseg),ci0(maxres0,maxseg),ngap(maxseg)
	integer trans(maxres0,maxseg),segmentrange0(2,maxseg),cut(maxdom)
	logical lnul,lpreali
	character*10 string,filnam
c
	integer iseg,ires,i,j,exclusion(maxseg,maxres0),ndom
c
        character*5 cd1,cd2
        real score,zscore,cutz(maxdom)
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
        integer ierr,k,restoseg(maxres),np,na
c
c
	!write(*,*) 'init_searchspace',nseg,nres2,bl,lnul,lpreali,string
	if(lpreali) then
		! restoseg
		do i=1,maxres
			restoseg(i)=0
		end do
		do iseg=1,nseg
			do i=segmentrange0(1,iseg),segmentrange0(2,iseg)
				restoseg(i)=iseg
			end do
		end do
		! initialize exclusion table
		do j=1,(nres2+bl-1)/bl
			do i=1,nseg
				exclusion(i,j)=0
			end do
		end do
		! open cd1.dccp
		if(string(5:5).eq.'_'.or.string(5:5).eq.' ') then
			filnam='code.cutz '
			filnam(1:4)=string(1:4)
		else
			filnam='codeX.cutz'
			filnam(1:5)=string(1:5)
		end if
		open(89,file=filnam,status='old',err=19)
c       	write(*,*) 'reading prealignment from ',filnam
		! accumulate prealignments in exclusion table
		ierr=0
		np=0
		na=0
		cd1=string(1:5)
		if(cd1(5:5).eq.'_') cd1(5:5)=' '
		cd2=string(6:10)
		if(cd2(5:5).eq.'_') cd2(5:5)=' '
		do while(ierr.eq.0)
			call read_nextcutz(cd1,cd2,score,zscore,
     $         			 nblock,l1,r1,l2,r2,ierr,89,cutz,ndom)
			if(ierr.eq.0) then
c
c				set cut to max(old,prealiscore*80 %)
c
				do i=1,ndom
c  JONG changed the following to avoid 'max' function
					if(8000*cutz(i).gt.cut(i)) then
						cut(i)=8000*cutz(i)
					end if
c                   cut(i)=max(cut(i),8000*cutz(i))
				end do
c
c				declump prealignment
c
				na=na+1
				do i=1,nblock
					do j=l1(i),r1(i)
						k=l2(i)+j-l1(i)
						iseg=restoseg(j)
						ires=(k+bl-1)/bl
						if(iseg.gt.0.and.ires.gt.0)then
						  exclusion(iseg,ires)=
     $							exclusion(iseg,ires)+1
						  np=np+1
c	type *,'excluded',j,k,iseg,ires,exclusion(iseg,ires)
						end if
					end do
				end do
			end if
		end do
		close(89)
19              continue
c               write(*,*) na,' prealignments',np,' equivalences used'
		! load search space minus exclusion table
		! (iseg -> segmentrange) * (ires -> 10-block)
		do iseg=1,nseg
			mi(iseg)=0
			ires=1
			do i=ngap(iseg),1,-1
				ires=ires-i*bl
c				write(*,*) 'add N-gap residue',ires,iseg
				if(ires.gt.0) then
				  j=(ires+bl-1)/bl
				  if(exclusion(iseg,j).lt.6) then
					mi(iseg)=mi(iseg)+1
					ci0(mi(iseg),iseg)=mi(iseg)
					trans(mi(iseg),iseg)=ires
				  end if
				end if
			end do
			do ires=1,nres2,bl
				j=(ires+bl-1)/bl
				if(exclusion(iseg,j).lt.6) then
					mi(iseg)=mi(iseg)+1
					ci0(mi(iseg),iseg)=mi(iseg)
					trans(mi(iseg),iseg)=ires
				end if
			end do
			if(lnul) then
				mi(iseg)=mi(iseg)+1
				ci0(mi(iseg),iseg)=mi(iseg)
				trans(mi(iseg),iseg)=nul
			end if
		end do
	else
	  do iseg=1,nseg
		mi(iseg)=0
		ires=1
		do i=ngap(iseg),1,-1
			ires=ires-i*bl
c			write(*,*) 'add N-terminal gap residue',ires,iseg
			mi(iseg)=mi(iseg)+1
			ci0(mi(iseg),iseg)=mi(iseg)
			trans(mi(iseg),iseg)=ires
		end do
		do ires=1,nres2,bl
			mi(iseg)=mi(iseg)+1
			ci0(mi(iseg),iseg)=mi(iseg)
			trans(mi(iseg),iseg)=ires
		end do
		if(lnul) then
			mi(iseg)=mi(iseg)+1
			ci0(mi(iseg),iseg)=mi(iseg)
			trans(mi(iseg),iseg)=nul
		end if
	  end do
	end if

	return
	end
c
c----------------------------------------------------------------------
c
       subroutine read_nextcutz(cd1,cd2,score,zscore,nblock,l1,r1,l2,r2,
     $          ierr,iunit,cutz,ndom)
        implicit none
        include 'parsizes.for'
c
c       reads dali-scores & Z-scores from domainparser | corescore output !
c
        character*5 cd1,cd2,c1,c2
        character*2 line
        real score,zscore,cutz(maxdom)
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
        integer i,ierr,iunit,ndom
        integer ide,lali
        real rmsd
c
c1bfg 1bfg               1604.8                30.9
c           11   531.2183       355.7889       28.14490       32.13352
c    208.2986       8.151642       8.163538       97.22334       32.09985
c    27.90044       28.12995
c            1           1         126           1         126
c//
c
100	read(iunit,500,err=100,end=999) c1,c2,score,zscore
	if(c1.ne.cd1.or.c2.ne.cd2) then
110		read(iunit,510,end=999) line
		if(line(1:2).ne.'//') goto 110
		goto 100
	end if
	read(iunit,*) ndom,(cutz(i),i=1,ndom)
	read(iunit,*) nblock,(l1(i),r1(i),i=1,nblock),(l2(i),r2(i),
     $ i=1,nblock)
500	format(2a5,2f20.1)
510	format(a2)

!       normal exit
        ierr=0
        return

!       error exit
999     ierr=-1
        return
        end
c
c------------------------------------------------------------------------------
c
