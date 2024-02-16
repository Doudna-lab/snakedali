c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
c	score calculation module
c
c
c	mark zero segment-segment scores with start(iseg,jseg)<0, lex-done
c	use as:
c		iwhere=start(iseg,jseg)
c		if(iwhere.lt.0) then
c			s=0
c		else
c			s=segsegscore()
c		end if
c
c	mi(iseg) is 'ncand(iseg)+1'
c	mi(jseg) is 'ncand(jseg)+1'
c
c----------------------------------------------------------------------
c
	function checkscore(ali1,bl,segmentrange,ns,domseglist,dist,nres1,
     $		dist2,nres2,dist2sum,s_beg,s_end,dist1sum,idom,start)
	implicit none
	include 'parsizes.for'
	integer ali1(maxseg),bl,segmentrange(2,maxseg),checkscore,ns
	integer domseglist(maxseg,maxdom),nres1,nres2,
     $ dist2sum(nres2,0:nres2)
	integer dist1sum(nres1,0:nres1),idom,segsegscore,
     $ start(maxseg,maxseg)
	integer*2 dist(nres1*nres1),dist2(nres2*nres2)
	integer*2 s_beg(-29:maxres2,maxseg),s_end(-29:maxres2,maxseg)
c
	integer is,js,iseg,jseg,x,ires,jres,ss(maxseg,maxseg),
     $ q,p,a1,a2,b1,b2
c
	x=0
	do is=1,ns
	  iseg=domseglist(is,idom)
	  ires=ali1(iseg)
	  if(ires.ne.nul) then
	   do js=1,ns
	    jseg=domseglist(js,idom)
	    jres=ali1(jseg)
	    if(jres.ne.nul.and.start(iseg,jseg).ge.0) then
		a1=segmentrange(1,iseg)
		a2=segmentrange(2,iseg)
		b1=segmentrange(1,jseg)
		b2=segmentrange(2,jseg)
		q=segsegscore(iseg,jseg,ires,jres,a1,a2,b1,b2,
     $		dist,nres1,dist2,nres2,0,0,dist2sum,bl,.true.,s_beg,
     $		s_end,dist1sum)
		p=segsegscore(jseg,iseg,jres,ires,b1,b2,a1,a2,
     $		dist,nres1,dist2,nres2,0,0,dist2sum,bl,.true.,s_beg,
     $		s_end,dist1sum)
		if(p.ne.q) write(*,*) 'asymmetry:',q,p,iseg,jseg,ires,jres
		ss(is,js)=q
		x=x+ss(is,js)
	    else
		ss(is,js)=0
	    end if
	   end do
	  else
	   do js=1,ns
		jseg=domseglist(js,idom)
		ss(is,js)=0
	   end do
	  end if
	  write(*,*) 'ss:',(ss(is,js),js=1,ns),ires,
     $		s_beg(ires,iseg),s_end(ires,iseg)
	end do
	checkscore=x
c
	return
	end
c
c----------------------------------------------------------------------
c
	function checkscore1(ali,ex,start,iseg,jseg,mi)
	implicit none
	include 'parsizes.for'
	integer ali(maxseg),ex(exdim),start(maxseg,maxseg)
	integer checkscore1,iseg,jseg,mi(maxseg)
c
	integer x,aseg,bseg
c
	aseg=min(iseg,jseg)
	bseg=max(iseg,jseg)
	if(start(aseg,bseg).lt.0) then
		x=0
	else if(aseg.eq.bseg) then
		x=ex(start(aseg,aseg)+ali(aseg))
	else
		x=ex(start(aseg,bseg)+(ali(aseg)-1)*mi(bseg)+ali(bseg))
	end if
	checkscore1=x
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine setstart1(iseg,jseg,nix,laststart,ixstart,
     $ mi,start,lexdone)
c
c	singlet(ires,iseg) at ex(start(iseg,iseg)+ires)
c	doublet(ires,jres,iseg,jseg) at ex(iwhere+jres),
c		where iwhere=start(iseg,jseg)+(ires-1)*mi(jseg) and iseg<jseg
c
c	globals: start(maxseg,maxseg),lexdone(maxseg,maxseg)
c
	implicit none
	include 'parsizes.for'
	integer nix,laststart,ixstart(2,maxseg*maxseg),mi(maxseg)
	integer iseg,jseg,start(maxseg,maxseg)
	logical lexdone(maxseg,maxseg)
c
	nix=nix+1
	ixstart(1,nix)=iseg
	ixstart(2,nix)=jseg
	start(iseg,jseg)=laststart
	start(jseg,iseg)=laststart
	if(iseg.eq.jseg) then
		laststart=laststart+mi(iseg)
	else
		laststart=laststart+mi(iseg)*mi(jseg)
	end if
	if(laststart.gt.exdim) then
		write(*,*) 'SEVERE ERROR: setstart overflow',laststart
		!write(97,*) 'SEVERE ERROR: setstart overflow',laststart
		laststart=0
	end if
	lexdone(iseg,jseg)=.true.
	lexdone(jseg,iseg)=.true.
c	write(*,*) 'setstart1',iseg,jseg,start(iseg,jseg),mi(iseg),mi(jseg)

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine update_ex(mi,trans,upper,lower,segmentrange,idom,
     $		domns,domseglist,ex,dist,nres1,dist2,nres2,dist2sum,nix,
     $		laststart,ixstart,bl,start,lexdone,lseqtl,s_beg,s_end,
     $		dist1sum,minseglen)
	implicit none
	include 'parsizes.for'
	integer mi(maxseg),trans(maxres0,maxseg),upper(maxseg,maxseg),idom
	integer lower(maxseg,maxseg),segmentrange(2,maxseg),domns(maxdom)
	integer domseglist(maxseg,maxdom),start(maxseg,maxseg)
	integer ex(exdim),nres1,nres2,dist2sum(nres2,0:nres2),bl
	integer*2 dist(nres1*nres1),dist2(nres2*nres2)
	integer nix,laststart,ixstart(2,maxseg*maxseg),
     $ dist1sum(nres1,0:nres1)
	logical lexdone(maxseg,maxseg),lseqtl
	integer*2 s_beg(-29:maxres2,maxseg),s_end(-29:maxres2,maxseg)
	integer minseglen(maxseg)
c
	integer is,iseg,js,jseg,aseg,bseg
c
c	check lexdone(iseg,jseg)
c
c	write(*,*) 'update_ex',idom,nres1,nres2,nix,laststart,lseqtl,bl
	do is=1,domns(idom)
		iseg=domseglist(is,idom)
		if(.not.lexdone(iseg,iseg)) then
			call setstart1(iseg,iseg,nix,laststart,ixstart,
     $				mi,start,lexdone)
			call singletex(iseg,
     $				upper(iseg,iseg),lower(iseg,iseg),
     $				segmentrange(1,iseg),segmentrange(2,iseg),
     $				mi(iseg),trans,dist2sum,nres2,dist2,nres1,dist,
     $				ex,bl,start,s_beg,s_end,dist1sum,
     $				minseglen(iseg))
			if(lskipflag) return
		end if
		do js=1,is-1
			jseg=domseglist(js,idom)
			if(.not.lexdone(iseg,jseg)) then
				aseg=min(iseg,jseg)
				bseg=max(iseg,jseg)
				call setstart1(aseg,bseg,nix,laststart,
     $				  ixstart,mi,start,lexdone)
				call doubletex(aseg,bseg,mi(aseg),mi(bseg),
     $				  upper(aseg,bseg),lower(aseg,bseg),
     $				  segmentrange(1,aseg),segmentrange(2,aseg),
     $				  segmentrange(1,bseg),segmentrange(2,bseg),
     $				  trans,dist2sum,nres2,dist2,nres1,dist,ex,bl,
     $				  start,lseqtl,s_beg,s_end,dist1sum)
			end if
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine singletex(iseg,upp,low,a1,a2,nir,trans,dist2sum,nres2,
     $		dist2,nres1,dist,ex,bl,start,s_beg,s_end,dist1sum,minlen)
	implicit none
	include 'parsizes.for'
	integer iseg,upp,low,a1,a2,trans(maxres0,maxseg),minlen
	integer nres2,dist2sum(nres2,0:nres2),nir,ex(exdim)
	integer nres1,start(maxseg,maxseg),dist1sum(nres1,0:nres1),bl
	integer*2 dist2(nres2*nres2),dist(nres1*nres1)
	integer*2 s_beg(-29:maxres2,maxseg),s_end(-29:maxres2,maxseg)
c
	integer iwhere,ir,transires,x,l,table(100,100),seglen,i,j,p,q
	integer ibeg,iend,t
	integer t1,segsegscore
	logical lcont
c
c	check segment has internal contacts !
c
	iwhere=start(iseg,iseg)
	if(iwhere.lt.0) return
c
c	check table() dimension
c
	seglen=a2-a1+1
	if(seglen.gt.100) then
	write(*,*) 'FATAL ERROR -- seglen overflow',a1,a2,seglen
		lskipflag=.true.
		return
	end if
c
c	do all candidates
c
	do ir=1,nir
		transires=trans(ir,iseg)
c
c		case: nul residue -> score=zero
c
		if(transires.eq.nul) then
			x=0
			goto 99
		end if
c
c		case: 	inspect multiple (10-blocks) candidates
c			remember highest score x of individual scores t
c
		x=-infinit
		do l=0,bl-1
			ibeg=0
			iend=0
c
c			case: beyond C-terminus
c
			if(transires+l.gt.nres2-minlen+1) then
				t=-infinit
				goto 19
			end if
c
c			case: N- or C-terminal gap
c
			if(transires+l.lt.1) ibeg=1-(transires+l)
			if(transires+l+seglen-1.gt.nres2) iend=
     $				seglen-(nres2-(transires+l)+1)
c
c			calculate table of perresiduescores
c
c
			if(ibeg.gt.seglen-minlen) then
			  t=-infinit
			  ibeg=infinit
			  iend=infinit
			else
			  t=0
			  do i=ibeg,seglen-1-iend
				do j=ibeg,seglen-1-iend
c!use fast indices!
				  p=(a1+i-1)*nres1+a1+j
				  q=(transires+l+i-1)*nres2+transires+l+j
				  table(i+1,j+1)=scoretable(dist2(q),dist(p))
c				  t=t+table(i+1,j+1)
				end do
			  end do
c
c			  trim segment ends of destabilizing residues
c
			  lcont=.true.
			  do while(lcont)
				call trimtable(table,ibeg,iend,seglen,lcont,
     $					t,minlen)
			  end do
			end if
c
c			remember individual end-gaps
c
19			if(transires.ne.nul.and.transires+l.le.nres2.and.
     $                          transires+l.ge.-29) then
				s_beg(transires+l,iseg)=ibeg
				s_end(transires+l,iseg)=iend
			end if
c
c			testing
c
c			t1=segsegscore(iseg,iseg,transires+l,transires+l,
c     $				a1,a2,a1,a2,dist,nres1,dist2,nres2,0,0,
c     $				dist2sum,1,.false.,s_beg,s_end,dist1sum)
c			if(t.ne.t1)write(*,*) 'chksing',t,t1,t-t1,
c     $				transires+l,ibeg,iend,iseg,seglen,minlen,l
c			if(transires.lt.0) write(*,*) 'singlet',
c     $				iseg,ir,transires+l,ibeg,iend,x,iwhere+ir
c
c			remember best score of shifted candidates
c
			x=max(t,x)
		end do
c
c		assign candidate-segment score
c
99		ex(iwhere+ir)=x
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine trimtable(table,ibeg,iend,seglen,lcont,totscore,minlen)
	implicit none
	include 'parsizes.for'
	integer table(100,100),ibeg,iend,seglen,totscore,minlen
	logical lcont
c
	integer rowsum(100),i,j
c
c	minimal segment length=4 residues
c
	if(seglen-ibeg-iend.le.minlen) then
		lcont=.false.
		goto 10
	end if
c
c	compute rowsum
c
	do i=1+ibeg,seglen-iend
		rowsum(i)=0
		do j=1+ibeg,seglen-iend
			rowsum(i)=rowsum(i)+table(i,j)
		end do
	end do
c
c	remove most negative end
c
	lcont=.false.
	if(rowsum(1+ibeg).lt.rowsum(seglen-iend)) then
		if(rowsum(1+ibeg).lt.0) then
			ibeg=ibeg+1
			lcont=.true.
		end if
	else
		if(rowsum(seglen-iend).lt.0) then
			iend=iend+1
			lcont=.true.
		end if
	end if
c
c	current total segsegscore
c
10	totscore=0
	do i=1+ibeg,seglen-iend
		do j=1+ibeg,seglen-iend
			totscore=totscore+table(i,j)
		end do
	end do
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine doubletex(aseg,bseg,nir,njr,upp,low,a1,a2,b1,b2,
     $		trans,dist2sum,nres2,dist2,nres1,dist,ex,bl,start,lseqtl,
     $		s_beg,s_end,dist1sum)
	implicit none
	include 'parsizes.for'
	integer nres2,dist2sum(nres2,0:nres2),ex(exdim),nres1,bl
	integer aseg,bseg,upp,low,a1,a2,b1,b2,nir,njr,
     $ dist1sum(nres1,0:nres1)
	integer trans(maxres0,maxseg),start(maxseg,maxseg)
	integer*2 dist2(nres2*nres2),dist(nres1*nres1)
	logical lseqtl
	integer*2 s_beg(-29:maxres2,maxseg),s_end(-29:maxres2,maxseg)
c
	integer segsegscore
c
	integer transires,transjres,ijstart,ir,jr,iwhere,x
c
	ijstart=start(aseg,bseg)
	if(ijstart.ge.0) then
		do ir=1,nir
			iwhere=ijstart+(ir-1)*njr
			transires=trans(ir,aseg)
			if(transires.eq.nul) then
			  do jr=1,njr
				if(iwhere+jr.gt.exdim) then
        write(*,*) 'FATAL ERROR exdim overflow'
        return ! loop
                                end if
				ex(iwhere+jr)=0
			  end do
			else
			  do jr=1,njr
				transjres=trans(jr,bseg)
				if(transjres.eq.nul) then
					x=0
				else
					x=segsegscore(
     $					 aseg,bseg,transires,transjres,
     $					 a1,a2,b1,b2,dist,nres1,dist2,nres2,
     $			 		 upp,low,dist2sum,bl,lseqtl,s_beg,
     $					 s_end,dist1sum)
				end if
				if(iwhere+jr.gt.exdim) then
        write(*,*) 'FATAL ERROR exdim overflow'
        return ! loop
                                end if
				ex(iwhere+jr)=x
c				write(*,*) 'doublet',aseg,bseg,ir,jr,
c     $					transires,transjres,x,iwhere+jr
			  end do
			end if
		end do
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine initlexdonestart(nseg,ss,start,lexdone)
	implicit none
	include 'parsizes.for'
	integer nseg,ss(maxseg,maxseg),start(maxseg,maxseg)
	logical lexdone(maxseg,maxseg)
c
	integer iseg,jseg
	do iseg=1,nseg
		do jseg=1,nseg
			if(ss(iseg,jseg).gt.0) then
				lexdone(iseg,jseg)=.false.
			else
				lexdone(iseg,jseg)=.true.
			end if
			start(iseg,jseg)=-infinit
		end do
c		write(*,*) 'lexdone:',iseg,(lexdone(iseg,jseg),jseg=1,nseg)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	function segsegscore(iseg,jseg,transires,transjres,a1,a2,b1,b2,
     $		dist,nres1,dist2,nres2,upp1,low1,dist2sum,bl,lseqtl,
     $		s_beg,s_end,dist1sum)
	implicit none
	include 'parsizes.for'
	integer iseg,jseg,transires,transjres,a1,a2,b1,b2,nres1,nres2,bl
	integer segsegscore,upp1,low1
	integer dist2sum(nres2,0:nres2),dist1sum(nres1,0:nres1)
	integer*2 dist(nres1*nres1),dist2(nres2*nres2)
	integer*2 s_end(-29:maxres2,maxseg),s_beg(-29:maxres2,maxseg)
	logical lseqtl
c
	integer x,s,a,b,ishift,jshift,k1,k2,j1,j2,d1,d2,ds2
	integer p,q,ibeg,iend,jbeg,jend,ds1,upp,low
	integer jfirst,jlast,ifirst,ilast,e1,e2,c
	logical lself,l0,l1
c
c	!check start,lseqtl in here!
c
	if(transires.eq.nul.or.transjres.eq.nul) then
		s=0
		goto 99
	end if
	s=-infinit
	if(lseqtl.and.iseg.gt.jseg.and.transires.lt.transjres) goto 99
	if(lseqtl.and.iseg.lt.jseg.and.transires.gt.transjres) goto 99
	lself=(iseg.eq.jseg)
	j1=transires-a1
	j2=transjres-b1
	do jshift=0,bl-1
		if(transjres+jshift.gt.nres2) goto 29
		jbeg=s_beg(transjres+jshift,jseg)
		jend=s_end(transjres+jshift,jseg)
		l0=(jbeg.eq.0.and.jend.eq.0)
		k2=j2+jshift
		d1=b1+k2-1+jbeg
		d2=b2+k2-jend
		e1=b1+jbeg-1
		e2=b2-jend
c		check chain ends !
		jfirst=transjres+jshift+jbeg
		jlast=transjres+jshift-jend+b2-b1
		if(jlast.gt.nres2) goto 29
c
		do ishift=0,bl-1
			if(lself.and.(jshift.ne.ishift)) goto 19
			if(transires+ishift.gt.nres2) goto 19
			ibeg=s_beg(transires+ishift,iseg)
			iend=s_end(transires+ishift,iseg)
			l1=(l0.and.ibeg.eq.0.and.jbeg.eq.0)
			k1=j1+ishift
c			check chain ends !
			ifirst=transires+ishift+ibeg
			ilast=transires+ishift-iend+a2-a1
			if(ilast.gt.nres2) goto 19
c
c			disallow segment overlaps
c
			if(iseg.lt.jseg.and.ilast.ge.jfirst) then
				goto 19
			else if(jseg.lt.iseg.and.jlast.ge.ifirst) then
				goto 19
			end if
c
c			check distancesum
c
c			why asymmetries if call with ishift ? ibeg/iend !
c			if(lself) then
c			  call checkdistsumself(iseg,jseg,a1+k1+ibeg,a2+k1-iend,
c     $				d1,d2,dist2sum,nres2,0,ds2)
c			  if(.not.l1)call checkdistsumself(iseg,jseg,
c     $				a1+ibeg,a2-iend,b1+jbeg-1,b2-jend,
c     $				dist1sum,nres1,0,ds1)
c			else
c			  call checkdistsumpair(iseg,jseg,a1+k1+ibeg,a2+k1-iend,
c     $				d1,d2,dist2sum,nres2,0,ds2)
c			  if(.not.l1)call checkdistsumpair(iseg,jseg,
c     $				a1+ibeg,a2-iend,b1+jbeg-1,b2-jend,
c     $				dist1sum,nres1,0,ds1)
c			end if
cinline
			ds2=0
			do c=a1+k1+ibeg,a2+k1-iend
				ds2=ds2+dist2sum(c,d2)-dist2sum(c,d1)
			end do
c  			! calculate low,upp for ibeg..iend/jbeg..jend rectangle !
			if(.not.l1) then
				ds1=0
				do c=a1+ibeg,a2-iend
					ds1=ds1+dist1sum(c,e2)-dist1sum(c,e1)
				end do

				low=0.7*ds1
				upp=1.3*ds1
				if(ds2.le.low.or.ds2.ge.upp) goto 19
			else
				if(ds2.le.low1.or.ds2.ge.upp1) goto 19
			end if
c
c			precise calculation
c
			x=0
			p=(a1+ibeg-1)*nres1
			q=(a1+ibeg+k1-1)*nres2+k2
			do a=a1+ibeg,a2-iend
c
c			q+b=(transires+ishift-1)*nres2+transjres+jshift+b-b1
c
				do b=b1+jbeg,b2-jend
				  x=x+scoretable(dist2(q+b),dist(p+b))
				end do
c!				if(x.lt.-200001) goto 19
				p=p+nres1
				q=q+nres2
			end do
			if(x.gt.s) s=x
19		end do
29	end do
99	segsegscore=s

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine checkdistsumpair(iseg,jseg,a1beg,a2end,d1,d2,
     $ dist2sum,nres2,ishift,ds2)
	implicit none
	include 'parsizes.for'
	integer iseg,jseg,a1beg,a2end,nres2,d1,d2,ishift
	integer dist2sum(nres2,0:nres2),ds2
c
	integer c
c
c	if(a2end.gt.nres2.or.d2.gt.nres2) write(*,*)
c     $		'checkdistsum WARNING',iseg,jseg,a1beg,a2end,d1,d2,ds2
	ds2=0
	do c=a1beg,a2end
		ds2=ds2+dist2sum(c,d2)-dist2sum(c,d1)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine checkdistsumself(iseg,jseg,a1beg,a2end,d1,d2,
     $ dist2sum,nres2,ishift,ds2)
	implicit none
	include 'parsizes.for'
	integer iseg,jseg,a1beg,a2end,nres2,d1,d2,ishift
	integer dist2sum(nres2,0:nres2),ds2
c
	integer c
c
	ds2=0
	do c=a1beg,a2end
		ds2=ds2+dist2sum(c,d2)-dist2sum(c,d1)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,
     $ trans,lseqtl)
	implicit none
	include 'parsizes.for'
	integer ess(maxseg,maxseg),ns,seglist(maxseg),start(maxseg,maxseg)
	integer ni(maxseg),ci(maxres0,maxseg),ex(exdim),est,mi(maxseg)
	integer trans(maxres0,maxseg)
	logical lseqtl
c
	integer is,iseg,e,js,jseg
c
	est=0
	do is=1,ns
		iseg=seglist(is)
		call get_estimate(e,iseg,iseg,ni,ci,ex,start,mi,trans,lseqtl)
c
c		local term added once
c
		ess(iseg,iseg)=e
		est=est+e
		do js=1,is-1
		  jseg=seglist(js)
		  call get_estimate(e,iseg,jseg,ni,ci,ex,start,mi,trans,lseqtl)
c
c		  cross-terms added twice
c
		  e=e+e
		  ess(iseg,jseg)=e
		  ess(jseg,iseg)=e
		  est=est+e
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine get_estimate1(est,iseg,jseg,ni,ci,ex,start,mi)
	implicit none
	include 'parsizes.for'
	integer iseg,jseg,ni(maxseg),ci(maxres0,maxseg)
	integer ex(exdim),est,mi(maxseg),start(maxseg,maxseg)
c
c	output: est=maximum of ex()
c
	integer i,j,x
	integer iwhere,jw,iwhere0,aseg,bseg
c
	if(iseg.lt.jseg) then
		aseg=iseg
		bseg=jseg
	else
		aseg=jseg
		bseg=iseg
	end if
	if(start(iseg,jseg).lt.0) then
c
c	 if lseqtl violations then est=-infinit!
c
	 est=0
	else
	 est=-infinit
	 iwhere0=start(aseg,bseg)
	 if(aseg.eq.bseg) then
	  do i=1,ni(aseg)
		x=ex(iwhere0+ci(i,aseg))
		if(x.gt.est) est=x
	  end do
	 else
	  jw=mi(bseg)
	  do i=1,ni(aseg)
		iwhere=iwhere0+(ci(i,aseg)-1)*jw
		do j=1,ni(bseg)
			x=ex(iwhere+ci(j,bseg))
			if(x.gt.est) est=x
		end do
	  end do
	 end if
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine get_estimate(est,iseg,jseg,ni,ci,ex,start,mi,
     $ trans,lseqtl)
	implicit none
	include 'parsizes.for'
	integer iseg,jseg,ni(maxseg),ci(maxres0,maxseg)
	integer start(maxseg,maxseg)
	integer ex(exdim),est,mi(maxseg),trans(maxres0,maxseg)
	logical lseqtl
c
c	output: est=maximum of ex()
c
	integer i,j,x
	integer iwhere,jw,iwhere0,aseg,bseg,ares,bres
c
	if(iseg.lt.jseg) then
		aseg=iseg
		bseg=jseg
	else
		aseg=jseg
		bseg=iseg
	end if
	if(start(iseg,jseg).lt.0) then
c
c	 if lseqtl violations then est=-infinit!
c
	 est=0
	 if(lseqtl.and.aseg.ne.bseg.and.ni(aseg).ge.1.and.ni(bseg).ge.1)
     $ then
c		hack: know ci is sorted, nul is last
c		case: either is nul -> return zero
		ares=trans(ci(ni(aseg),aseg),aseg)
		if(ares.eq.nul) return
		bres=trans(ci(ni(bseg),bseg),bseg)
		if(bres.eq.nul) return
c		case: neither is nul, check order
		ares=trans(ci(1,aseg),aseg)
		if(ares.gt.bres) est=-infinit
	 end if
	else
c	 callcount=callcount+ni(aseg)*ni(bseg)
c	 rankcount=rankcount+mi(aseg)*mi(bseg)
c	 testcount=testcount+1
	 est=-infinit
	 iwhere0=start(aseg,bseg)
	 if(aseg.eq.bseg) then
	  do i=1,ni(aseg)
		x=ex(iwhere0+ci(i,aseg))
		if(x.gt.est) est=x
	  end do
	 else
	  jw=mi(bseg)
	  do i=1,ni(aseg)
		iwhere=iwhere0+(ci(i,aseg)-1)*jw
		do j=1,ni(bseg)
			x=ex(iwhere+ci(j,bseg))
			if(x.gt.est) est=x
		end do
	  end do
	 end if
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine multmat(a,adim1,adim2,b,bdim1,bdim2,c,cdim1,cdim2)
c
c	a kertaa b = c
c
	integer adim1,adim2,bdim1,bdim2,cdim1,cdim2
	real a(adim1,adim2),b(bdim1,bdim2),c(cdim1,cdim2)
c
	integer i,j,k
	real s
	logical lskipflag
c
	! check dimensions are valid
	if(adim1.ne.bdim2.or.cdim1.ne.bdim2.or.cdim2.ne.adim2) then
	write(*,*) 'FATAL ERROR -- multmat: incompatible dimensions'
		lskipflag=.true.
		return
	end if
	do i=1,adim2
		do j=1,bdim1
			s=0.0
			do k=1,adim1
				s=s+a(k,i)*b(j,i)
			end do
			c(j,i)=s
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
