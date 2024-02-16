c This module/program is part of DaliLite (c) L. Holm 1999
c
c
c----------------------------------------------------------------------
c
c	stack module
c		getnextbest: return best unique alignment or exception flag
c		clearstack: after overflow, make space by removing baddies
c		swapstack: chunk1 <-> chunk2
c		putchunk, putstack: pack, save chunk
c		getchunk, getstack: read, unpack chunk
c		declump: blank out alignment trace from stack
c		split: divide search space
c
c	variables
c		chunk: 1+[sum(is=1,ns)mi(seglist(is))+30]/31; est+candidates
c		link(): linked list of score estimates, points to stack
c			link(1,*) lower, link(2,*) higher score
c		bestest: points to beginning of linked list
c		stack: in chunks
c		top: starts from 1, points to estimate of chunk
c		ex(): energy table
c		start(iseg,jseg): pointer to ex(), index<0 means score is zero
c		ns: number of segments
c		mi: maximal number of residues in search space, including nul
c		ni,ci: work search space [packed]
c		ci0: initial search space [packed]
c
c----------------------------------------------------------------------
c
	subroutine split(est,ns,seglist,ex,top,stack,scorecutoff,ni,ci,mi,
     $		chunk,bestest,link,start,link_pointer,forward_pointer,
     $		backward_pointer,firstest,trans,lseqtl)
	implicit none
	include 'parsizes.for'
	integer re_estimate,start(maxseg,maxseg),trans(maxres0,maxseg)
	integer ns,seglist(maxseg),mi(maxseg),chunk,bestest,
     $ link(2,maxstack2)
	integer est,ex(exdim),scorecutoff,top,stack(maxstack)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
	logical lseqtl
c
c	split xseg according to ex(xseg,yseg,...)
c		find emax(xres)= max of ex(xseg,yseg,xres,...)
c		erange=smallest to largest emax(...)
c		lcand1(xres)=lcand(xseg,xres).and.emax(xres).le.midpoint
c		lcand2(xres)=lcand(xseg,xres).and.emax(xres).gt.midpoint
c		special case: smallest.eq.largest => assign every second to 1,2
c
	integer mini,maxi,xres,e,emax(maxres0),midpoint,
     $ ess(maxseg,maxseg)
	integer f,n1,ni(maxseg),ci(maxres0,maxseg),xr,yr
	integer iwhere,jw,ijstart
c
	integer emaxim,is,js,iseg,jseg,xseg,yseg,seed
	integer ni1,ni2,ci1(maxres0),ci2(maxres0),ci3(maxres0),
     $ ci4(maxres0)
	logical lclosed(maxseg)
	integer bestxres,ni3,ni4,ni5,ci5(maxres0),pr,ni6,ci6(maxres0)
c!!
	seed=1234567
c
	do is=1,ns
		iseg=seglist(is)
		lclosed(iseg)=(ni(iseg).le.1)
	end do
c
!	write(*,*) 'split',est,ns,top,scorecutoff
	call get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,trans,lseqtl)
	if(est.le.scorecutoff) return
c
	xseg=0
	yseg=0
	emaxim=-infinit
	do is=1,ns
	  iseg=seglist(is)
	  call getemax(ess(iseg,iseg),emaxim,lclosed,iseg,iseg,
     $ xseg,yseg,seed)
	  do js=1,is-1
	    jseg=seglist(js)
	    call getemax(ess(iseg,jseg),emaxim,lclosed,iseg,jseg,xseg,
     $ yseg,seed)
	  end do
	end do
	if((xseg.eq.0).and.(yseg.eq.0)) return
c
	if(start(xseg,yseg).lt.0) then
		maxi=0
		mini=0
		n1=0
	else
		maxi=-infinit
		mini=+infinit
		if(xseg.eq.yseg) then
			iwhere=start(xseg,xseg)
			do xr=1,ni(xseg)
				xres=ci(xr,xseg)
				emax(xres)=ex(iwhere+xres)
			end do
		else if(xseg.lt.yseg) then
			ijstart=start(xseg,yseg)
			jw=mi(yseg)
			do xr=1,ni(xseg)
				e=-infinit
				xres=ci(xr,xseg)
				iwhere=ijstart+(xres-1)*jw
				do yr=1,ni(yseg)
					e=max(e,ex(iwhere+ci(yr,yseg)))
				end do
				emax(xres)=e
			end do
		else if(xseg.gt.yseg) then
			do xr=1,ni(xseg)
				emax(ci(xr,xseg))=-infinit
			end do
			ijstart=start(yseg,xseg)
			jw=mi(xseg)
			do yr=1,ni(yseg)
				iwhere=ijstart+(ci(yr,yseg)-1)*jw
				do xr=1,ni(xseg)
					xres=ci(xr,xseg)
					e=emax(xres)
					emax(xres)=max(e,ex(iwhere+xres))
				end do
			end do
		end if
		n1=ni(xseg)
	end if
	bestxres=ci(1,xseg)
	do xr=1,n1
		xres=ci(xr,xseg)
		e=emax(xres)
		if(xr.gt.2) then
			if(e.gt.emax(bestxres)) bestxres=xres
		end if
		maxi=max(maxi,e)
		mini=min(mini,e)
	end do
c
c	split 90 % - 10 %
c
	midpoint=mini+(maxi-mini)*9/10
	ni1=0
	ni2=0
	if(mini.eq.maxi) then
		do xr=1,ni(xseg),2
			ni1=ni1+1
			ci1(ni1)=ci(xr,xseg)
		end do
		do xr=2,ni(xseg),2
			ni2=ni2+1
			ci2(ni2)=ci(xr,xseg)
		end do
	else
		do xr=1,ni(xseg)
			xres=ci(xr,xseg)
			if(emax(xres).le.midpoint) then
				ni2=ni2+1
				ci2(ni2)=ci(xr,xseg)
			else
				ni1=ni1+1
				ci1(ni1)=ci(xr,xseg)
			end if
		end do
	end if
!	if(ni1.le.0.or.ni2.le.0.or.ni1+ni2.ne.ni(xseg)) then
!	  write(*,*) 'mid:',midpoint,mini,maxi,xseg,yseg,ni(xseg),ni1,ni2
!	  stop
!	end if
c
c	stack better half
c
	if(lseqtl) then 	! split better half left&nul-top&right ci3,c4
		ni3=0
		ni4=0
		ni5=0
		ni6=0
		do xres=1,ni1
			pr=trans(ci1(xres),xseg)
			if(pr.eq.nul) then
				ni6=ni6+1
				ci6(ni6)=ci1(xres)
			else if(ci1(xres).lt.bestxres) then
				ni3=ni3+1
				ci3(ni3)=ci1(xres)
			else
				ni4=ni4+1
				ci4(ni4)=ci1(xres)
			end if
		end do
		if(ni6.gt.0) call copyandput(ni6,ci6,ni,ci,est,ns,
     $			seglist,ess,xseg,ex,start,mi,trans,lseqtl,
     $          	scorecutoff,top,stack,chunk,bestest,link,
     $			link_pointer,forward_pointer,backward_pointer,
     $			firstest)
		if(lskipflag) return
		if(ni5.gt.0) call copyandput(ni5,ci5,ni,ci,est,ns,
     $			seglist,ess,xseg,ex,start,mi,trans,lseqtl,
     $          	scorecutoff,top,stack,chunk,bestest,link,
     $			link_pointer,forward_pointer,backward_pointer,
     $			firstest)
		if(lskipflag) return
		if(ni3.gt.0) call copyandput(ni3,ci3,ni,ci,est,ns,
     $			seglist,ess,xseg,ex,start,mi,trans,lseqtl,
     $          	scorecutoff,top,stack,chunk,bestest,link,
     $			link_pointer,forward_pointer,backward_pointer,
     $			firstest)
		if(lskipflag) return
		if(ni4.gt.0) call copyandput(ni4,ci4,ni,ci,est,ns,
     $			seglist,ess,xseg,ex,start,mi,trans,lseqtl,
     $          	scorecutoff,top,stack,chunk,bestest,link,
     $			link_pointer,forward_pointer,backward_pointer,
     $			firstest)
		if(lskipflag) return
	else
		call copyandput(ni1,ci1,ni,ci,est,ns,
     $			seglist,ess,xseg,ex,start,mi,trans,lseqtl,
     $          	scorecutoff,top,stack,chunk,bestest,link,
     $			link_pointer,forward_pointer,backward_pointer,
     $			firstest)
 		if(lskipflag) return
	end if
c
c       stack lower half
c
	call copyandput(ni2,ci2,ni,ci,est,ns,
     $		seglist,ess,xseg,ex,start,mi,trans,lseqtl,
     $         	scorecutoff,top,stack,chunk,bestest,link,
     $		link_pointer,forward_pointer,backward_pointer,
     $		firstest)
	if(lskipflag) return
!       write(*,*) 'top',top

	return
	end
c
c----------------------------------------------------------------------
c
        subroutine copyandput(ni1,ci1,ni,ci,est,ns,seglist,ess,xseg,ex,
     $          start,mi,trans,lseqtl,scorecutoff,top,stack,chunk,
     $		bestest,link,link_pointer,forward_pointer,
     $          backward_pointer,firstest)
        implicit none
        include 'parsizes.for'
        integer re_estimate,start(maxseg,maxseg),trans(maxres0,maxseg)
        integer ns,seglist(maxseg),mi(maxseg),chunk,bestest
	integer link(2,maxstack2)
        integer est,ex(exdim),scorecutoff,top,stack(maxstack)
        integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $          backward_pointer(0:boxdim)
        real firstest
        logical lseqtl
        integer xres,ess(maxseg,maxseg),f,ni(maxseg),ci(maxres0,maxseg)
        integer xseg,ni1,ci1(maxres0)
c
        if(ni1.eq.0) return
        ni(xseg)=ni1
        do xres=1,ni1
                ci(xres,xseg)=ci1(xres)
        end do
!        call excludeimpossible(ni,ci,ns,seglist,trans,lseqtl)
        f=re_estimate(est,ns,seglist,ess,xseg,ex,ni,ci,start,
     $          mi,trans,lseqtl)
        if(f.gt.scorecutoff) then
                top=top+chunk
                call putrequest(f,ni,ci,ns,seglist,mi,top,stack,chunk,
     $                  bestest,link,link_pointer,
     $                  forward_pointer,backward_pointer,firstest)
		if(lskipflag) return
        end if

        return
        end
c
c----------------------------------------------------------------------
c
	subroutine excludeimpossible(ni,ci,ns,seglist,trans,lseqtl)
	implicit none
	include 'parsizes.for'
	integer ni(maxseg),ci(maxres0,maxseg),trans(maxres0,maxseg)
	integer ns,seglist(maxseg)
	logical lseqtl
c
	integer is,iseg,leftmost,ritemost,ir,irx,i,nx
	logical lhasnul(maxseg)
c
	if(.not.lseqtl) return
	nx=0
	do is=1,ns
		iseg=seglist(is)
		lhasnul(iseg)=(trans(ci(ni(iseg),iseg),iseg).eq.nul)
!		type *,'before',iseg,(trans(ci(ir,iseg),iseg),ir=1,ni(iseg))
	end do
	leftmost=-infinit
	do is=1,ns
	  iseg=seglist(is)
	  if(.not.lhasnul(iseg)) then
		irx=0
		do ir=1,ni(iseg)
			i=trans(ci(ir,iseg),iseg)
!			if(i.le.leftmost.and.i.ne.nul) irx=ir
			if(i.lt.leftmost.and.i.ne.nul) irx=ir
		end do
		nx=nx+irx
		! exclude 1...irx
		if(irx.gt.0) then
			i=0
			do ir=irx+1,ni(iseg)
				i=i+1
				ci(i,iseg)=ci(ir,iseg)
			end do
			ni(iseg)=i
		end if
		if(ni(iseg).gt.0) then
			if(.not.lhasnul(iseg)) then
				leftmost=trans(ci(1,iseg),iseg)
			else if (ni(iseg).gt.1) then
				leftmost=trans(ci(1,iseg),iseg)
			end if
		else
			leftmost=+infinit
			goto 99
		end if
		i=ci(1,iseg)
!		type *,'leftmost:',lhasnul(iseg),iseg,leftmost,i,
!     $			trans(i,iseg),irx
	  end if
	end do
	if(leftmost.eq.infinit) goto 99
10	ritemost=+infinit
	do is=ns,1,-1
		iseg=seglist(is)
		if(.not.lhasnul(iseg)) then
			irx=0
			do ir=ni(iseg),1,-1
				i=trans(ci(ir,iseg),iseg)
!				if(i.ge.ritemost.and.i.ne.nul) irx=ir
				if(i.gt.ritemost.and.i.ne.nul) irx=ir
			end do
			if(irx.gt.0) nx=nx+ni(iseg)-irx
			! exclude irx...ni(iseg) except nul
			if(irx.gt.0) then
				if(lhasnul(iseg)) then
					ci(irx,iseg)=ci(ni(iseg),iseg)
					ni(iseg)=irx
			  	else
					ni(iseg)=irx-1
				end if
			end if
			if(ni(iseg).gt.0) then
				if(.not.lhasnul(iseg)) then
				  ritemost=trans(ci(ni(iseg),iseg),iseg)
				else if(ni(iseg).gt.1) then
				  ritemost=trans(ci(ni(iseg)-1,iseg),iseg)
				end if
			else
				ritemost=-infinit
				goto 99
			end if
		end if
		i=ci(ni(iseg),iseg)
!		type *,'ritemost:',iseg,ritemost,i,trans(i,iseg),irx
	end do

99	continue
!	if(nx.gt.0) then
!	  do is=1,ns
!		iseg=seglist(is)
!		type *,'after',iseg,(trans(ci(ir,iseg),iseg),ir=1,ni(iseg))
!	  end do
!	  type *,'nx',nx
!	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getnextbest(scorecutoff,stack,top,chunk,ns,seglist,
     $		mi,bestest,ali,flag,link,ex,start,est,link_pointer,
     $		forward_pointer,backward_pointer,firstest,trans,lseqtl)
c
c	flags:	1 = empty stack
c		2 = top score below cutoff
c		4 = no candidates
c		5 = unique alignment, returned in ali()
c		3 = overflow
c
	implicit none
	include 'parsizes.for'
	integer scorecutoff,flag,stack(maxstack),top,chunk,ns,mi(maxseg)
	integer link(2,maxstack2),ali(maxseg),bestest,
     $ seglist(maxseg),ex(exdim)
	integer start(maxseg,maxseg),trans(maxres0,maxseg)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
	logical lseqtl
c
	integer is,est,nmin,nmax,ni(maxseg),ci(maxres0,maxseg),ir,iseg,i,m
c	write(*,*) 'getnextbest',scorecutoff,top,chunk,ns,bestest,firstest
	est=-infinit
	flag=0
c	check that bestest is!
c	call checkstack(stack,bestest,top,chunk,'getnext BUG ')
	do while(flag.eq.0)
c
c		check exceptions
c
		if(top.le.0) then
			flag=1
			est=-infinit
			return
		end if
		if(bestest.le.0) then
	write(*,*) scorecutoff,(stack(i),i=1,10),top,chunk,bestest
	write(*,*) (link(1,i),link(2,i),i=1,top)
	write(*,*) 'FATAL ERROR -- bestest out of range'
			lskipflag=.true.
			return
		end if
		if(stack(bestest).lt.scorecutoff) flag=2
		if(top+2*chunk.gt.maxstack) flag=3
c!!
		if(top.gt.maxstack-100) flag=3
		if(flag.ne.0) return
c
c		continue splitting / export unique alignment
c
		call toprequest(bestest,top,stack,chunk,link,mi,ns,seglist,
     $			ni,ci,est,link_pointer,forward_pointer,
     $			backward_pointer,firstest)
		if(lskipflag) return
c		write(*,*) 'top search space'
c		do is=1,ns
c			iseg=seglist(is)
c			write(*,*) is,iseg,mi(iseg),ni(iseg),
c     $				(ci(ir,iseg),ir=1,ni(iseg))
c		end do
		call getnminnmax(ni,ns,seglist,nmin,nmax)
		if(nmin.eq.0) flag=4
		if((nmin.eq.1).and.(nmax.eq.1)) flag=5
		if(flag.eq.0) call split(est,ns,seglist,ex,top,stack,
     $			scorecutoff,ni,ci,mi,chunk,bestest,link,start,
     $			link_pointer,forward_pointer,backward_pointer,firstest,
     $			trans,lseqtl)
		if(lskipflag) return
		if(flag.eq.5) then
			do is=1,ns
				ali(seglist(is))=ci(1,seglist(is))
			end do
		end if
c		write(*,*) 'flag is',flag
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine lpack(lwork,i,l)
	implicit none
	include 'parsizes.for'
	logical lwork(maxres0*maxseg),l(31)
	integer i,j

	do j=0,30
		l(j+1)=lwork(i+j)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine putrequest(est,ni,ci,ns,seglist,mi,to,stack,
     $ chunk,bestest,link,link_pointer,forward_pointer,
     $ backward_pointer,firstest)
	implicit none
	include 'parsizes.for'
	integer ni(maxseg),ci(maxres0,maxseg),ns,mi(maxseg)
	integer seglist(maxseg)
	integer to,stack(maxstack),chunk,est,bestest,link(2,maxstack2)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer i,j,is,ir,irray((1+(maxres0+30)/31)*maxseg),s,iseg,n
	logical lwork(maxres0*maxseg),l(31),lw(maxres0)
c
c	write(*,*) 'putrequest',est,ns,to,chunk,bestest,firstest
c	call checkstack(stack,bestest,to,chunk,'put BUG 1')
	n=0
	do is=1,ns
		iseg=seglist(is)
		do ir=1,mi(iseg)
			lw(ir)=.false.
		end do
		do ir=1,ni(iseg)
			lw(ci(ir,iseg))=.true.
		end do
		do ir=1,mi(iseg)
			n=n+1
			lwork(n)=lw(ir)
		end do
	end do
	do i=n+1,min(n+30,maxres0*maxseg)
		lwork(i)=.false.
	end do
	j=1
	irray(1)=est
	do i=1,n,31
		call lpack(lwork,i,l)
		call bitpack(s,l,31)
		j=j+1
		irray(j)=s
c		write(*,*) 'l',i,j,s
	end do
c	write(*,*) 'put irray',(irray(i),i=1,j)
	if(to+chunk.gt.maxstack) then
	write(*,*)  'FATAL ERROR in putrequest'
		lskipflag=.true.
		return
	end if
	call putchunk(to,irray,stack,chunk)
	if(lskipflag) return
	call createlink(to,bestest,link,stack,link_pointer,
     $		forward_pointer,backward_pointer,firstest,to,chunk)
	if(lskipflag) return
c	call checkstack(stack,bestest,to,chunk,'put BUG 2')

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine declump(scorecutoff,stack,top,chunk,link,bestest,ns,
     $		seglist,mi,ali,ex,start,singletcutoff,link_pointer,
     $		forward_pointer,backward_pointer,firstest,trans,lseqtl)
	implicit none
	include 'parsizes.for'
	integer stack(maxstack),top,chunk,link(2,maxstack2),bestest
	integer scorecutoff,ns,mi(maxseg),ali(maxseg),seglist(maxseg)
	integer ex(exdim),start(maxseg,maxseg),singletcutoff
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim),trans(maxres0,maxseg)
	real firstest
	logical lseqtl
c
	integer i,irray((1+(maxres0+30)/31)*maxseg),est,ci(maxres0,maxseg)
	integer ni(maxseg),oldtop,is,ir,jr,ess(maxseg,maxseg),
     $ ires,iseg,nold
	logical lchange
	integer m
c
c
c	write(*,500) scorecutoff,top,chunk,bestest,ns,(ali(seglist(i)),i=1,ns)
	oldtop=top
	top=1-chunk
	bestest=0
	call init_link_pointer(link_pointer,forward_pointer,
     $ backward_pointer)
	do i=1,oldtop,chunk
		call getchunk(i,irray,stack,chunk)
		call translatechunk(irray,mi,ns,seglist,ni,ci,est,chunk)
c
c		only declump chunks below singletcutoff !
c
c????
		if(est.gt.singletcutoff) goto 20
c
c		blank out ali()
c
		lchange=.false.
		do is=1,ns
			iseg=seglist(is)
			nold=ni(iseg)
			m=ali(iseg)
			do ir=1,nold
			  ires=ci(ir,iseg)
			  if(ires.eq.m) then
			    if(trans(ires,iseg).ne.nul) then
c				write(*,*) 'remove',ires,iseg
				do jr=ir+1,ni(iseg)
					ci(jr-1,iseg)=ci(jr,iseg)
				end do
				ni(iseg)=ni(iseg)-1
c				no candidates left !
				if(ni(iseg).eq.0) then
					est=-infinit
					goto 29
				end if
				lchange=.true.
				goto 19
			    end if
			  end if
			end do
19			continue
		end do
c
c		re_estimate score if changes made
c
		if(lchange) call get_ess(ess,est,ns,seglist,ni,ci,ex,start,mi,
     $			trans,lseqtl)
20		if(est.gt.scorecutoff) then
			top=top+chunk
			call putrequest(est,ni,ci,ns,seglist,mi,
     $ top,stack,chunk,bestest,link,link_pointer,
     $ forward_pointer,backward_pointer,firstest)
			if(lskipflag) return
		else
	end if
29	end do
500	format('declump',4i12,85i4)
c   500 format('declump',4i12,<maxseg+5>i4)
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine clearstack(cutoff,stack,top,chunk,link,bestest,
     $ link_pointer,forward_pointer,backward_pointer,firstest)
c
c	remove all chunks that have a score below cutoff
c
	implicit none
	include 'parsizes.for'
	integer cutoff,stack(maxstack),top,chunk,link(2,maxstack2),bestest
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer i,dn,oldtop,x,getk,e
c
c	compress chunks
c
c	write(*,*) 'clearstack in ',cutoff,top,chunk,bestest
	oldtop=top
	top=1-chunk
	do i=1,oldtop,chunk
		if((stack(i).gt.cutoff).and.(top.lt.i)) then
			top=top+chunk
			call swapstack(i,top,stack,chunk,link,bestest,
     $			 link_pointer,forward_pointer,backward_pointer,firstest)
			if(lskipflag) return
		end if
	end do
	call init_link_pointer(link_pointer,forward_pointer,
     $ backward_pointer)
	do i=1,top,chunk
		e=stack(i)
		x=getk(e,i,stack,link_pointer,forward_pointer,
     $			backward_pointer,firstest)
		if(lskipflag) return
	end do
c
c	curtail linked list
c
	do i=1,top,chunk
		dn=link(1,i)
		if((dn.gt.0).and.(stack(dn).le.cutoff)) then
			link(1,i)=0
c			write(*,*) 'clearstack out',cutoff,top,chunk,bestest
			return
		end if
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine init_link_pointer(link_pointer,forward_pointer,
     $		backward_pointer)
	implicit none
	include 'parsizes.for'
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
c
	integer i
c
	do i=0,boxdim
		link_pointer(i)=0
		forward_pointer(i)=boxdim
		backward_pointer(i)=0
	end do
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine createlink_ok(i,bestest,link,stack,link_pointer,
     $		forward_pointer,backward_pointer,firstest,top,chunk)
c
c	link in new chunk 'i'
c
	implicit none
	include 'parsizes.for'
	integer i,bestest,link(2,maxstack2),stack(maxstack),top,chunk
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer j,k,est
c
	if(i.gt.maxstack2) then
	write(*,*)  'FATAL ERROR link overflow - increase maxstack2 '
		lskipflag=.true.
		return
	end if
	link(1,i)=0
	link(2,i)=0
	est=stack(i)
	j=bestest
	k=0
	do while(j.gt.0)
		if(stack(j).gt.est) k=j
		j=link(1,j)
	end do
	if(k.eq.0) then
		if(bestest.ne.0) then
			link(1,i)=bestest
			link(2,bestest)=i
		end if
		bestest=i
	else
		j=link(1,k)
		link(1,i)=j
		if(j.ne.0) link(2,j)=i
		link(1,k)=i
		link(2,i)=k
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine createlink(i,bestest,link,stack,link_pointer,
     $		forward_pointer,backward_pointer,firstest,top,chunk)
c
c	link in new chunk 'i'
c
	implicit none
	include 'parsizes.for'
	integer i,bestest,link(2,maxstack2),stack(maxstack)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim),top,chunk
	real firstest
c
	integer j,k,est,getk,x
c
c	write(*,*) 'createlink',i,bestest,firstest
	if(i.gt.maxstack2) then
	write(*,*)  'FATAL ERROR link overflow - increase maxstack2 '
		lskipflag=.true.
		return
	end if
	link(1,i)=0
	link(2,i)=0
	est=stack(i)
!	j=bestest
!	k=0
!	do while(j.gt.0)
!		if(stack(j).gt.est) k=j
!		j=link(1,j)
!	end do
!	x=k
	j=getk(est,i,stack,link_pointer,forward_pointer,
     $		backward_pointer,firstest)
	if(lskipflag) return
c	write(*,*) 'create',bestest,stack(bestest),j,stack(j),i,est
	k=0
	if(j.eq.0) j=bestest
	do while(j.gt.0)
		if(stack(j).gt.est) k=j
		j=link(1,j)
	end do
!	if(k.ne.x) write(*,*) 'BUG',x,stack(max(1,x)),bestest,
!     $		stack(max(1,bestest)),i,est,k,stack(max(1,k)),firstest
!	k=x
	if(k.eq.0) then
		if(bestest.ne.0) then
			link(1,i)=bestest
			link(2,bestest)=i
		end if
		bestest=i
c		call checkstack(stack,bestest,top,chunk,'create BUG 1')
	else
		j=link(1,k)
		link(1,i)=j
		if(j.ne.0) link(2,j)=i
		link(1,k)=i
		link(2,i)=k
c		call checkstack(stack,bestest,top,chunk,'create BUG 2')
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	function getk(est,i,stack,link_pointer,forward_pointer,
     $		backward_pointer,firstest)
	implicit none
	include 'parsizes.for'
	integer getk,est,i,stack(maxstack)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer box,k,b,funbox
c
c	returns stack-pointer to a next larger estimate, 0 if new best
c
c	box=funbox(est,firstest)
	box=max(1,int(float(est)/firstest*float(boxdim-1)))
	if(box.gt.boxdim) then
	write(*,*)  'FATAL ERROR -- boxfun overflow'
		lskipflag=.true.
		return
	end if
	if(link_pointer(box).eq.0) then
		do b=box+1,forward_pointer(box)
			backward_pointer(b)=box
		end do
		do b=backward_pointer(box),box-1
			forward_pointer(b)=box
		end do
		link_pointer(box)=i
		box=forward_pointer(box)
	else
		if(est.ge.stack(link_pointer(box))) then
			link_pointer(box)=i
			box=forward_pointer(box)
		end if
	end if
	if(box.eq.boxdim) then
		k=0
	else
		k=link_pointer(box)
	end if

	getk=k
c	write(*,*) 'getk',est,i,k,box,(link_pointer(b),b=0,boxdim)

	return
	end
c
c----------------------------------------------------------------------
c
	function funbox(est,firstest)
	implicit none
	include 'parsizes.for'
	integer est,funbox
	real firstest,r
c
	r=float(est)/firstest
	funbox=max(1,int(r*r*(boxdim-1)))
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine removebox(est,i,stack,link,link_pointer,
     $ forward_pointer,backward_pointer,firstest)
	implicit none
	include 'parsizes.for'
	integer est,i,link(2,maxstack2),stack(maxstack)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer b,box,k,bwp,fwp,funbox
c
c	box=funbox(est,firstest)
	box=max(1,int(float(est)/firstest*(boxdim-1)))
c	case: middle occupant of box -> no action
	if(link_pointer(box).ne.i) return
	k=link(1,i)
c	case: top is last in box, last in stack
	if(k.eq.0) then
		link_pointer(box)=0
		bwp=backward_pointer(box)
		fwp=forward_pointer(box)
		do b=box+1,fwp
			backward_pointer(b)=bwp
		end do
		do b=bwp,box-1
			forward_pointer(b)=fwp
		end do
		return
	end if
	if(k.gt.0) b=max(1,int(float(stack(k))/firstest*(boxdim-1)))
c	case: last in box
	if(b.ne.box) then
		link_pointer(box)=0
		bwp=backward_pointer(box)
		fwp=forward_pointer(box)
		do b=box+1,fwp
			backward_pointer(b)=bwp
		end do
		do b=bwp,box-1
			forward_pointer(b)=fwp
		end do
		return
	end if
c	case: remainders in box
	if(b.eq.box) then
		link_pointer(box)=k
		return
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine toprequest(bestest,top,stack,chunk,link,mi,ns,seglist,
     $		ni,ci,est,link_pointer,forward_pointer,backward_pointer,
     $		firstest)
	implicit none
	include 'parsizes.for'
	integer bestest,top,stack(maxstack),chunk,link(2,maxstack2),ns
	integer ni(maxseg),ci(maxres0,maxseg),est,mi(maxseg),
     $ seglist(maxseg)
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer irray((1+(maxres0+30)/31)*maxseg),s,a
c
c	move best chunk to top of stack
c	translate chunk contents to ni,ci
c	decrement stack counter
c
c	call checkstack(stack,bestest,top,chunk,'top BUG 1')
	a=bestest
	if(top.ne.bestest) call swapstack(a,top,stack,chunk,link,bestest,
     $		link_pointer,forward_pointer,backward_pointer,firstest)
	if(lskipflag) return
c	check that top is really maximum!
c	call checkstack(stack,bestest,top,chunk,'top BUG 2')
	call getchunk(top,irray,stack,chunk)
	call translatechunk(irray,mi,ns,seglist,ni,ci,est,chunk)
	s=stack(bestest)
	call removebox(s,bestest,stack,link,link_pointer,
     $		forward_pointer,backward_pointer,firstest)
c	call checkstack(stack,bestest,top,chunk,'top BUG 3')
	call removetoplink(top,chunk,link,bestest)
c	call checkstack(stack,bestest,top,chunk,'top BUG 4')
c	write(*,*) 'toprequest',bestest,top,chunk,ns,est
c	write(*,*) 'ni',(ni(seglist(is)),is=1,ns)

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine checkstack(stack,bestest,top,chunk,text)
	implicit none
	include 'parsizes.for'
	integer bestest,top,stack(maxstack),chunk
	character*40 text
c
	integer i,m
c
!	type *,'checkstack',bestest,top,chunk,text
	if(bestest.gt.0) then
	  m=stack(bestest)
	  do i=1,top-1,chunk
		if(stack(i).gt.m) then
			write(*,*) top,bestest,stack(top),stack(bestest)
			write(*,*) (stack(m),m=1,top,chunk)
			write(*,*) text
	write(*,*) ' FATAL ERROR in checkstack'
			lskipflag=.true.
			return
		end if
	  end do
	end if
c
	return
	end
c
c----------------------------------------------------------------------
c
	subroutine removetoplink(top,chunk,link,bestest)
	implicit none
	include 'parsizes.for'
	integer link(2,maxstack2),bestest,top,chunk
c
	bestest=link(1,bestest)
	if(bestest.gt.0) link(2,bestest)=0
	top=top-chunk

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine translatechunk(irray,mi,ns,seglist,ni,ci,est,chunk)
	implicit none
	include 'parsizes.for'
	integer irray((1+(maxres0+30)/31)*maxseg),mi(maxseg)
	integer ns,est,ni(maxseg)
	integer chunk,seglist(maxseg),ci(maxres0,maxseg)
c
	integer i,j,ires,is,iseg,k
	logical l(31)
c
	est=irray(1)
	do is=1,ns
		ni(seglist(is))=0
	end do
	is=0
	ires=0
	do i=2,chunk
		call bitunpack(irray(i),l,1073741823,31)
		do j=1,31
			if(ires.eq.0) then
				is=is+1
				if(is.gt.ns) return
				iseg=seglist(is)
			end if
			ires=ires+1
			if(l(j)) then
				ni(iseg)=ni(iseg)+1
				ci(ni(iseg),iseg)=ires
			end if
c			zero ires counter before next segment
			if(ires.eq.mi(iseg)) ires=0
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine swaplinks_new(a,b,link,bestest)
	implicit none
	include 'parsizes.for'
	integer a,b,link(2,maxstack2),bestest
c
	integer dn1,dn2,up1,up2
c
	dn1=link(1,a)
	up1=link(2,a)
	dn2=link(1,b)
	up2=link(2,b)

	if(up1.eq.b.or.up2.eq.a) then	! adjacent --a-b--
		link(1,a)=up2
		link(1,b)=a
		if(dn1.ne.0) link(1,dn1)=b
		link(2,a)=b
		link(2,b)=dn1
		if(up2.ne.0) link(2,up2)=a
	else			! non-adjacent --a-x-b--
		link(1,b)=dn1
		link(1,a)=dn2
		link(2,b)=up1
		link(2,a)=up2
		if(dn1.ne.0) link(2,dn1)=a
		if(dn2.ne.0) link(2,dn2)=b
		if(up1.ne.0) link(1,up1)=b
		if(up2.ne.0) link(1,up2)=a
	end if

	if(up1.eq.0) bestest=b
	if(up2.eq.0) bestest=a

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine swaplinks(a,b,link,bestest)
	implicit none
	include 'parsizes.for'
	integer bestest,a,b,link(2,maxstack2)
c
	integer dn1,dn2,up1,up2
c
	dn1=link(1,a)
	up1=link(2,a)
	if(dn1.eq.b) dn1=a
	dn2=link(1,b)
	up2=link(2,b)
	if(up2.eq.a) up2=b
	link(1,b)=dn1
	link(2,b)=up1
	link(1,a)=dn2
	link(2,a)=up2
	if(dn1.ne.0) link(2,dn1)=b
	if(up1.ne.0) link(1,up1)=b
	if(dn2.ne.0) link(2,dn2)=a
	if(up2.ne.0) link(1,up2)=a
	if(up1.eq.0) bestest=b
	if(up2.eq.0) bestest=a

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine swapstack(from,to,stack,chunk,link,bestest,
     $ link_pointer,forward_pointer,backward_pointer,firstest)
	implicit none
	include 'parsizes.for'
	integer from,to,stack(maxstack),chunk,link(2,maxstack2),bestest
	integer link_pointer(0:boxdim),forward_pointer(0:boxdim),
     $		backward_pointer(0:boxdim)
	real firstest
c
	integer i,a,b,e,getk
c
c	write(*,*) 'swap',from,stack(from),to,stack(to),firstest
	if(from.eq.to) return
	e=stack(from)
	call removebox(e,from,stack,link,link_pointer,
     $		forward_pointer,backward_pointer,firstest)
	e=stack(to)
	call removebox(e,to,stack,link,link_pointer,
     $		forward_pointer,backward_pointer,firstest)
	a=from
	b=to
	do i=1,chunk
		e=stack(b)
		stack(b)=stack(a)
		stack(a)=e
		a=a+1
		b=b+1
	end do
	call swaplinks(from,to,link,bestest)
	e=stack(from)
	i=getk(e,from,stack,link_pointer,
     $		forward_pointer,backward_pointer,firstest)
	if(lskipflag) return
	e=stack(to)
	i=getk(e,to,stack,link_pointer,
     $		forward_pointer,backward_pointer,firstest)
	if(lskipflag) return

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine putchunk(to,irray,stack,chunk)
	implicit none
	include 'parsizes.for'
	integer to,stack(maxstack),chunk,irray(chunk)
c
	integer i,a
c
	a=to
	if(a+chunk.gt.maxstack) then
		write(*,*) to,chunk
	write(*,*)  'FATAL ERROR stack overflow - increase maxstack'
		lskipflag=.true.
		return
	end if
	do i=1,chunk
		stack(a)=irray(i)
		a=a+1
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getchunk(from,irray,stack,chunk)
	implicit none
	include 'parsizes.for'
	integer from,stack(maxstack),chunk,irray(chunk)
c
	integer i,a
c
	a=from
	do i=1,chunk
		irray(i)=stack(a)
		a=a+1
	end do

	return
	end
c
c----------------------------------------------------------------------
c
	subroutine getnminnmax(ni,ns,seglist,nmin,nmax)
	implicit none
	include 'parsizes.for'
	integer ns,ni(maxseg),nmin,nmax,seglist(maxseg)
c
	integer is,n
c
	nmin=infinit
	nmax=-infinit
	do is=1,ns
		n=ni(seglist(is))
		nmin=min(nmin,n)
		nmax=max(nmax,n)
	end do

	return
	end
c
c----------------------------------------------------------------------
c
