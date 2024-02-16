	program puutos
	implicit none
	include 'parsizes-puutos.for'
        integer ndom,node_child(2,maxdom),node_range(maxdom,2,maxseg)
        integer node_parent(maxdom),node_size(maxdom),node_nseg(maxdom)
	integer node_size_s(maxdom),node_nseg_s(maxdom),
     $		node_range_s(maxdom,2,maxseg),acc(maxres)
        character chainid,node_type(maxdom),secstr(maxseg)
        character*4 code
        integer ierr,nres,nseg,bot(maxseg)
	integer domns(maxdom),domseglist(maxseg,maxdom),na,nb,countchr
	character struc(maxres)
	character*80 findfile,dsspfile,dalidatpath,puufile
	real ca(3,maxres)
	integer segmentrange(2,maxseg),checkrange(2,maxseg)
        integer checkx(maxseg)
	character seq(maxres)
	character*80 compnd
        character*5 resno(0:maxres)
c
c	assign segment units using linearpuu tree & segment midpoints
c	1. read units.puu
c	2. read segment data
c	3. assign bottom unit to residues
c	4. find segment midpoints
c	5. assign segments to bottom units
c	6. prune tree of empty / redundant units
c
c	do SINGLE protein
c
        if(iargc().lt.3) stop "USAGE: puutos units.puu dsspfile dat_1"
        call getarg(1,puufile)
        call getarg(2,dsspfile)
        call getarg(3,dalidatpath)
c
	open(90,file=puufile,status='old')
	  call readnexttree(90,ndom,node_child,node_range,
     $      node_nseg,node_parent,node_size,code,chainid,node_type,ierr)
	  write(*,*) 'readnexttree returned ',ierr
	  if(ierr.ne.0) goto 19
	  call readdssp1(dsspfile,chainid,91,nres,struc,ca,seq,compnd,
     $          acc,resno)
	  call getsecstr(nres,struc,nseg,secstr,segmentrange,6,8)
	  na=countchr('H',secstr,nseg)
	  nb=countchr('E',secstr,nseg)
	  call getcheckrange(nres,nseg,segmentrange,secstr,checkrange,
     $		checkx)
	  call assignsegments(nseg,segmentrange,ndom,node_nseg,
     $          node_range,bot)
	  if(ndom.gt.0) call prunetree(nseg,bot,ndom,node_child,
     $          node_parent,
     $		node_size_s,node_nseg_s,node_range_s,node_type,
     $		node_size,node_nseg,node_range)
	  call getdomseglist(domns,domseglist,ndom,node_nseg_s,
     $          node_range_s)
c
c	  write database for parsi:
c		>>>> code chainid nres nseg na nb secstr
c		iseg segmentrange checkrange checkx
c		ca
c		>>>> code chainid ndom
c		idom node_type node_child domns domseglist
c
	  call writeit(code,chainid,nres,nseg,na,nb,secstr,segmentrange,
     $	checkrange,checkx,ca,ndom,node_type,node_child,domns,
     $	domseglist,92,node_size,node_nseg,node_range,dalidatpath,
     $	seq,compnd,struc,acc,resno)
19	close(90)

500	format(a80)

	end
c
c-------------------------------------------------------------------------------
c
	subroutine writeit(code,chainid,nres,nseg,na,nb,secstr,
     $  segmentrange,
     $	checkrange,checkx,ca,ndom,node_type,node_child,domns,
     $	domseglist,iunit,node_size_o,node_nseg_o,node_range_o,
     $	dalidatpath,seq,compnd,struc,acc,resno)
	implicit none
	include 'parsizes-puutos.for'
        integer ndom,node_child(2,maxdom),nres,nseg,iunit,acc(maxres)
        character chainid,node_type(maxdom),secstr(maxseg),seq(maxres),
     $          struc(maxres)
	character*80 compnd
        character*4 code
	integer domns(maxdom),domseglist(maxseg,maxdom),na,nb
	real ca(3,maxres)
	integer segmentrange(2,maxseg),checkrange(2,maxseg)
	integer node_size_o(maxdom),node_nseg_o(maxdom),checkx(maxseg),
     $		node_range_o(maxdom,2,maxseg)
c
	character*80 filnam,constructfilnam,dalidatpath
	integer i,j
	character*5 cd,resno(0:maxres)
c
	cd(1:4)=code
	cd(5:5)=chainid
	filnam=constructfilnam(cd,dalidatpath,'.dat')
c	open(iunit,file=filnam,status='new',err=99)
        open(iunit,file=filnam,err=99)
        write(iunit,500) code,chainid,nres,nseg,na,nb,
     $          (secstr(i),i=1,nseg)
	do i=1,nseg
	  write(iunit,510) i,(segmentrange(j,i),j=1,2),
     $          (checkrange(j,i),j=1,2),checkx(i)
	end do
	write(iunit,520) ((ca(j,i),j=1,3),i=1,nres)
	write(iunit,500) code,chainid,ndom
	do i=1,ndom
	  write(iunit,530) i,node_type(i),
     $          (node_child(j,i),j=1,2),domns(i),
     $	        (domseglist(j,i),j=1,domns(i))
	end do
	write(iunit,500) code,chainid,ndom
	do i=1,ndom
	  write(iunit,530) i,node_type(i),(node_child(j,i),j=1,2),
     $	node_size_o(i),node_nseg_o(i),
     $	(node_range_o(i,1,j),node_range_o(i,2,j),j=1,node_nseg_o(i))
	end do
        ! blank > underscore in dssp-string
        do i=1,nres
                if(struc(i) .eq. ' ') struc(i)='L'
        end do
        write(iunit,532) (struc(i),i=1,nres)
	write(iunit,540) (seq(i),i=1,nres)
	write(iunit,550) compnd(11:80)
        write(iunit,560) (i,resno(i),acc(i),(ca(j,i),j=1,3),i=1,nres)
	close(iunit)

500	format('>>>> ',a4,a1,4i5,2x,10000a1)
510	format(6i10)
520	format(10f8.1)
530	format(i4,1x,a1,1x,3i4,10000i4)
532     format('-dssp     "',64000a1,'"')
540	format('-sequence "',64000a1,'"')
550	format('-compnd   "',a70,'"')
560     format('-acc ',i4,2x,a5,i6,3f12.3)
c	normal exit
	return

c	error exit
99	write(*,*) 'error opening file -- skipped',filnam

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine getcheckrange(nres,nseg,segmentrange,secstr,
     $		checkrange,checkx)
	implicit none
	include 'parsizes-puutos.for'
	integer nseg,segmentrange(2,maxseg),checkrange(2,maxseg)
	integer nres,checkx(maxseg)
	character secstr(maxseg)
c
	integer iseg,l,x
c
c	define checkranges: helix -> 9-10 residues; strand -> 5-6 residues
c
	do iseg=1,nseg
	  l=segmentrange(2,iseg)-segmentrange(1,iseg)+1
	  x=0
	  if(secstr(iseg).eq.'H') then
		if(l.ge.9) x=0
		if ((l.eq.8).or.(l.eq.7)) x=1
		if(l.eq.6) x=2
	  else
		if(l.ge.5) x=0
		if((l.eq.3).or.(l.eq.4)) x=1
		if(l.eq.2) x=2
	  end if
	  checkx(iseg)=x
	  checkrange(1,iseg)=max(1,segmentrange(1,iseg)-x)
	  checkrange(2,iseg)=min(nres,segmentrange(2,iseg)+x)
	end do

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine prunetree(nseg,bot,ndom,node_child,node_parent,
     $		node_size,node_nseg,node_range,node_type,
     $		node_size_o,node_nseg_o,node_range_o)
        implicit none
        include 'parsizes-puutos.for'
	integer nseg,bot(maxseg),ndom,node_child(2,maxdom)
	integer node_parent(maxdom)
	character node_type(maxdom)
	integer node_size_o(maxdom),node_nseg_o(maxdom),
     $		node_range_o(maxdom,2,maxseg)
c
	integer iseg,idom,node_nseg(maxdom),node_range(maxdom,2,maxseg)
	logical lactive(maxres)
	integer pnode,inode,jnode,i,j,k,n,list(0:maxres*2)
	integer tsil(0:maxres*2)
	integer parent(maxdom),child(maxdom,2),node_size(maxdom)
c
c	write(*,*) 'prunetree:',(bot(iseg),iseg=1,nseg)
c
c	copy _o to work segment arrays
c
	do idom=1,ndom
		node_size(idom)=node_size_o(idom)
		node_nseg(idom)=node_nseg_o(idom)
		do j=1,node_nseg_o(idom)
			node_range(idom,1,j)=node_range_o(idom,1,j)
			node_range(idom,2,j)=node_range_o(idom,2,j)
		end do
	end do
c
c	check consistency
c
	if(node_range(1,2,1).lt.nseg) then
	write(*,*) ' inconsistent puu/dssp data',node_range(1,2,1),nseg
	write(*,*) 'resetting nseg to ',node_range(1,2,1)
		nseg=node_range(1,2,1)
	end if
c
c	print tree
c
c	write(*,*) 'tree before pruning'
c	do idom=1,ndom
c	 	write(*,*) idom,node_parent(idom),
c     $	  		(node_child(i,idom),i=1,2),node_nseg(idom),
c     $			((node_range(idom,i,j),i=1,2),j=1,node_nseg(idom))
c	end do
c
c	given bot, reconstruct node_nseg,node_range
c
c	initialize bottom
c
c	write(*,*) 'initialize bottom ',ndom,nseg,(bot(iseg),iseg=1,nseg)
	do idom=1,ndom
		node_nseg(idom)=0
		lactive(idom)=.true.
	end do
	do iseg=1,nseg
		idom=bot(iseg)
		if(idom.gt.0) then
		  node_nseg(idom)=node_nseg(idom)+1
		  node_range(idom,1,node_nseg(idom))=iseg
		  node_range(idom,2,node_nseg(idom))=iseg
		end if
	end do
c
c	build parent ranges
c
	write(*,*) 'build parent ranges'
	do idom=ndom,1,-1
	    pnode=node_parent(idom)
	    if(pnode.gt.0) then
		inode=node_child(1,pnode)
		jnode=node_child(2,pnode)
	    else
		inode=1
		jnode=0
	    end if
	    if(inode.eq.idom.and.pnode.gt.0) then
		do iseg=1,node_nseg(inode)
	  node_nseg(pnode)=node_nseg(pnode)+1
	  node_range(pnode,1,node_nseg(pnode))=node_range(inode,1,iseg)
	  node_range(pnode,2,node_nseg(pnode))=node_range(inode,2,iseg)
		end do
		if(jnode.gt.0) then
		 do iseg=1,node_nseg(jnode)
	  node_nseg(pnode)=node_nseg(pnode)+1
	  node_range(pnode,1,node_nseg(pnode))=node_range(jnode,1,iseg)
	  node_range(pnode,2,node_nseg(pnode))=node_range(jnode,2,iseg)
		 end do
		end if
	    end if
c     write(*,499) idom,' ',pnode,inode,jnode,node_nseg(idom),
c  $			(node_range(idom,1,iseg),iseg=1,node_nseg(idom))
	end do
c
c	delete empty nodes
c
	write(*,*) 'delete empty nodes'
	do idom=ndom,1,-1
	  pnode=node_parent(idom)
	  if(lactive(idom).and.pnode.gt.0) then
	    inode=node_child(1,pnode)
	    jnode=node_child(2,pnode)
	    if(jnode.eq.idom) then
		jnode=inode
		inode=idom
	    end if
	    if(node_nseg(idom).eq.0) then
		if(pnode.eq.0) then
			write(*,*) 'no segments at root'
			goto 19
		else if(jnode.gt.0) then
			node_child(1,pnode)=node_child(1,jnode)
			node_child(2,pnode)=node_child(2,jnode)
			lactive(jnode)=.false.
		end if
	    end if
	  end if
	end do
19	continue
c
	do idom=1,ndom
		if(.not.lactive(idom).or.node_nseg(idom).eq.0) then
			node_parent(idom)=0
			node_child(1,idom)=0
			node_child(2,idom)=0
		end if
	end do
c
c	rebuild parent from children
c
	write(*,*) 'rebuild parent from children'
	do idom=1,ndom
		inode=node_child(1,idom)
		jnode=node_child(2,idom)
		if(inode.gt.0) node_parent(inode)=idom
		if(jnode.gt.0) node_parent(jnode)=idom
	end do
c
c	print tree
c
	do idom=1,ndom
	    write(*,499) idom,' ',node_parent(idom),node_child(1,idom),
     $			node_child(2,idom),node_nseg(idom),
     $			(node_range(idom,1,iseg),iseg=1,node_nseg(idom))
	end do
c
c	minimal tree top-down
c
	n=1
	list(1)=1
	tsil(1)=1
	list(0)=0
	tsil(0)=0
	do idom=2,ndom
		tsil(idom)=0
		if(node_parent(idom).gt.0) then
			n=n+1
			list(n)=idom
			tsil(idom)=n
		end if
	end do
	ndom=n
	do i=1,n
		idom=list(i)
		if(idom.gt.0) node_type(i)=node_type(idom)
		parent(i)=tsil(node_parent(idom))
		child(i,1)=tsil(node_child(1,idom))
		child(i,2)=tsil(node_child(2,idom))
		node_nseg(i)=node_nseg(idom)
c		segments are given as one at a time,so copy i,1,iseg to i,2,iseg
		do iseg=1,node_nseg(i)
			node_range(i,1,iseg)=node_range(idom,1,iseg)
			node_range(i,2,iseg)=node_range(i,1,iseg)
		end do
		node_nseg_o(i)=node_nseg_o(idom)
		node_size_o(i)=node_size_o(idom)
		do iseg=1,node_nseg_o(i)
			node_range_o(i,1,iseg)=node_range_o(idom,1,iseg)
			node_range_o(i,2,iseg)=node_range_o(idom,2,iseg)
		end do
	end do
	do i=1,n
		node_parent(i)=parent(i)
		node_child(1,i)=child(i,1)
		node_child(2,i)=child(i,2)
	end do
c
c	sort ranges, find blocks
c
	do i=1,n
		do iseg=1,nseg+1
			lactive(iseg)=.false.
		end do
		do j=1,node_nseg(i)
			iseg=node_range(i,1,j)
			lactive(iseg)=.true.
		end do
		k=0
		do iseg=1,nseg+1
			if(lactive(iseg)) then
				if((iseg.eq.1).or.
     $	  (iseg.gt.1.and..not.lactive(iseg-1))) then
					k=k+1
					node_range(i,1,k)=iseg
				end if
			end if
			if(.not.lactive(iseg)) then
				if(iseg.gt.1.and.lactive(iseg-1))
     $					node_range(i,2,k)=iseg-1
			end if
		end do
		node_size(i)=node_nseg(i)
		node_nseg(i)=k
	end do
c
c	print tree
c
	do idom=1,n
	    write(*,499) idom,node_type(idom),
     $	node_parent(idom),node_child(1,idom),
     $	node_child(2,idom),node_size(idom),node_nseg(idom),
     $	(node_range(idom,1,iseg),node_range(idom,2,iseg),
     $	iseg=1,node_nseg(idom))
	end do

499	format(i4,1x,a1,1x,100i3)

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine assignsegments(nseg,range,ndom,node_nseg,node_range,
     $          bot)
        implicit none
        include 'parsizes-puutos.for'
	integer nseg,range(2,maxseg),ndom,node_nseg(maxdom),bot(maxseg)
	integer node_range(maxdom,2,maxseg)
c
	integer ires,idom,iseg,bottom(maxres)
c
c	bottom node per each residue
c
	do idom=1,ndom
		do iseg=1,node_nseg(idom)
c			write(*,*) 'assign ',idom,iseg,
c     $				node_range(idom,1,iseg),node_range(idom,2,iseg)
	do ires=node_range(idom,1,iseg),node_range(idom,2,iseg)
		bottom(ires)=idom
	end do
		end do
	end do
c
c	segment midpoints -> bottom unit
c
	do iseg=1,nseg
		bot(iseg)=bottom((range(1,iseg)+range(2,iseg))/2)
c		write(*,*) 'assignsegments:',iseg,bot(iseg),nseg,
c     $			range(1,iseg),range(2,iseg)
	end do
c	write(*,*) (bottom(ires),ires=1,434)

	return
	end
c
c-------------------------------------------------------------------------------
c
       subroutine readnexttree(iunit,ndom,node_child,node_range,
     $    node_nseg,node_parent,node_size,code,chainid,node_type,ierr)
        implicit none
        include 'parsizes-puutos.for'
        integer iunit,ndom,node_child(2,maxdom)
        integer node_parent(maxdom),node_size(maxdom)
        real node_tau(maxdom)
        integer node_nseg(maxdom)
        character chainid,node_type(maxdom)
        character*4 code
        character*32 node_label(maxdom)
        real node_tert(maxdom)
        integer node_nbp(maxdom),ierr,node_range(maxdom,2,maxseg)
c
		integer i,j,k,idom
	character*80 line,separator
c
c	iunit must be opened outside
c       read next protein entry; skip until '>>>> '
c	returns ierr=1 on EOF; ierr=0 otherwise; ndom=0 on error
c
	ndom=0
10      read(iunit,500,err=10,end=99) line
c	write(*,*) line
	if(line(1:5).ne.'PDBID') goto 10
	read(line,600,end=99) code
	chainid=' '
20	read(iunit,500,end=99) line
	if(line(1:5).eq.'CHAIN') then
		read(line,610) chainid
		goto 20
	end if


c The following is the change made by Liisa to make puu3.f compilation work
c Aug 99.
		if(line(1:5).eq.'NUNIT') then
                read(line,620,end=99) ndom
                read(iunit,500,end=99) line
c                write(*,*) code,chainid,ndom
                do j=1,ndom
                 read(iunit,630,end=99) i,
     $ node_child(1,i),node_child(2,i),node_type(i),
     $ node_label(i),
     $ node_tau(i),node_tert(i),
     $ node_nbp(i),node_size(i),node_nseg(i)
				do k=1,node_nseg(i)
	read(iunit,631,end=99)  node_range(i,1,k),node_range(i,2,k)
				end do
		end do
                read(iunit,*) separator

		end if
c
c	rebuild node_parent
c
	do i=1,ndom
		node_parent(i)=0
	end do
	do i=1,ndom
		j=node_child(1,i)
		if(j.gt.0) node_parent(j)=i
		j=node_child(2,i)
		if(j.gt.0) node_parent(j)=i
	end do
c
c	print tree
c
c	write(*,*) 'imported tree'
c	do idom=1,ndom
c		write(*,640) idom,node_parent(idom),
c     $		  (node_child(i,idom),i=1,2),node_nseg(idom),
c     $         	  ((node_range(idom,i,j),i=1,2),j=1,node_nseg(idom))
c	end do
c
c	normal exit
c
	ierr=0
	return
c
500	format(a80)
600	format(10x,a4)
610	format(10x,a1)
620	format(10x,i5)
630     format(i5,1x,2i4,1x,a1,a32,2f9.3,i8,2i7,4x,100i4)
631	format(2i10)
640	format(7i5,100i4)
c
c	EOF exit
c
99	write(*,*) line
        ierr=1
	return
        end
c
c-------------------------------------------------------------------------------
c
	subroutine getdomseglist(domns,domseglist,ndom,node_nseg,
     $          node_range)
	implicit none
	include 'parsizes-puutos.for'
	integer domns(maxdom),domseglist(maxseg,maxdom),ndom
	integer node_range(maxdom,2,maxseg),node_nseg(maxdom)
c
	integer idom,ns,seglist(maxseg),is
c
	do idom=1,ndom
		call getseglist(ns,seglist,idom,node_nseg,node_range)
		domns(idom)=ns
		do is=1,ns
			domseglist(is,idom)=seglist(is)
		end do
	end do

	return
	end
c
c----------------------------------------------------------------------
c
        subroutine getsecstr(nres,struc,nseg,secstr,range,lene,lenh)
c
c       returns list of H and E segments only
c
        implicit none
        include 'parsizes-puutos.for'
        integer nres,nseg,range(2,maxseg),lene,lenh
        character struc(maxres),c,secstr(maxseg)
c
        integer i,j,k,l,n,ierr
c
        character simp
        character*5 fromstring,tostring
        data fromstring/'GIHE*'/
        data tostring  /'HHHEL'/
c
	if(nres.lt.1) then
		nseg=0
		return
	end if
        do i=1,nres
                struc(i)=simp(struc(i),5,fromstring,tostring)
        end do
        nseg=1
        range(1,nseg)=1
        c=struc(1)
        do i=2,nres
               if(struc(i).ne.c) then
                        range(2,nseg)=i-1
                        nseg=nseg+1
                        range(1,nseg)=i
                        c=struc(i)
                end if
        end do
        range(2,nseg)=nres
c
c       delete loop segments
c
        n=0
        do i=1,nseg
                k=range(1,i)
                l=range(2,i)
                if(struc(k).ne.'L') then
                        n=n+1
                        range(1,n)=k
                        range(2,n)=l
                end if
        end do
        nseg=n
c
c       delete short helix segments
c
        n=0
        do i=1,nseg
                k=range(1,i)
                l=range(2,i)
		c=struc(k)
                if((c.eq.'H').and.(l-k.lt.6)) then
			write(*,*) ' exclude short helix',k,l
		else if(c.eq.'H'.and.l-k.lt.lenh-1) then
			write(*,*) ' grow short helix',k,l
			ierr=0
			do while(l-k.lt.lenh-1.and.ierr.eq.0)
				call growleft(k,struc,ierr)
				call growrite(l,struc,nres,ierr)
			end do
			if(ierr.eq.0) then
                        	n=n+1
                        	range(1,n)=k
                        	range(2,n)=l
				secstr(n)=c
			else
				write(*,*) 'crowded segment excluded!'
			end if
		else if(c.eq.'E'.and.l-k.lt.lene-1) then
			write(*,*) ' grow short strand',k,l
			ierr=0
			do while(l-k.lt.lene-1.and.ierr.eq.0)
				call growleft(k,struc,ierr)
				call growrite(l,struc,nres,ierr)
			end do
			if(ierr.eq.0) then
  	    	                n=n+1
         	                range(1,n)=k
         	                range(2,n)=l
				secstr(n)=c
			else
				write(*,*) 'crowded segment excluded!'
			end if
		else if(n.eq.maxdom) then
	write(*,*) ' WARNING: too many segments ! skip at ',n
	goto 100
		else
                        n=n+1
                        range(1,n)=k
                        range(2,n)=l
    			secstr(n)=c
              end if
        end do
100     nseg=n
c        do i=1,nseg
c                k=range(1,i)
c                l=range(2,i)
c                write(*,*) i,k,l,'  >',(struc(j),j=k,l),'< ',secstr(i)
c        end do

        return
        end
c
c-----------------------------------------------------------------------
c
	subroutine growleft(k,struc,ierr)
        implicit none
        include 'parsizes-puutos.for'
        integer k,ierr
	character struc(maxres)
c
	character c
c
	ierr=-1
	if(k.eq.1) return
	c=struc(k-1)
	if(c.ne.'L') return
	k=k-1
	struc(k)=' '
	ierr=0

        return
        end
c
c-----------------------------------------------------------------------
c
	subroutine growrite(l,struc,nres,ierr)
        implicit none
        include 'parsizes-puutos.for'
        integer l,nres,ierr
	character struc(maxres)
c
	character c
c
	ierr=-1
	if(l.eq.nres) return
	c=struc(l+1)
	if(c.ne.'L') return
	l=l+1
	struc(l)=' '
	ierr=0

        return
        end
c
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c
	subroutine readdssp1(filnam,chainid,iunit,nres,struc,ca,seq,
     $          compnd,acc,resno)
c
c       frontend to getdssp to transport wanted subset of DSSP information
c       >>> uses 'parsizes-puutos.for' <<<
c
        implicit none
        include 'parsizes-puutos.for'
        integer iunit,nres,nchain,nss,nssintra,nssinter,nhb,acc(maxres)
        character*80 filnam,header,compnd,source,author
        character chainid,seq(maxres),struc(maxres)
        real nhbp,tco(maxres),kappa(maxres),alfa(maxres),phi(maxres),
     $          psi(maxres),ca(3,maxres)
        character*5 resno(0:maxres)
        character*33 dsspinfo(maxres)
        integer area,ierr,dsspresno(maxres),fsspresno(0:maxres)
c
        call getdssp(filnam,chainid,iunit,maxres,header,compnd,source,
     $          author,nres,nchain,nss,nssintra,nssinter,area,
     $          resno,seq,struc,acc,dsspinfo,dsspresno,fsspresno,
     $          tco,kappa,alfa,phi,psi,ca,ierr)
        write(*,*) header,compnd,source,author
        write(*,*) ierr
        write(*,*) nres
        if(ierr.eq.0) write(*,*) nres,' residues read ',filnam(1:60)
        if(ierr.eq.-1) write(*,*) ' DSSP file not found ',filnam(1:60)
        if(ierr.eq.-2) write(*,*) ' DSSP file truncated at ',
     $		maxres,'residues'

        return
        end
c
c-----------------------------------------------------------------------
c
	function simp(c,len,fromstring,tostring)
c
c	converts character c in fromstring to its equivalent in tostring
c	len-th is wildcard
c
	implicit none
	integer len
	character*(*) fromstring,tostring
	character simp,c
c
	integer i
c
	simp=tostring(len:len)
	do i=1,len-1
		if(c.eq.fromstring(i:i)) then
			simp=tostring(i:i)
			goto 99
		end if
	end do

99	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine getseglist(ns,seglist,inode,node_nseg,node_range)
	implicit none
	include 'parsizes-puutos.for'
	integer ns,seglist(maxseg),inode,node_nseg(maxdom)
	integer node_range(maxdom,2,maxseg)
c
	integer j,k
c
	ns=0
	if(inode.gt.0) then
		do j=1,node_nseg(inode)
			do k=node_range(inode,1,j),node_range(inode,2,j)
				ns=ns+1
				seglist(ns)=k
			end do
		end do
	end if

	return
	end
c
c----------------------------------------------------------------------
c
c
	function countchr(c,secstr,nseg)
	implicit none
	integer countchr,nseg
	character secstr(nseg)
	character c
c
	integer i,n

	n=0
	do i=1,nseg
		if(secstr(i).eq.c) n=n+1
	end do
	countchr=n

	return
	end
c
c----------------------------------------------------------------------
c

