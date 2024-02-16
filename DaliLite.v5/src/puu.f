c This module/program is part of DaliLite (c) L. Holm 1999
c
CCCCCCCC puu.f - standalone version 18-Aug-1993 -- linearized 25-Jul-1994
	subroutine puu1(cd1,brkfile,dsspfile)
	implicit none
	include 'sizes-puu.for'
        character(len=5) cd1
        character(len=80) brkfile,dsspfile
c
c	f77 main.f puu.f
c	(main.f = call puu)
c	on GOLD link with ~/dccp/rubber/seed/GOL/ran.o
c
	integer nsearch,search(3,1000)
        character*80 dssplist,pdblist,defaultsfilename
        character*4 code
        character chainid
c
	integer minseglen,printdepth
	real breakdist,hbond,compact,taucutoff
	character*80 treefilename,prettyfilename,subunitfilename
	character*60 basedefault,basestring,line,emptystring,
     $          versionstring
	integer*4 seed,i
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
c
c
	versionstring='2.1 (February 2016)' ! first round = chain-subunits
c
c
c
	write(*,*) ' PUU: Parser for Protein Unfolding Units '
	write(*,*) ' Ref: L. Holm, C. Sander, Proteins 19:256-268, 1994'
	write(*,*) ' Ver: ',versionstring
	write(*,*) ' Dimensioned for ',maxres,' residues.'
	write(*,*) ' '
c
c
c
        defaultsfilename='puu.default'
	call defaults(dssplist,pdblist,breakdist,minseglen,hbond,
     $  compact,taucutoff,seed,treefilename,prettyfilename,
     $  defaultsfilename,lhetatm,lsplit,printdepth,
     $  basedefault,subunitfilename)
c
c	initialization
c
 	call init_search(search,nsearch,4.0)
c
c	open output files
c
	open(15,file=treefilename)
	open(32,file=prettyfilename)
	open(33,file=subunitfilename)
c
c	write headers
c
	call writeheadertobase(15,
     $	'Protein Unfolding Trees',versionstring,
     $		basestring)
	call writeparams(15,breakdist,minseglen,hbond,compact,
     $		taucutoff,seed,printdepth,lsplit,lhetatm,.true.)
	call writeunitsnotation(15)

	call writeheadertobase(32,
     $	'Protein Structural Domains',versionstring,
     $		basestring)
	call writeparams(32,breakdist,minseglen,hbond,compact,
     $		taucutoff,seed,0,lsplit,lhetatm,.false.)
	call writedomainnotation(32)

	call writeheadertobase(33,
     $	'Protein Subunits',versionstring,basestring)
	write(33,530) '//'

c
c	input loop
c
c10      write(*,*) ' enter code+chainid ? (END to quit) '
c        read(*,500) code,chainid
        code=cd1(1:4)
        chainid=cd1(5:5)
        if(code(1:3).eq.'END') goto 999
	call doprotein(code,chainid,nsearch,search,dssplist,lhetatm,
     $  pdblist,
     $	breakdist,minseglen,seed,hbond,compact,taucutoff,lsplit,
     $	printdepth,brkfile,dsspfile)
c	goto 10
c
500	format(a4,a1)
510	format(a80)
520	format(a60)
530	format(a2)
c
999	close(15)
	close(32)
	close(33)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine doprotein(code,chainid,nsearch,search,dssplist,
     $  lhetatm,
     $	pdblist,breakdist,minseglen,seed,hbond,compact,taucutoff,
     $	lsplit,printdepth,brkfile,dsspfile)
	implicit none
	include 'sizes-puu.for'
	integer nsearch,search(3,1000),printdepth
	real ca0(3,maxres),cont(500,500),compact,taucutoff
        character*80 dssplist,pdblist,confile,defaultsfilename
        character*4 code
        character(len=80) brkfile,dsspfile
        character chainid
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
	integer nres0,ndom,natom0(maxres)
	integer*4 seed
	real hbond
c
	integer minseglen
	real breakdist
	real node_tau(maxres)
	integer node_child(maxres,2),node_nseg(maxres)
	integer node_range(maxres,2,maxnseg),node_nbp(maxres)
	integer bp0(2,maxres)
	integer node_parent(maxres),node_size(maxres)
	real node_tert(maxres)
c
	integer npoint0,point0(2,maxres*20)
	real occx0(maxres*20)
c
	integer nstack,stack(maxres),set(maxres)
	character node_type(maxres)
c
	integer i,j,k,n
	character c,struc(maxres)
	character*5 resno(0:maxres)
	character*80 protname
c
	logical lonechain
	integer nchain,seg(2,maxres),restoseg(maxres),natom(maxres)
	character oldc
	real sscont(100,100)
	integer node_r(maxres,2,maxnseg),tmp(maxres)
	real ca01(3,maxres)
	integer npoint01,point01(2,maxres*20),minseglen0
	real occx01(maxres*20)
c
	character chain(maxres)
	integer subunit(maxres)
        logical lwritecontacts
        lwritecontacts=.false.
c
c	read protein
c
	nres0=0
	call setup_protein(code,chainid,nres0,natom0,occx0,npoint0,
     $  point0,
     $	ca0,nsearch,search,dssplist,pdblist,bp0,hbond,resno,struc,
     $	protname,lhetatm,brkfile,dsspfile)

        if(lwritecontacts) then
c       write contacts to file
                write(*,*) 'enter name of contact file (output)'
                read(*,510) confile
510             format(a80)
                open(90,file=confile)
	        do i=1,npoint0
		write(90,500) code,chainid,point0(1,i),point0(2,i),
     $                  int(occx0(i))
	        end do
                close(90)
        end if
c
c	>>> if multichain then collapse chains <<<
c
c	set minseglen=1 !
c
	minseglen0=minseglen
	minseglen=1
	i=1
	lonechain=.true.
	do while(i.lt.nres0.and.lonechain)
		i=i+1
		lonechain=(resno(i)(1:1).eq.resno(i-1)(1:1))
	end do
	if(.not.lonechain) then
		nchain=1
		oldc=resno(1)(1:1)
		chain(1)=oldc
		seg(1,1)=1
		do i=2,nres0
			c=resno(i)(1:1)
			if(c.ne.oldc) then
				seg(2,nchain)=i-1
				nchain=nchain+1
				seg(1,nchain)=i
				chain(nchain)=c
				subunit(nchain)=0
				oldc=c
			end if
		end do
		seg(2,nchain)=nres0
c		type *,'nchain=',nchain
		do i=1,nchain
c			type *,resno(seg(1,i))(1:1),seg(1,i),seg(2,i)
			natom(i)=0
		end do
		do i=1,nchain
			do j=seg(1,i),seg(2,i)
				restoseg(j)=i
				natom(i)=natom(i)+natom0(j)
			end do
			do j=1,nchain
				sscont(i,j)=0.0
			end do
		end do
		do i=1,npoint0
			j=restoseg(point0(1,i))
			k=restoseg(point0(2,i))
			sscont(j,k)=sscont(j,k)+occx0(i)
		end do
c		do i=1,nchain
c			type *,(sscont(i,j),j=1,nchain)
c		end do
		npoint01=0
		do i=1,nchain
			do j=1,nchain
				npoint01=npoint01+1
				point01(1,npoint01)=i
				point01(2,npoint01)=j
				occx01(npoint01)=sscont(i,j)
			end do
			do j=1,3
				ca01(j,i)=0.0
			end do
		end do
c
c		!! first round: decompose with "chain-residues"
c
		call initialize_one(nchain,ndom,nstack,stack,set,
     $	node_range,node_nseg,node_child,node_tau)
		call dostack(nchain,ca01,occx01,npoint01,point01,natom,
     $			breakdist,
     $			minseglen,seed,node_range,node_nseg,node_child,
     $			node_tau,ndom,cont,nstack,stack,set)
		call patch(ndom,nchain,node_child,node_parent,node_nseg,
     $			node_range,node_size)
		do i=1,ndom
			do j=1,node_nseg(i)
				node_r(i,1,j)=seg(1,node_range(i,1,j))
				node_r(i,2,j)=seg(2,node_range(i,2,j))
			end do
c	write(*,*) 'node_r:',i,((node_r(i,k,j),k=1,2),j=1,node_nseg(i))
		end do
		call contstren(occx0,npoint0,point0,node_r,node_nseg,
     $			natom0,ndom,nres0,node_tert,cont,4)
c		write(*,*) 'tert:',(node_tert(i),i=1,ndom)
c		write(*,*) 'natom:',(natom(i),i=1,nchain)
		call betasheets(ndom,node_child,node_nbp,bp0,node_r,
     $			node_nseg)
		call findstrucdoms(ndom,node_child,node_size,node_tau,
     $			node_tert,node_nbp,compact,taucutoff,lsplit,
     $			node_type,0)
c
		write(*,*) 'type:',(node_type(i),i=1,ndom)
c
c		remember subunits (inode) of each chain !
c
		do i=1,ndom
			if(node_type(i).eq.'*') then
				do j=1,node_nseg(i)
	  do k=node_range(i,1,j),node_range(i,2,j)
		subunit(k)=i
	  end do
				end do
			end if
		end do
c
c		!! second round: initialize stack with subunits
c
		do i=1,maxres
			set(i)=0
		end do
		n=0
		nstack=0
		do i=1,ndom
		  if(node_type(i).eq.'*') then ! copy and stack
			n=n+1
			do j=1,node_nseg(i)
				do k=node_r(i,1,j),node_r(i,2,j)
					set(k)=n
				end do
			end do
			call addrequest(n,nstack,stack)
			node_nseg(n)=node_nseg(i)
			node_tau(n)=node_tau(i)
			node_child(n,1)=0
			node_child(n,2)=0
			do j=1,node_nseg(i)
				node_range(n,1,j)=node_r(i,1,j)
				node_range(n,2,j)=node_r(i,2,j)
			end do
			tmp(i)=n
		  else if(node_type(i).eq.'+') then
			n=n+1
			tmp(i)=n
		  end if
		end do
		k=n
		do i=1,ndom
		  if(node_type(i).eq.'+') then ! copy but do not stack
			n=tmp(i)
			node_nseg(n)=node_nseg(i)
			node_tau(n)=node_tau(i)
			node_child(n,1)=tmp(node_child(i,1))
			node_child(n,2)=tmp(node_child(i,2))
			do j=1,node_nseg(i)
				node_range(n,1,j)=node_r(i,1,j)
				node_range(n,2,j)=node_r(i,2,j)
			end do
		  end if
		end do
		ndom=k
c		write(*,*) 'nstack',nstack,(stack(i),i=1,nstack)
c		write(*,*) (set(i),i=1,nres0)
	else
		nchain=1
		chain(1)=chainid
		subunit(1)=1
		call initialize_one(nres0,ndom,nstack,stack,set,
     $		 node_range,node_nseg,node_child,node_tau)
	end if
c
c	output subunit decomposition
c
        write(33,551) code
        write(33,552) nchain
        write(33,553) (subunit(i),i=1,nchain)
        write(33,554) (chain(i),i=1,nchain)
        write(33,555) '//'
c
c	hack
c
	minseglen=minseglen0
c
	call dostack(nres0,ca0,occx0,npoint0,point0,natom0,breakdist,
     $	minseglen,seed,node_range,node_nseg,node_child,node_tau,
     $	ndom,cont,nstack,stack,set)
	write(*,*) 'dostack done'
c
c###	process further: patch_nodes, tertcon, patch_solp
c
	call patch(ndom,nres0,node_child,node_parent,node_nseg,
     $	 node_range,node_size)
c	write(*,*) 'patch done'
c
c	tertiary contacts
c
	call contstren(occx0,npoint0,point0,node_range,node_nseg,natom0,
     $		ndom,nres0,node_tert,cont,4)
c	write(*,*) 'constren done'
c
c	beta-sheets
c
	call betasheets(ndom,node_child,node_nbp,bp0,node_range,
     $          node_nseg)
c	write(*,*) 'beta-sheets done'
c
	call findstrucdoms(ndom,node_child,node_size,node_tau,
     $	 node_tert,node_nbp,compact,taucutoff,lsplit,node_type,80)
	write(*,*) 'findstrucdoms done'
c
c###	write some output !
c
	call pret2(ndom,node_tau,node_child,node_range,node_nseg,
     $		node_parent,node_size,code,chainid,
     $		node_tert,node_nbp,resno,struc,compact,
     $		protname,printdepth,node_type)
        write(*,*) 'pret2 done'
c
500     format(a4,a1,1x,3i5)
520	format(a4,a1,' best cut was ',i5,f10.3,10i5)
530	format(20i4)
540	format(i4,f8.3,20i4)
551     format('PDBID     ',a4)
552     format('NCHAIN    ',i5)
553     format('UNITID    ',i5)
554     format('CHAINID   ',300(4x,a1))
555     format(a2)
600	format('>>>> ',a4,a1,i10)
610     format(5i4,2x,a32,f8.3,f8.1,f8.3,20i4)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine betasheets(ndom,node_child,node_nbp,bp0,node_range,
     $		node_nseg)
	implicit none
	include 'sizes-puu.for'
	integer ndom
	integer node_child(maxres,2),node_nseg(maxres),node_nbp(maxres)
	integer node_range(maxres,2,maxnseg)
	integer bp0(2,maxres)
	integer i,j,k,nsplitsheet
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
	do i=1,ndom
		j=node_child(i,1)
		k=node_child(i,2)
		if((j.gt.0).and.(k.gt.0)) then
	node_nbp(i)=nsplitsheet(j,k,bp0,node_range,node_nseg)
		else
			node_nbp(i)=0
		end if
	end do

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine labels(code,chainid,ndom,node_child,node_label)
	implicit none
	include 'sizes-puu.for'
        character*4 code
        character chainid
	integer ndom
	integer node_child(maxres,2)
	character*32 node_label(maxres)
	integer node_depth(maxres)
c
	integer i,j,m,depth
	character c
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
c	generate label codeA.1.1.2.2.1.2.1
c
	m=depth(ndom,node_child,node_depth)
c!	write(*,*) ' depth of tree is ',m
	c=chainid
	if(c.eq.' ') c='_'
	do i=1,ndom
		node_label(i)(1:4)=code
		node_label(i)(5:5)=c
		do j=6,32
			node_label(i)(j:j)=' '
		end do
	end do
	do i=1,ndom
	j=node_child(i,1)
	if(j.gt.0)call extendlabel('1',j,i,node_label,node_depth(j))
	j=node_child(i,2)
	if(j.gt.0) call extendlabel('2',j,i,node_label,node_depth(j))
	end do
c	write(*,*) 'labels done'
c	do i=1,ndom
c	  write(*,500) i,node_depth(i),(node_child(i,j),j=1,2),node_label(i)
c	end do
c500	format(4i5,2x,a32)

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine initialize_one(nres0,ndom,nstack,stack,set,
     $   node_range,node_nseg,node_child,node_tau)
	implicit none
	include 'sizes-puu.for'
	integer nres0,ndom,set(maxres)
c
	integer nstack,stack(maxres)
	real node_tau(maxres)
	integer node_child(maxres,2),node_nseg(maxres)
	integer node_range(maxres,2,maxnseg)
c
	integer i
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
c	initialize search with protein
c
	do i=1,maxres
		if(i.le.nres0) set(i)=1
		if(i.gt.nres0) set(i)=0
	end do
	ndom=1
	nstack=0
c
c	put set 1 in stack
c
	call addrequest(1,nstack,stack)
	node_range(1,1,1)=min(1,nres0)
	node_range(1,2,1)=nres0
	node_nseg(1)=1
	node_child(1,1)=0
	node_child(1,2)=0
	node_tau(1)=0.0

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine dostack(nres0,ca0,occx0,npoint0,point0,natom0,
     $  breakdist,
     $	minseglen,seed,node_range,node_nseg,node_child,node_tau,
     $	ndom,cont,nstack,stack,set)
	implicit none
	include 'sizes-puu.for'
	real ca0(3,maxres),ca(3,maxres),cont(500,500)
	integer nres,nres0,ndom,set(maxres),natom0(maxres), n(maxres)
	logical ltest(maxres)
	integer*4 seed
c
	integer bottom,bestcut(maxres),nstack,stack(maxres),iset
	integer resix(maxres),minseglen
	real taubest,breakdist
	real node_tau(maxres)
	integer node_child(maxres,2),node_nseg(maxres)
	integer node_range(maxres,2,maxnseg)
c
	integer npoint0,point0(2,maxres*20),npoint,point(2,maxres*20)
	real occx0(maxres*20),occx(maxres*20)
	integer i
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
c
c	binary decomposition
c
        write(*,*) '# This is dostack',nstack
	do while(nstack.gt.0)
c
c	  copy selected set to work arrays: ca, nres, n, cont
c
	  call toprequest(iset,nstack,stack)
          write(*,*) '#popped',iset,nstack
	  call selectset(iset,set,ca0,occx0,npoint0,point0,occx,npoint,
     $		point,nres0,natom0,ca,cont,nres,n,ltest,
     $		resix,breakdist,minseglen)
	  write(*,*) 'selectset done',nstack
c
c	  find best cut for selected set
c
	  taubest=0.0
	  do i=1,maxres
		bestcut(i)=0
	  end do
	  call eigen(bestcut,bottom,taubest,nres,n,occx,npoint,point,
     $		cont,ltest,minseglen,seed)
	  write(*,*) 'eigen done',nstack
	  write(*,*) 'tau',ndom,taubest,(resix(bestcut(i)),i=1,bottom)
c
c	  remember domains: quit if taubest=0, i.e. no subdivision was found
c
	  if(taubest.gt.0.0) call savecut(iset,bestcut,bottom,nres,
     $  set,resix,
     $	taubest,ndom,nstack,stack,node_tau,node_child,node_range,
     $	node_nseg,minseglen)
	end do

	return
	end
c
c-----------------------------------------------------------------------
c
	function nsplitsheet(inode,jnode,bp0,node_range,node_nseg)
	implicit none
	include 'sizes-puu.for'
	integer nsplitsheet
	integer inode,jnode,node_range(maxres,2,maxnseg)
	integer bp0(2,maxres),node_nseg(maxres)
c
	integer i,j,k1,k2,idom(maxres),n
c
c	count no of non-local inter-domain H-bonds
c
	do i=1,maxres
		idom(i)=0
	end do
	do i=1,node_nseg(inode)
		do j=node_range(inode,1,i),node_range(inode,2,i)
			idom(j)=1
		end do
	end do
	n=0
	do i=1,node_nseg(jnode)
		do j=node_range(jnode,1,i),node_range(jnode,2,i)
			k1=bp0(1,j)
			k2=bp0(2,j)
			if((k1.gt.0).and.(k2.gt.0)) then
				if(idom(k1).eq.1) n=n+1
				if(idom(k2).eq.1) n=n+1
			end if
		end do
	end do
c	write(*,*) 'inter-H-bonds: ',n
	nsplitsheet=n
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine fastcont(nres,occx,npoint,point,cont)
	implicit none
	include 'sizes-puu.for'
	integer nres,npoint,point(2,maxres*20)
	real occx(maxres*20),cont(500,500)
c
	integer ires,jres,i
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
	if(nres.le.500) then
		do ires=1,nres
			do jres=1,nres
				cont(ires,jres)=0.0
			end do
		end do
		do i=1,npoint
			ires=point(1,i)
			jres=point(2,i)
			cont(ires,jres)=cont(ires,jres)+occx(i)
		end do
		do ires=1,nres
			do jres=1,ires-1
	 if(cont(ires,jres).ne.cont(jres,ires))write(*,*)
     $	 'mismatching:',ires,jres,cont(ires,jres),cont(jres,ires)
			end do
		end do
		do ires=1,nres
			do jres=2,ires
	cont(ires,jres)=cont(ires,jres-1)+cont(ires,jres)
			end do
		end do
	else
	write(*,*) 'ERROR: fastcont called with nres=',nres
	end if

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine contstren(occx,npoint,point,node_range,node_nseg,
     $	natom0,ndom,nres,node_tert,cont,four)
	implicit none
	include 'sizes-puu.for'
	integer natom0(maxres),node_range(maxres,2,maxnseg)
	integer ndom,nres,node_nseg(maxres)
	real occx(maxres*20)
	integer npoint,point(2,maxres*20),four
	real node_tert(maxres)
	real cont(500,500)
c
	integer i,j,k,l,iseg,natom,n(maxres),ires,jres,iset,jset
        integer set(maxres),ix
	real x
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
	n(1)=natom0(1)
	do i=2,nres
		n(i)=n(i-1)+natom0(i)
	end do
c
c	!!!!use fastcont for small proteins
c
	if(nres.le.500) call fastcont(nres,occx,npoint,point,cont)
	do i=1,ndom
c
c		# atoms in domain
c
		natom=0
		do iseg=1,node_nseg(i)
	natom=natom+n(node_range(i,2,iseg))-n(node_range(i,1,iseg))
		end do
		if(natom.eq.0) natom=1
c
c		calculate tertiary intra-domain contact strength
c
		if(nres.le.500) then
		  x=1e-9
		  do iseg=1,node_nseg(i)
			k=node_range(i,1,iseg)
			l=node_range(i,2,iseg)
			do j=k+four,l
				x=x+cont(j,j-four)
				if(k.gt.1) x=x-cont(j,k-1)
			end do
		  end do
		else
		  do ires=1,nres
			set(ires)=2
		  end do
		  do iseg=1,node_nseg(i)
	do ires=node_range(i,1,iseg),node_range(i,2,iseg)
		set(ires)=1
	end do
		  end do
		  x=1e-9
		  do ix=1,npoint
			ires=point(1,ix)
			jres=point(2,ix)
			if(abs(ires-jres).ge.four) then
				iset=set(ires)
				jset=set(jres)
	if((iset.eq.1).and.(jset.eq.1)) x=x+occx(ix)
			end if
		  end do
		end if
		node_tert(i)=x/natom*100.0
c		write(*,*) 'tert ',i,x,natom
	end do
c
	return
	end
c
c-----------------------------------------------------------------------------
c
	function depth(ndom,node_child,d)
	implicit none
	include 'sizes-puu.for'
	integer depth,ndom,node_child(maxres,2),d(maxres)
c
	integer i,j,k,t
c
	d(1)=1
	do i=2,ndom
		d(i)=0
	end do
	do i=1,ndom
	  do t=1,2
		j=node_child(i,t)
		if(j.gt.0) then
			d(j)=d(i)+1
		end if
	  end do
	end do
	k=0
	do i=1,ndom
		if(d(i).gt.k) k=d(i)
	end do
	depth=k
c
	end
c
c-------------------------------------------------------------------------------
c
	subroutine extendlabel(c,inode,pnode,node_label,depth)
	implicit none
	include 'sizes-puu.for'
	integer inode,pnode,depth
	character*32 node_label(maxres)
	character c
        logical lhetatm,lsplit,lsource,lauthor,lcompnd,lheader
c
	integer i
c
	i=5+(depth-2)*2+1
c
c	error exit
c
	if(i.gt.32) then
		write(*,*) ' label overflow'
		return
	end if
c
c	normal exit
c
	if(i.gt.5) node_label(inode)(6:i)=node_label(pnode)(6:i)
	node_label(inode)(i:i)='.'
	i=i+1
	node_label(inode)(i:i)=c
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine  patch(ndom,nres,node_child,node_parent,node_nseg,
     $		node_range,node_size)
	implicit none
	include 'sizes-puu.for'
	integer ndom,node_parent(maxres),node_nseg(maxres),nres
	integer node_range(maxres,2,maxnseg),node_size(maxres)
	integer node_child(maxres,2)
c
	integer i,resix(maxres),j,k,p
c
	do i=1,ndom
	  do k=1,2
		j=node_child(i,k)
		if(j.gt.0) node_parent(j)=i
	  end do
	end do
	do i=1,nres
		resix(i)=1
	end do
	do i=2,ndom
		p=node_parent(i)
		do j=1,node_nseg(i)
			do k=node_range(i,1,j),node_range(i,2,j)
				if(resix(k).eq.p) resix(k)=i
			end do
		end do
		p=0
		if(resix(1).eq.i) then
			p=1
			node_range(i,1,p)=1
		end if
		do j=1,nres-1
		if((resix(j).ne.i).and.(resix(j+1).eq.i)) then
			p=p+1
			if(p.gt.maxnseg) then
				write(*,*) 'ERROR: nseg overflow'
				p=maxnseg
			end if
			node_range(i,1,p)=j+1
		end if
		if((resix(j).eq.i).and.(resix(j+1).ne.i)) then
			node_range(i,2,p)=j
		end if
        	end do
		if(resix(nres).eq.i) node_range(i,2,p)=nres
		node_nseg(i)=p
	end do
c
c	remember size of node
c
	node_size(1)=nres
	do i=2,ndom
		k=0
		do j=1,node_nseg(i)
			k=k+node_range(i,2,j)-node_range(i,1,j)+1
		end do
		node_size(i)=k
	end do
c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine savecut(iset,bestcut,bottom,nres,set,resix,taubest,
     $		ndom,nstack,stack,node_tau,node_child,node_range,
     $		node_nseg,minseglen)
	implicit none
	include 'sizes-puu.for'
	integer bestcut(maxres),bottom,nres,set(maxres),resix(maxres)
	integer ndom,nstack,stack(maxres),iset,minseglen
	real node_tau(maxres)
	integer node_child(maxres,2)
	integer node_range(maxres,2,maxnseg),node_nseg(maxres)
	real taubest
c
	integer iseg,i,j,n1,n2,nseg
c
c	fill tau, child records for parent node
c
	node_tau(iset)=taubest
	node_child(iset,1)=ndom+1
	node_child(iset,2)=ndom+2
	bestcut(bottom+1)=nres
c
c	child-2 node
c		update set(), stack request
c		fill range record
c
	n2=0
	nseg=0
	do iseg=2,bottom+1,2
		j=bestcut(iseg-1)+1
		do i=j,bestcut(iseg)
			set(resix(i))=ndom+2
			n2=n2+1
		end do
		nseg=nseg+1
		node_range(ndom+2,1,nseg)=resix(j)
		node_range(ndom+2,2,nseg)=resix(bestcut(iseg))
	end do
	node_child(ndom+2,1)=0
	node_child(ndom+2,2)=0
	node_tau(ndom+2)=0.0
	node_nseg(ndom+2)=nseg
c	write(*,*) 'save2:',n2,(bestcut(iseg),iseg=2,bottom+1,2)
	if(n2.ge.2*max(1,minseglen)) call addrequest(ndom+2,nstack,stack)
c
c	child-1 node
c		update set(), stack request
c		fill range record
c
	n1=0
	nseg=0
	do iseg=1,bottom+1,2
		if(iseg.eq.1) then
			j=1
		else
			j=bestcut(iseg-1)+1
		end if
		do i=j,bestcut(iseg)
			set(resix(i))=ndom+1
			n1=n1+1
		end do
		nseg=nseg+1
		node_range(ndom+1,1,nseg)=resix(j)
		node_range(ndom+1,2,nseg)=resix(bestcut(iseg))
	end do
	node_child(ndom+1,1)=0
	node_child(ndom+1,2)=0
	node_tau(ndom+1)=0.0
	node_nseg(ndom+1)=nseg
c	write(*,*) 'save1:',n1,(bestcut(iseg),iseg=1,bottom+1,2)
	if(n1.ge.2*max(1,minseglen)) call addrequest(ndom+1,nstack,stack)
c
c	update ndom
c
	ndom=ndom+2
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine eigen(bestcut,bottom,taubest,nres,n,occx,npoint,point,
     $		cont,ltest,minseglen,seed)
	implicit none
	include 'sizes-puu.for'
	real taubest,cont(500,500)
	logical ltest(maxres)
	integer*4 seed
	integer bestcut(maxres),bottom,nres,n(maxres),minseglen
	integer npoint,point(2,maxres*20)
	real occx(maxres*20)
c
	integer i,k,icut(maxres+1),ires,bestres,ei,choice
	real x(3,maxres),tau,t,r2(3),tauslow
	integer eigenvalues
c
c	!!!! call with shortlist
c
	call correspondence_shortlist(nres,nres,occx,point,npoint,x,r2,
     $          seed,2)
c
c	!!! ordinate sequence, do 1D sweep for best cutting-point
c
	bottom=0
	bestres=0
	ei=2
c	write(*,*) 'eigenvector:',(x(ei,ires),ires=1,nres)
	do ires=1,nres
	  if(ltest(ires)) then
	    do choice=1,2
		t=x(ei,ires)
	call definesets(x,nres,icut,k,minseglen,ltest,t,ei,choice)
c
c		calculate score for selected set in icut(1...k)
c
		if(k.gt.0) then
			icut(k+1)=nres
			if(nres.le.500) then
				t=tau(icut,k+1,nres,cont,n)
			else
		t=tauslow(icut,k+1,nres,occx,npoint,point,n)
			end if
			if(t.gt.taubest) then
				bottom=k
				taubest=t
				bestres=ires
				do i=1,k
					bestcut(i)=icut(i)
				end do
			end if
		end if
c		write(*,*) 'sets',ires,(icut(i),i=1,k),choice,k,t
	    end do
	  end if
	end do
c	if(bestres.gt.0) then
c		write(*,*) 'eigen:',bestres,taubest,nres,x(2,bestres)
c		write(*,500) (bestcut(i),i=1,bottom)
c	else
c		write(*,*) 'eigen: bestres=0'
c	end if
c
500	format(20i4)
520	format(2i5,20(f8.4))
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine addrequest(i,nstack,stack)
	implicit none
	include 'sizes-puu.for'
	integer i,nstack,stack(maxres)
c
	if(nstack.lt.maxres) then
		nstack=nstack+1
		stack(nstack)=i
	else
		stop ' stack overflow'
	end if
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine toprequest(i,nstack,stack)
	implicit none
	include 'sizes-puu.for'
	integer i,nstack,stack(maxres)
c
	i=stack(nstack)
	nstack=nstack-1
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine definesets(x,nres,icut,k,minseglen,ltest,zerolevel,
     $		ei,choice)
	implicit none
	include 'sizes-puu.for'
	integer k,icut(maxres),nres,minseglen,ei,choice
	real x(3,maxres),zerolevel
	logical ltest(maxres)
c
	integer i,n,j,l,set(maxres)
c
c	determine segment ends as where y changes sign
c
	do i=1,nres
		set(i)=1
	end do
	if(choice.eq.1) then
		do i=1,nres
			if(x(ei,i).lt.zerolevel) set(i)=2
		end do
	else
		do i=1,nres
			if(x(ei,i).gt.zerolevel) set(i)=2
		end do
	end if
	n=0
	do i=2,nres
		if(set(i).ne.set(i-1)) then
			n=n+1
			icut(n)=i-1
		end if
	end do
c	last cut added anyway
	if(icut(n).eq.nres) n=n-1
c	if(n.eq.0) write(*,500) zerolevel,(x(ei,i),i=1,nres)
c500	format('zero    ',9f8.3,50(/10f8.3))
c	write(*,*) n,'defined1:',(icut(i),i=1,n)
c
c	remove too short loops, and cuts at a forbidden residue
c
	do i=1,n
c		distance to previous active cut
		l=icut(i)
		if(i.gt.1) then
			j=i-1
			do while(j.gt.0)
				if(icut(j).eq.0) then
					j=j-1
				else
					goto 10
				end if
			end do
10			if(j.gt.0) l=icut(i)-icut(j)
		end if
		if((l.lt.minseglen).or.(.not.ltest(icut(i)))) then
			icut(i)=0
			if(i.gt.1) icut(i-1)=0
		end if
	end do
c	write(*,*) n,'defined:',(icut(i),i=1,n)
c
c	return clean list
c
	k=0
	do i=1,n
		if(icut(i).gt.0) then
			k=k+1
			icut(k)=icut(i)
		end if
	end do
	if(k.eq.1.and.icut(k).eq.nres) k=0
c
c
c
c	write(*,*) k,'defined2:',(icut(i),i=1,k)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	function tau(icut,bottom,nres,cont,n)
c
c	returns tau-squared = reduced mass / interface contact strength
c
	implicit none
	include 'sizes-puu.for'
	real tau,cont(500,500)
	integer icut(maxres+1),bottom,n(maxres),nres
c
	integer iseg,jseg,i
	real x,m,taux
c
c	contact sum over interface
c
	x=0.0
	do iseg=1,bottom
		do jseg=iseg+1,bottom,2
			do i=icut(jseg-1)+1,icut(jseg)
				x=x+cont(i,icut(iseg))
				if(iseg.gt.1) x=x-cont(i,icut(iseg-1))
			end do
		end do
	end do
	x=max(x,1e-9)
c
c	# atoms in domain 2
c
	m=0.0
	do jseg=2,bottom,2
		m=m+n(icut(jseg))-n(icut(jseg-1))
	end do
c
c	get tau-squared
c
c	type *,'tau',nres,n(nres),m,x
	taux=m*(n(nres)-m)/max(1,n(nres))/x
	tau=min(taux,999.999)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	function tauslow(icut,bottom,nres,occx,npoint,point,n)
c
c	returns tau-squared = reduced mass / interface contact strength
c
	implicit none
	include 'sizes-puu.for'
	real tauslow,occx(maxres*20)
	integer npoint,point(2,maxres*20)
	integer icut(maxres+1),bottom,n(maxres),nres
c
	integer iseg,jseg,i,iset,jset,set(maxres),ires
	real x,m,tauslowx
c
c	two sets
c
	do ires=1,nres
		set(ires)=1
	end do
	do iseg=2,bottom,2
		do ires=icut(iseg-1)+1,icut(iseg)
			set(ires)=2
		end do
	end do
c
c	contact sum over interface
c
	x=0.0
	do i=1,npoint
		iset=set(point(1,i))
		jset=set(point(2,i))
		if(iset.ne.jset) x=x+occx(i)
	end do
	x=max(x,1e-9)
c
c	# atoms in domain 2
c
	m=0.0
	do jseg=2,bottom,2
		m=m+n(icut(jseg))-n(icut(jseg-1))
	end do
c
c	get tau-squared
c
c	type *,'tauslow',nres,n(nres),m,x
	tauslowx=m*(n(nres)-m)/max(1,n(nres))/x
	tauslow=min(tauslowx,999.999)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine selectset(iset,set,ca0,
     $		occx0,npoint0,point0,occx,npoint,point,
     $		nres0,natom0,ca,cont,nres,n,
     $		ltest,resix,breakdist,minseglen)
	implicit none
	include 'sizes-puu.for'
	integer iset,set(maxres),nres0,nres,natom0(maxres),n(maxres)
	integer resix(maxres),minseglen
	real ca0(3,maxres),ca(3,maxres)
	real cont(500,500),breakdist
	logical ltest(maxres)
	integer npoint0,point0(2,maxres*20),npoint,point(2,maxres*20)
	real occx0(maxres*20),occx(maxres*20)
c
	integer i,ix,j,serix(maxres),ires,jres
c
c	extract iset to work arrays:
c		ca0 -> ca		CA coordinates for chainbreakcheck
c		cs0 -> cs -> cont	cumulative rowsums of contact matrix
c		nres0 -> nres		# residues
c		ltest			residues to test for cut
c
	nres=0
	do i=1,nres0
		serix(i)=0
		if(set(i).eq.iset) then
			nres=nres+1
			resix(nres)=i
			serix(i)=nres
		end if
	end do
c
c	extract iset from ca0,cs0 to ca,cs
c
	do ix=1,nres
		i=resix(ix)
		do j=1,3
			ca(j,ix)=ca0(j,i)
		end do
	end do
c
c	contacts
c
        write(*,*) '#selectset contacts',npoint0
	npoint=0
	do i=1,npoint0
		ires=point0(1,i)
		jres=point0(2,i)
		if((set(ires).eq.iset).and.(set(jres).eq.iset)) then
			npoint=npoint+1
			point(1,npoint)=serix(ires)
			point(2,npoint)=serix(jres)
			occx(npoint)=occx0(i)
		end if
	end do
	if(nres.le.500) call fastcont(nres,occx,npoint,point,cont)
c
c	cumulative number of atoms
c
	do ix=1,nres
		i=resix(ix)
		if(ix.eq.1) n(ix)=natom0(i)
		if(ix.gt.1) n(ix)=n(ix-1)+natom0(i)
	end do
c
c
	call checkchain(ltest, nres, ca, breakdist, minseglen)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine gethelix(nres,struc,nseg,range)
c
c       returns list of H and E segments only
c
        implicit none
        include 'sizes-puu.for'
        integer nres,nseg,range(2,maxres)
        character struc(maxres),c
c
        integer i,k,l,n
c
	do i=1,nres
		c=struc(i)
		if((c.ne.'H').and.(c.ne.'H').and.(c.ne.'I')) struc(i)='L'
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
        do i=1,nseg
                k=range(1,i)
                l=range(2,i)
c                write(*,*) i,k,l,'  ',(struc(j),j=k,l)
        end do

        return
        end
c
c------------------------------------------------------------------------------
c
	subroutine setup_protein(code,chainid,nres,natom0,occx,npoint,
     $  point,
     $	ca,nsearch,search,dssplist,pdblist,bp0,hbond,resno,struc,
     $	protname,lhetatm,pdbfilnam,dsspfile)

	implicit none
	include 'sizes-puu.for'
c
	integer nres,nsearch,search(3,1000),natom0(maxres)
	character*4 code
	character*5 resno(0:maxres)
        character chainid,struc(maxres)
        character*80 dssplist,pdblist,dsspfile
	integer bp0(2,maxres)
	real hbond,ca(3,maxres)
	real occx(maxres*20)
	integer npoint,point(2,maxres*20)
	character*80 protname
	logical lhetatm
c
	character seq(maxres)
        character*33 dsspinfo(maxres)
        integer acc(maxres)
	real xyz(3,maxatm)
	integer natom,atmres(maxatm)
	character*3 resnam(maxres),atmnam(maxatm)
c
	integer i,j,nres1
	character*80 findfile,pdbfilnam,compnd
c
c	setup protein, contact matrix
c
c	write(*,*) ' Reading coordinates from PDB file'
c        pdbfilnam=findfile(code,pdblist,90)
        call readpdb(90,pdbfilnam,chainid,xyz,
     $                  natom,atmnam,atmres,nres,resnam,resno,lhetatm)
c       write(*,*) ' nres: ',nres,' natom: ',natom
c
c	remember C(alpha)s
c
	call grepatom(nres,natom,atmnam,xyz,atmres,ca,'CA ')
c
c	remember #atoms per residue in natom0
c
	do i=1,nres
		natom0(i)=0
	end do
	do i=1,natom
		j=atmres(i)
		natom0(j)=natom0(j)+1
	end do
c
c	calculate contact matrix
c
	call getcontact_shortlist(occx,npoint,point,
     $		nres,natom,atmnam,xyz,atmres,nsearch,search,.true.)
c
c
c
	compnd='unknown'
	do i=1,nres
		do j=1,33
			dsspinfo(i)(j:j)=' '
		end do
		struc(i)=' '
		seq(i)='?'
		acc(i)=0
	end do
	write(*,*) ' Reading DSSP file'
        call readdssp1(dsspfile,chainid,90,
     $                  nres1,seq,struc,acc,dsspinfo,compnd)
	if(nres1.ne.nres) then
		write(*,*) 'WARNING: PDB-DSSP discrepancy',nres,nres1
		do i=1,nres
			bp0(1,i)=0
			bp0(2,i)=0
		end do
	end if
	protname=compnd(11:80)
c
c	sheet info
c
        do i=1,nres1 ! residues read from DSSP
                read(dsspinfo(i)(21:24),*) bp0(1,i)
                read(dsspinfo(i)(25:28),*) bp0(2,i)
c
c               ignore H-bonds to other chain
c
                if(bp0(1,i).eq.-99) bp0(1,i)=0
                if(bp0(2,i).eq.-99) bp0(2,i)=0
c		if(bp0(1,i).gt.0) write(*,*) 'bp0:',i,bp0(1,i),bp0(2,i)
        end do
c	write(*,*) 'sheet info done'
c
c	check and correct consistency of bp0
c
	do i=1,nres1 ! residues read from DSSP
		call checkdssp(i,bp0(1,i),bp0)
		call checkdssp(i,bp0(2,i),bp0)
	end do
c	write(*,*) 'checkdssp done'
c
c	add H-bonds to contact map
c
	do i=1,nres1 ! residues read from DSSP
		j=bp0(1,i)
		if(j.gt.0) call addhbond(i,j,npoint,point,occx,hbond)
c		if(j.gt.0) call addhbond(j,i,npoint,point,occx,hbond)
		j=bp0(2,i)
		if(j.gt.0) call addhbond(i,j,npoint,point,occx,hbond)
c		if(j.gt.0) call addhbond(j,i,npoint,point,occx,hbond)
	end do
c	write(*,*) 'addhbond done'

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine checkdssp(i,j,bp0)
	implicit none
	include 'sizes-puu.for'
	integer i,j,bp0(2,maxres)
c
	if(j.gt.0) then
		if((bp0(1,j).ne.i).and.(bp0(2,j).ne.i)) then
			write(*,*) 'bp0 mismatch in DSSP',i,j
			if(bp0(1,j).eq.0) then
				bp0(1,j)=i
				write(*,*) '... corrected '
			else if(bp0(2,j).eq.0) then
				bp0(2,j)=i
				write(*,*) '... corrected '
			else
				write(*,*) '... unable to correct'
			end if
		end if
	end if

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine addhbond(i,j,npoint,point,occx,hbond)
	implicit none
	include 'sizes-puu.for'
	integer i,j,npoint,point(2,maxres*20)
	real occx(maxres*20),hbond
c
	integer k
c

	if(j.eq.0) return
	k=0
	do while(k.lt.npoint)
		k=k+1
		if((point(1,k).eq.i.and.point(2,k).eq.j)) then
			occx(k)=occx(k)+hbond
c			write(*,*) 'addhbond',i,j,k
			return
		end if
	end do
	npoint=npoint+1
	point(1,npoint)=i
	point(2,npoint)=j
	occx(npoint)=hbond
c	write(*,*) 'new hbond',i,j,npoint

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine checkchain(ltest,nres,ca,distcutoff,minseglen)
c
c	ltest: FALSE within minseglen-1 of chain end
c	chain break = CA-CA distance > distcutoff
c
	implicit none
	include 'sizes-puu.for'
	logical ltest(maxres)
	integer nres,minseglen
	real ca(3,maxres),distcutoff
c
	integer i,j
	logical lbreak(0:maxres)
	real distance
c
c	check for chain breaks
c
	do i=1,maxres
		lbreak(i)=.false.
	end do
	lbreak(0)=.true.
	do i=1,nres-1
		lbreak(i)=(distance(ca(1,i),ca(2,i),ca(3,i),
     $			ca(1,i+1),ca(2,i+1),ca(3,i+1)).gt.distcutoff)
c		if(lbreak(i)) write(*,*) ' chain break detected at ',i
	end do
	lbreak(nres)=.true.
c
c	mark chain ends in ltest
c
	do i=1,maxres
		ltest(i)=(i.le.nres)
	end do
	do i=0,nres
		if(lbreak(i)) then
			do j=i-minseglen+1,i+minseglen-1 ! cut points forward !
				if(j.ge.1.and.j.le.nres) ltest(j)=.false.
			end do
c			chain break residue is allowed at domain border !
			if(i.ge.1) ltest(i)=.true.
		end if
	end do
c	write(*,500) (ltest(i),i=1,nres)
c500	format(80l1)
c
	return
	end
c
c-----------------------------------------------------------------------
c
        subroutine readdssp1(filnam,chainid,iunit,nres,seq,struc,
     $          acc,dsspinfo,compnd)
c
c       frontend to getdssp to transport wanted subset of DSSP information
c
        implicit none
	include 'sizes-puu.for'
        integer iunit,nres,nchain,nss,nssintra,nssinter,nhb,acc(maxres)
        character*80 filnam,header,compnd,source,author
        character chainid,seq(maxres),struc(maxres)
        real nhbp,tco(maxres),kappa(maxres),alfa(maxres),phi(maxres),
     $          psi(maxres),ca(3,maxres)
        character*5 resno(0:maxres)
        character*33 dsspinfo(maxres)
        integer area,ierr,dsspresno(0:maxres),fsspresno(0:maxres)
c
        call getdssp(filnam,chainid,iunit,maxres,header,compnd,source,
     $          author,nres,nchain,nss,nssintra,nssinter,area,
     $          resno,seq,struc,acc,dsspinfo,dsspresno,fsspresno,
     $          tco,kappa,alfa,phi,psi,ca,ierr)
        write(*,*) header,compnd,source,author
c        write(*,*) ierr
c        write(*,*) nres
        if(ierr.eq.0) write(*,*) nres,' residues read '
        if(ierr.eq.-1) write(*,*) ' DSSP file not found '
        if(ierr.eq.-2) write(*,*) ' DSSP file truncated at ',maxres,
     $          'residues'

        return
        end
c
c------------------------------------------------------------------------------
c
        function file_exists(infile,iunit)
        implicit none
        logical file_exists
        character*80 infile
        character machine
        integer iunit
c
        file_exists=.false.
        machine='U'
        if(machine.eq.'V') then
                open(iunit,file=infile,err=19,status='old')
        else
                open(iunit,file=infile,err=19,status='old')
        end if
        file_exists=.true.
19      close(iunit)

        return
        end
c
c------------------------------------------------------------------------------
	subroutine writeparams(iunit,breakdist,minseglen,hbond,compact,
     $		taucutoff,seed,printdepth,lsplit,lhetatm,lseed)
	implicit none
	integer iunit,minseglen,seed,printdepth
	real hbond,compact,taucutoff,breakdist
	logical lsplit,lhetatm,lseed
c
	write(iunit,500) hbond,taucutoff,compact
	if(lsplit) then
		write(iunit,510) minseglen,breakdist,'YES'
	else
		write(iunit,510) minseglen,breakdist,'NO '
	end if
	if(lhetatm) then
		if(.not.lseed) then
			write(iunit,520) 'DSSP+HETATM',printdepth
		else
			write(iunit,530) 'DSSP+HETATM',printdepth,seed
		end if
	else
		if(.not.lseed) then
			write(iunit,520) 'DSSPclassic',printdepth
		else
			write(iunit,530) 'DSSPclassic',printdepth,seed
		end if
	end if
c
500	format('PARAMETER hbond =',f5.1,'; taucutoff =',f5.1,
     $ '; compact =',f5.1)
510	format('PARAMETER minseglen =',i5,'; breakdist =',f5.1,
     $ '; split sheets: ',a3)
520	format('PARAMETER pdbreader = ',a11,'; depth =',i5)
530	format('PARAMETER pdbreader = ',a11,'; depth =',i5,
     $  '; seed = ',i20)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine writedomainnotation(iunit)
	implicit none
	integer iunit
c
	write(iunit,500) ' '
	write(iunit,510) '//'
c
500	format('NOTATION  PDBID: PDB code of coordinate set.',
     $ /'NOTATION  COMPND: echoed from the PDB file.',
     $ /'NOTATION  CHAIN: chain identifier (* is wildcard).',
     $ /'NOTATION  NDOM: number of domains.  ',
     $  'Units > 40 residues are putative autonomous',
     $ /'NOTATION   units (structural domains) and numbered ',
     $ 'sequentially.  Small units',
     $ /'NOTATION   are labelled alphabetically.',
     $ /'NOTATION  Secondary structural class: based on the number of',
     $ ' residues in helix',
     $ /'NOTATION   and in strand (taken from DSSP, Kabsch & Sander, ',
     $ 'Biopolymers ',
     $ /'NOTATION   22:2577-2637, 1983):',
     $ /'NOTATION   all-alpha:   > 40 % helix and < 15 % strand',
     $ /'NOTATION   all-beta:    < 15 % helix and > 30 % strand',
     $ /'NOTATION   alpha+beta: N-term half is all-alpha and C-term',
     $ /'NOTATION                  half is all-beta, or vice versa',
     $  /'NOTATION   alpha/beta:  > 15 % helix and > 15 % strand',a1)
510	format('NOTATION   mixed:       otherwise',
     $ /'NOTATION  Shape: two classes (compact/nonglobular) based on',
     $ ' intraunit',
     $ /'NOTATION   tertiary (long-sequence-range) ',
     $ 'interatomic contact density.',
     $ /'NOTATION  Seqtl resnos: sequential numbering.',
     $ /'NOTATION  PDB resnos: residue ranges echoed from the ',
     $ 'PDB file. ',
     $ /'NOTATION   Residue ranges of discontinuous domains are given ',
     $ 'on adjacent lines.',
     $ /'NOTATION  Protein entries are separated by double ',
     $ 'slashes (//).',
     $ /a2)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine writeunitsnotation(iunit)
	implicit none
	integer iunit
c
	write(iunit,500) ' '
	write(iunit,510) '//'
c
500	format('NOTATION  PDBID: PDB code of coordinate set.',
     $ /'NOTATION  CHAIN: chain identifier (* is wildcard).',
     $ /'NOTATION  NUNITS:total number of listed unfolding units.',
     $ /'NOTATION  Children: point to units in the unfolding tree',
     $ /'NOTATION  Types of unfolding unit:',
     $ /'NOTATION   * structural domain',
     $ /'NOTATION   = short, flexible piece of chain (<40 residues)',
     $ /'NOTATION   + above structural domains in unfolding tree',
     $ /'NOTATION   - below structural domains in unfolding tree',
     $ /'NOTATION  Tau-squared: interdomain fluctuation time [ps**2]',
     $ /'NOTATION  Compactness: intraunit tertiary (long-sequence-',
     $ 'range) interatomic contact density',a1)
510	format('NOTATION  Hbonds: number of H-bonds across sub-domain ',
     $ 'interface',
     $ /'NOTATION  Size: number of residues in unit',
     $ /'NOTATION  Segments: number of chain segments in unit',
     $ /'NOTATION  Residues: sequentially numbered (PDB numbers in',
     $ ' parentheses)',
     $ /'NOTATION  Protein entries are separated by double ',
     $ 'slashes (//).',
     $ /a2)
c
	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine writeheadertobase(iunit,headerstring,versionstring,
     $		basestring)
	implicit none
	integer iunit
	character*(*) headerstring,versionstring,basestring
c
	integer l
	character cdate*24
c
	l=len(headerstring)
	write(iunit,500) headerstring
	call getdate(cdate)
        write(iunit,510) cdate
	l=len(versionstring)
	write(iunit,520) versionstring
	l=len(basestring)
	write(iunit,530) basestring
c
500	format('PUU       ',a240)
510	format('DATE      File generated on ',a24)
520	format('VERSION   ',a240)
530	format(
     $ 'REFERENCE L. Holm, C. Sander, Proteins 19:256-268, 1994.',
     $ /'AVAILABLE Free academic use. No commercial use. No inclusion ',
     $ /'AVAILABLE  in other databases without permission.',
     $ /'CONTACT   Holm@EMBL-Heidelberg.DE, Sander@EMBL-Heidelberg.DE',
     $ /'CONTACT    / phone +49-6221-387361 / fax +49-6221-387306',
     $ /'BASE      ',a240)
c
	return
	end
c
c-----------------------------------------------------------------------
c
c
CCCCCCC grid.f
c
c------------------------------------------------------------------------------
c
	function lgrid(gx,gy,gz,maxgrid)
c
c	check that g* are inside -maxgrid..+maxgrid
c	.true. if ok, .false. if outside
c
	implicit none
	logical lgrid
	integer gx,gy,gz,maxgrid
c
	lgrid=.true.
	if(gx.gt.maxgrid) lgrid=.false.
	if(gx.lt.-maxgrid)lgrid=.false.
	if(gy.gt.maxgrid) lgrid=.false.
	if(gy.lt.-maxgrid)lgrid=.false.
	if(gz.gt.maxgrid) lgrid=.false.
	if(gz.lt.-maxgrid)lgrid=.false.

	return
	end
c
c------------------------------------------------------------------------------
c
	function fung(x)
	implicit none
	integer fung
	real x

	fung=nint(x/2)

	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine init_search(search,nsearch,realcutoff)
	implicit none
	integer nsearch,fung,search(3,1000)
	real realcutoff
c
	integer dx,dy,dz,gx,gy,gz,cutoff,cutoff2,d2
	logical test(-20:20,-20:20,-20:20)
	integer minx,miny,minz,maxx,maxy,maxz
	integer x2,y2,z2
c
	minx=0
	miny=0
	minz=0
	maxx=0
	maxy=0
	maxz=0
	nsearch=0
	cutoff=nint(realcutoff+0.5)
	cutoff2=nint(realcutoff*realcutoff)
c	write(*,*) cutoff2
	do dx=-(cutoff+2),cutoff+2
		do dy=-(cutoff+2),cutoff+2
			do dz=-(cutoff+2),cutoff+2
				gx=fung(float(dx))
				gy=fung(float(dy))
				gz=fung(float(dz))
				test(gx,gy,gz)=.false.
			end do
		end do
	end do
	do dx=-(cutoff+2),cutoff+2
		do dy=-(cutoff+2),cutoff+2
			do dz=-(cutoff+2),cutoff+2
				gx=fung(float(dx))
				gy=fung(float(dy))
				gz=fung(float(dz))
				if(gx.lt.minx) minx=gx
				if(gy.lt.miny) miny=gy
				if(gz.lt.minz) minz=gz
				if(gx.gt.maxx) maxx=gx
				if(gy.gt.maxy) maxy=gy
				if(gz.gt.maxz) maxz=gz
				gx=gx*2
				gy=gy*2
				gz=gz*2
				x2=min((gx-2)*(gx-2),(gx+2)*(gx+2),gx*gx)
				y2=min((gy-2)*(gy-2),(gy+2)*(gy+2),gy*gy)
				z2=min((gz-2)*(gz-2),(gz+2)*(gz+2),gz*gz)
				d2=x2+y2+z2
				if(d2.lt.cutoff2) test(gx/2,gy/2,gz/2)=.true.
			end do
		end do
	end do
	do dx=minx,maxx
		do dy=miny,maxy
			do dz=minz,maxz
				if(test(dx,dy,dz)) then
					nsearch=nsearch+1
					search(1,nsearch)=dx
					search(2,nsearch)=dy
					search(3,nsearch)=dz
				end if
			end do
		end do
	end do
c	write(*,510) nsearch,realcutoff

500	format(5i10,f10.2)
510	format(i10,' grid cells to be searched per atom at cutoff radius',
     $    f8.2)

	return
	end
c
c------------------------------------------------------------------------------
c
	subroutine getcontact_shortlist(occx,npoint,point,
     $		nres,natom,atmnam,xyz,atmres,nsearch,search,lbackbone)
c
c	writes shortlisted contacts in occx, point
c
c	lbackbone: .true. -> include backbone atoms; .false. -> exclude backbone
c
	implicit none
	include 'sizes-puu.for'
	real contact_strength,cont(maxres)
	real xyz(3, maxatm)
	integer natom, atmres(maxatm), nres
	character*3 atmnam(maxatm)
	logical lbackbone
	real occx(maxres*20)
	integer npoint,point(2,maxres*20)
c
	real waterscale
	parameter(waterscale=0.31)
c
	integer grid(-maxgrid:maxgrid,-maxgrid:maxgrid,-maxgrid:maxgrid)
	integer overlay(maxatm),nocc,occ(3,maxatm)
	logical activeatm(maxatm)
c
	real distance
	integer gx,gy,gz,nsearch,search(3,1000),fung
	logical lgrid,lskip(maxres)
c
	integer i,j,ires,jres,icell,gx1,gy1,gz1
c
	integer oldires
	logical lset2
	integer firstatom(maxres),lastatom(maxres)
	real x
c
c	write(*,*) ' Calculating contact matrix'
c
c	center molecule !
c
	call center(xyz,natom)
c
c	initialize skiplist: set1=in box; set2=outside
c
	do ires=1,nres
		lskip(ires)=.false.
	end do
	do i=1,natom
		ires=atmres(i)
		if(.not.lskip(ires)) then
			gx=fung(xyz(1,i))
			gy=fung(xyz(2,i))
			gz=fung(xyz(3,i))
			if(.not.lgrid(gx,gy,gz,maxgrid)) lskip(ires)=.true.
		end if
	end do
c
c	set flag for non-empty set2
c
	i=0
	do ires=1,nres
		if(lskip(ires)) i=i+1
	end do
	lset2=(i.gt.0)
c	write(*,*) i,' residues are outside grid '
c
c	put set 1 in box
c
	do i=1,natom
		ires=atmres(i)
		activeatm(i)=(.not.lskip(ires))
		if(.not.lbackbone.and.activeatm(i)) then
			if(atmnam(i).eq.'C  ') activeatm(i)=.false.
			if(atmnam(i).eq.'N  ') activeatm(i)=.false.
			if(atmnam(i).eq.'O  ') activeatm(i)=.false.
			if(atmnam(i).eq.'OXT') activeatm(i)=.false.
		end if
	end do
c
	call fillgrid(grid,overlay,xyz,natom,nocc,occ,activeatm)
c
c	calculate contact matrix
c
c	get pairs: set1 * set1 + set 1 * set 2
c
	oldires=0
	npoint=0
	do i=1,natom+1
c	  special case: save last residue
	  if(i.gt.natom) then
		if(oldires.gt.0) then
		  if(.not.lskip(oldires)) call
     $		    savecontacts(oldires,cont,nres,npoint,point,occx)
		end if
	  else if(activeatm(i)) then
		ires=atmres(i)
		if(ires.ne.oldires) then
c
c			save contacts of previous residue
c
			if(oldires.gt.0) then
			  if(.not.lskip(oldires)) call
     $			    savecontacts(oldires,cont,nres,npoint,point,occx)
			end if
c
c			initialize contacts of current residue
c
			do jres=1,nres
				cont(jres)=0.0
			end do
			oldires=ires
		end if
		gx=fung(xyz(1,i))
		gy=fung(xyz(2,i))
		gz=fung(xyz(3,i))
		do icell=1,nsearch
			gx1=gx+search(1,icell)
			gy1=gy+search(2,icell)
			gz1=gz+search(3,icell)
			if(.not.lgrid(gx1,gy1,gz1,maxgrid)) goto 19
			j=grid(gx1,gy1,gz1)
			do while(j.gt.0)
				jres=atmres(j)
c
c				exclude contacts from ires to ires,...,ires+-3
c
				if(abs(jres-ires).ge.0)
     $		cont(jres)=cont(jres)+contact_strength(
     $		distance(xyz(1,i),xyz(2,i),xyz(3,i),
     $			 xyz(1,j),xyz(2,j),xyz(3,j)))
				j=overlay(j)
			end do
19		end do
c
c		add contacts with set 2
c
		if(lset2) then
			if(.not.lbackbone) stop ' dunno how to exclude backbone'
			do j=1,natom
				jres=atmres(j)
				if(lskip(jres)) then
					if(abs(jres-ires).ge.0)
     $		cont(jres)=cont(jres)+contact_strength(
     $		distance(xyz(1,i),xyz(2,i),xyz(3,i),
     $			 xyz(1,j),xyz(2,j),xyz(3,j)))
				end if
			end do
		end if
	  end if
	end do
c
c	do contacts of residues in skiplist: set 2 * (set 1+ set 2)
c
	if(lset2) then
	  do ires=1,nres
		firstatom(ires)=0
		lastatom(ires)=0
	  end do
	  do i=1,natom
		ires=atmres(i)
		if(firstatom(ires).eq.0) firstatom(ires)=i
		lastatom(ires)=i
	  end do
	  do ires=1,nres
	    if(lskip(ires)) then
		write(*,*) ' calculate residue outside grid:',ires
		do jres=1,nres
			cont(jres)=0.0
		end do
		do jres=1,nres
c		  if(abs(jres-ires).ge.0) then
			i=firstatom(ires)
			j=firstatom(jres)
			x=distance(xyz(1,i),xyz(2,i),xyz(3,i),
     $				   xyz(1,j),xyz(2,j),xyz(3,j))
			if(x.lt.25.0) then
				do i=firstatom(ires),lastatom(ires)
				  do j=firstatom(jres),lastatom(jres)
				    cont(jres)=cont(jres)+contact_strength(
     $		distance(xyz(1,i),xyz(2,i),xyz(3,i),
     $			 xyz(1,j),xyz(2,j),xyz(3,j)))
				  end do
				end do
			end if
c		  end if
		end do
		call savecontacts(ires,cont,nres,npoint,point,occx)
	    end if
	  end do
	end if

	return
	end
c
c-----------------------------------------------------------------------------
c
	subroutine  savecontacts(ires,cont,nres,npoint,point,occx)
	implicit none
	include 'sizes-puu.for'
	integer ires,nres,npoint,point(2,maxres*20)
	real occx(maxres*20),cont(maxres)
c
	integer jres
	real cs
c
	do jres=1,nres
		cs=cont(jres)
		if(cs.gt.0.0) then
			npoint=npoint+1
			if(npoint.gt.maxres*20) stop 'npoint overflow'
			point(1,npoint)=ires
			point(2,npoint)=jres
			occx(npoint)=cs
		end if
	end do
c
	return
	end
c
c-----------------------------------------------------------------------------
c
	function contact_strength(r)
	real r, contact_strength
c
	real rmin,tworw,rhi
c!!!	count atomic contacts at less than 4 Angstrom
	parameter(rmin=4.0,tworw=0.0,rhi=rmin+tworw)
c
	if(r.le.rmin)then
		contact_strength=1.0
c	else if(r.lt.rhi) then
c		contact_strength=1.0-(r-rmin)/tworw
	else
		contact_strength=0.0
	end if

	return
	end
c
c-----------------------------------------------------------------------------
c
	subroutine fillgrid(grid, overlay, xyz, natom,
     $			  nocc, occ, activeatm)
c
c	center molecule at (0,0,0) by shifting coordinates x(),y(),z()
c	put atom numbers to respective grid cells or overlay array
c	clear grid of old atoms
c	insert only active atoms (activeatm(iatm)=.true.)
c
  	implicit none
   	include 'sizes-puu.for'
    	integer natom
      	real xyz(3, maxatm)
      	integer grid(-maxgrid:maxgrid, -maxgrid:maxgrid,
     $			   -maxgrid:maxgrid)
        integer overlay(maxatm), nocc, occ(3, maxatm)
	logical activeatm(maxatm)
c
	logical lgrid
c
 	integer fung, gx,gy,gz
        integer hist(0:maxatm), i,j,m,n,k,nx
c
c	initialize grid locally
c
	nocc=0
	do i=-maxgrid,maxgrid
		do j=-maxgrid,maxgrid
			do k=-maxgrid,maxgrid
				grid(i,j,k)=0
			end do
		end do
	end do
c
	call center(xyz,natom)
c
c	zero grid and overlay -- all nocc atoms in occ()
c
	do i=1,nocc
		grid(occ(1,i), occ(2,i), occ(3,i))=0
	end do
	hist(0)=0
	do i=1,natom
		overlay(i)=0
		hist(i)=0
	end do
c
c	put atoms in grid -- only actives
c	(or pointers to overlay)
c	!!! warn if atoms are outside grid !!!
c
	nx=0
	m=0
	nocc=0
	do i=1,natom
	  if(activeatm(i)) then
		gx=fung(xyz(1,i))
		gy=fung(xyz(2,i))
		gz=fung(xyz(3,i))
c		check inside
		if (.not.lgrid(gx,gy,gz,maxgrid)) then
			nx=nx+1
c			write(*,*) 'outside:',nx,i,(xyz(j,i),j=1,3),gx,gy,gz
			goto 19
		end if
		j=grid(gx,gy,gz)
c		write(*,*) i,j,gx,gy,gz
		if (j.eq.0) then
			grid(gx,gy,gz)=i
			nocc=nocc+1
			occ(1,nocc)=gx
			occ(2,nocc)=gy
			occ(3,nocc)=gz
			hist(0)=hist(0)+1
		else
			n=0
10			k=overlay(j)
			n=n+1
			if (n.gt.m) m=n
			if (k.eq.0) then
				overlay(j)=i
				hist(n)=hist(n)+1
			else
				j=k
				go to 10
			end if
		end if
19	  end if
	end do

c	write(iout,40) m, nocc, natom, nx
40 	format (/' maximum overlap of atoms in one grid cell: ', i6,
     $        /' number of occupied grid cells: ',i6,
     $        /' natom: ',i6, /' atoms outside grid: ', i6)
c	s=0
c	do i=0,m
c		hist(i)=hist(i)-hist(i+1)
c		s=s+hist(i)*(i+1)
c	end do
c	write(iout,*) ' sum: ',s
c	histogram of cell occupancy
c	write(iout,42) (hist(i),i=0,m)
42	format(/' cell density histogram: ',/10(i8))

      return
      end
c
c-----------------------------------------------------------------------------
c
	subroutine center(xyz,natom)
  	implicit none
   	include 'sizes-puu.for'
    	integer natom
      	real xyz(3, maxatm)
c
        real xyz0(3)
	integer i,j
c
c	center molecule at 0,0,0
c
	do j=1,3
		xyz0(j)=0
	end do
c	write(*,*) natom
	if (natom.eq.0) then
		write(*,*) 'GRID -- No atoms !'
		return
	end if
	do i=1,natom
	  do j=1,3
		xyz0(j)=xyz0(j)+xyz(j,i)
	  end do
	end do
	do j=1,3
		xyz0(j)=xyz0(j)/natom
	end do
c
c	shift all atoms by -x0,-y0,-z0
c
	do i=1,natom
	  do j=1,3
		xyz(j,i)=xyz(j,i)-xyz0(j)
	  end do
	end do

      return
      end
c
c-----------------------------------------------------------------------------
c
CCCCCCCC excerpt from generic.f
c
c------------------------------------------------------------------------------
c
	function findfile(code,listfile,iunit)
	implicit none
	character*80 listfile,findfile
	character*4 code
	integer iunit
c
	integer i
	character*4 c
	character*80 f
	logical lfound
c
c	generic function to scan listfile for code and read complete filename
c	format of listfile is a4,5x,a80
c	use e.g. ~/lib/dssp.list or ~/lib/pdb.list
c
	if(listfile.eq.'NONE') then
		write(*,*) 'enter file name'
		read(*,540) findfile
		return
	end if
c
	lfound=.false.
	do i=1,80
		findfile(i:i)=' '
	end do
	write(*,500) code, listfile
	open(iunit,file=listfile,status='old')
	i=0
10	read(iunit,510,end=19) c,f
	i=i+1
	if(c.eq.code) then
		findfile(1:60)=f
		write(*,520) i,f
		lfound=.true.
		goto 19
	end if
	goto 10
19	if(.not.lfound) then
		write(*,530) i
		findfile='?'
	end if
	close(iunit)

500	format(/'Looking for ',a4,' in ',a60)
510	format(a4,5x,a80)
520	format('position:',i4,2x,a60)
530	format('No match found in ',i4,' entries')
540	format(a80)

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine readpdb(iunit,filnam,chainid,
     $		xyz,natom,atmnam,atmres,nres,resnam,resno,lhetatm)
c
c	frontend to getcoor COMPATIBLE WITH DSSP
c	(classical version lhetatm=.false.)
c	>>> uses 'sizes-puu.for' <<<
c	lhetatm: .true. compatible m. scharf's modification to dssp at embl !
c
c	extracts residues with chainid (* is wildcard)
c
c	parameter:
c		maxres
c		maxatm
c
c	input:	filnam
c		chainid
c		iunit
c
c	output:	header, compnd, source, author -- echo pdb cards
c		nres,natom
c		resno -- pdb residue number string with appended insertion
c			 indicator, e.g. 1acx 63a
c		resnam -- 'ala'
c		atmnam -- 'ca '
c		atmres -- residue index of atom
c		bvalue -- real
c		xyz -- atomic coordinates
c		ierr = 0 -- normal
c		     = -1 -- file not found
c		     = -2 -- too many residues, chain truncated to maxres
c
c	acceptable residue: onelettercode.ne.'-'
c	new residue: different resno from previous one
c	ignore atoms which start 'h' or 'd'
c	ignore residues with incomplete backbone (n,ca,c,o)
c
	implicit none
	include 'sizes-puu.for'
	integer iunit
	character*80 filnam,header,compnd,source,author
	character chainid
	real xyz(3,maxatm),bvalue(maxatm)
	integer natom,atmres(maxatm),nres,ierr
	character*3 atmnam(maxatm),resnam(maxres)
	character*5 resno(0:maxres)
	logical lhetatm
c
c	write(*,*) 'reading ',filnam(1:60)
	call getcoor(iunit,maxres,maxatm,filnam,chainid,header,compnd,
     $		source,author,xyz,natom,atmnam,atmres,nres,
     $		resno,resnam,ierr,lhetatm)
c	write(*,*) ierr
	write(*,*) header
	write(*,*) compnd
	write(*,*) source
	write(*,*) author

	return
	end
c
c-----------------------------------------------------------------------
c
	function distance(x1,y1,z1,x2,y2,z2)
	implicit none
	real distance,x1,y1,z1,x2,y2,z2

	distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))

	return
	end
c
c------------------------------------------------------------------------------
c

	subroutine correspondence_shortlist(individuals,attributes,
     $		occx,point,npoint,x,r2,
     $		seed0,eigenvalues)
c
c	sparse matrix version, shortlisted contacts given as input
c
	implicit none
	include 'sizes-puu.for'
	integer eigenvalues
	integer individuals, attributes,npoint,point(2,maxres*20)
	real x(3,maxres),occx(maxres*20)
c
c	input data:
c		inames = names of individuals
c		individuals = their #
c		anames = names of attributes
c		attributes = their #
c		occ(attrib,individ) = incidence counts
c		eigenvalues = 3
c	local:
c		r2 = eigenvalues
c		x,y = eigenvaluevectors
c		rsums,csums
c
	real r2(3),xc(maxres)
c!!!
	real y(3,maxres)
	real rsums(maxres),csums(maxres)
c
	integer i,j,index,round
	real prev_norm,norm
	integer*4 seed0,seed
	real ran
c
	real xx,f(maxres)
	integer ix
c
123	format(20i4)
c
c	initialize all variables
c
	seed=seed0
	do i=1,eigenvalues
		r2(i)=0.0
		do j=1,attributes
			x(i,j)=0.0
		end do
		do j=1,individuals
			y(i,j)=0.0
		end do
	end do
	do j=1,attributes
		rsums(j)=0.0
	end do
	do j=1,individuals
		csums(j)=0.0
	end do
	prev_norm=0.0
c
c
c	write(*,*) '*Computation of eigenvalues begins',individuals,attributes
c	write(*,*) 'sparse matrix: ',npoint,' points',individuals,' residues'
c
c	use shortlisted contacts
c
	do ix=1,npoint
		i=point(1,ix)
		j=point(2,ix)
		xx=occx(ix)
		rsums(j)=rsums(j)+xx
		csums(i)=csums(i)+xx
	end do
c
c	check for zero csum or rsum
c
	do i=1,individuals
		if(csums(i).eq.0) write(*,*) ' FATAL: individual ',
     $			i,' has zero incidence !'
	if(csums(i).ne.rsums(i)) write(*,*) 'mismatch',i,csums(i),rsums(i)
	end do
	do i=1,attributes
		if(rsums(i).eq.0) write(*,*) ' FATAL: attribute ',
     $			i,' has zero incidence !'
	end do
c
	do index=1,eigenvalues
c		define trivial solution of random trial for the components
c			of x:
		do i=1,attributes
			if(index.eq.1) then
				x(1,i)=0.5
			else
				x(index,i)=ran(seed)
			end if
		end do
		call schmidt(index,attributes,x,xc,norm)
		round=0
c		Test for convergence criterion:
		do while((abs(norm-prev_norm).gt.0.00001*norm).and.
     $			(round.lt.100))
			round=round+1
c			if(mod(round,100).eq.0) write(*,*) ' round:',round,
c     $	'; residual:',norm-prev_norm

			do j=1,individuals
				f(j)=0.0
			end do
			do ix=1,npoint
				j=point(2,ix)
				i=point(1,ix)
				f(j)=f(j)+x(index,i)*occx(ix)
			end do
			do j=1,individuals
				y(index,j)=f(j)/csums(j)
			end do
			do i=1,attributes
				f(i)=0.0
			end do
			do ix=1,npoint
				j=point(2,ix)
				i=point(1,ix)
				f(i)=f(i)+y(index,j)*occx(ix)
			end do
			do i=1,attributes
				x(index,i)=f(i)/rsums(i)
			end do
			prev_norm=norm
			call schmidt(index,attributes,x,xc,norm)
		end do
c		write(*,500) index,norm,round,individuals
		r2(index)=norm
	end do
c	write(*,*) '*Computation of eigenvalues ends'
c
500	format(i10,'. eigenvalue equal to ',f7.4,' computed in',i5,
     $ ' iterations',i5)
c	return
	end
c
c-----------------------------------------------------------------------
c
c
c--------------------------------------------------------------------------
c
	subroutine schmidt(index,attributes,x,xc,norm)
c
c	Schmidt orthogonalization and normalization of attribute scores
c
	implicit none
	integer index,attributes
	real x(3,attributes),xc(attributes),norm
c
	integer i,j
	real ssq,e,dotp(100)
c
	do j=1,100
		dotp(j)=0.0
	end do
	do i=1,attributes
		do j=1,index-1
			dotp(j)=dotp(j)+x(index,i)*x(j,i)
		end do
	end do
	ssq=1.0e-9
	do i=1,attributes
		e=x(index,i)
		do j=1,index-1
			e=e-dotp(j)*x(j,i)
		end do
		xc(i)=e
		ssq=ssq+e*e
	end do
	norm=sqrt(ssq)
	do i=1,attributes
		x(index,i)=xc(i)/norm
	end do

	return
	end
c
c-----------------------------------------------------------------------
c
	subroutine defaults(dssplist,pdblist,breakdist,minseglen,hbond,
     $		compact,taucutoff,seed,treefilename,prettyfilename,
     $		defaultsfilename,lhetatm,lsplit,printdepth,basedefault,
     $		subunitfilename)
	implicit none
	include 'sizes-puu.for'
	character*80 defaultsfilename,subunitfilename
c
	integer*4 seed
        character*80 dssplist,pdblist,treefilename,prettyfilename
	character*60 basedefault
	integer minseglen,printdepth
	real breakdist,hbond,compact,taucutoff
	logical lhetatm,lsplit
c
	integer nparam
	parameter(nparam=15)
	character*20 key(nparam)
	data key/'dssplist','pdblist','breakdist','minseglen','hbond',
     $		'compact','taucutoff','seed','treefilename',
     $		'prettyfilename','lhetatm','lsplit','printdepth',
     $		'base','subunitfilename'/
	character c
	character*20 word
c
	integer i,k
c
c	hardwired defaults
c
	dssplist='NONE'
	pdblist='NONE'
	breakdist=5.0
	minseglen=10
	hbond=15.0
	compact=80.0
	taucutoff=2.6
	treefilename='units.puu'
	prettyfilename='domains.puu'
	subunitfilename='subunits.puu'
	seed=1234567811
c       dsspCMBI reads HETATMs (e.g. MSE) as amino acids
	lhetatm=.true.
	lsplit=.false.
	printdepth=-99
	basedefault='unknown'
c
c	check for puu.default
c
	write(*,*) ' read defaults from ',defaultsfilename
	open(90,file=defaultsfilename,status='old',err=99)
10	read(90,500,end=19) c,word
	if(c.eq.'#') then
		k=0
		i=0
		do while(i.lt.nparam)
			i=i+1
			if(word.eq.key(i)) then
				k=i
	if(i.eq.1) read(90,510) dssplist
	if(i.eq.2) read(90,510) pdblist
	if(i.eq.3) read(90,*) breakdist
	if(i.eq.4) read(90,*) minseglen
	if(i.eq.5) read(90,*) hbond
	if(i.eq.6) read(90,*) compact
	if(i.eq.7) read(90,*) taucutoff
	if(i.eq.8) read(90,*) seed
	if(i.eq.9) read(90,510) treefilename
	if(i.eq.10) read(90,510) prettyfilename
	if(i.eq.11) read(90,*) lhetatm
	if(i.eq.12) read(90,*) lsplit
	if(i.eq.13) read(90,*) printdepth
	if(i.eq.14) read(90,520) basedefault
	if(i.eq.15) read(90,510) subunitfilename
			end if
		end do
		if(k.gt.0)write(*,*) key(k),' reset from puu.default'
	end if
	goto 10
c
c	'error' exit
c
99	write(*,*) ' default file not found: that''s ok '
	goto 20
c
c	'normal' exit
c
19	close(90)
c	write(*,*) ' closing file'
	goto 20
c
c	echo parameters
c
20	write(*,*) ' current parameters:'
	do i=1,nparam
		if(i.eq.1) write(*,*) key(i),dssplist(1:40)
		if(i.eq.2) write(*,*) key(i),pdblist(1:40)
		if(i.eq.3) write(*,*) key(i),breakdist
		if(i.eq.4) write(*,*) key(i),minseglen
		if(i.eq.5) write(*,*) key(i),hbond
		if(i.eq.6) write(*,*) key(i),compact
		if(i.eq.7) write(*,*) key(i),taucutoff
		if(i.eq.8) write(*,*) key(i),seed
		if(i.eq.9) write(*,*) key(i),treefilename(1:40)
		if(i.eq.10) write(*,*) key(i),prettyfilename(1:40)
		if(i.eq.11) write(*,*) key(i),lhetatm
		if(i.eq.12) write(*,*) key(i),lsplit
		if(i.eq.13) write(*,*) key(i),printdepth
c		if(i.eq.14) write(*,*) key(i),basedefault
		if(i.eq.15) write(*,*) key(i),subunitfilename(1:40)
	end do
	write(*,*) ' '
c
500	format(a1,a20)
510	format(a80)
520	format(a60)
c
	return
	end
c
c--------------------------------------------------------------------------
c
        subroutine grepatom(nres,natom,atmnam,xyz,atmres,cb,atomname)
c
c       scan xyz-atmnam, return only atoms with 'atomname' in cb
c
        implicit none
        include 'sizes-puu.for'
        integer nres,natom,atmres(natom)
        character*3 atmnam(natom),atomname
        real xyz(3,natom),cb(3,nres)
c
        integer ires, iatm, i
        logical cbfound(maxres)
        character*3 atnm
c
        do ires=1,nres
                cbfound(ires)=.false.
        end do
        do iatm=1, natom
                ires=atmres(iatm)
                atnm=atmnam(iatm)
                if(atnm.eq.atomname) then
                        cbfound(ires)=.true.
                        do i=1,3
                                cb(i,ires)=xyz(i,iatm)
                        end do
                end if
        end do
        do ires=1, nres
                if(.not.cbfound(ires)) write(*, 150) ires, atomname
        end do
150     format(' Warning: cannot grep ',i5,2x,a3,' missing')

        return
        end
c
c-----------------------------------------------------------------------
c
	subroutine pret2(ndom,node_tau,node_child,node_range,node_nseg,
     $		node_parent,node_size,code,chainid,
     $		node_tert,node_nbp,resno,struc,compact,
     $		protname,printdepth,node_type)
c
c	write out compact domains
c
	implicit none
	include 'sizes-puu.for'
	integer ndom,node_child(maxres,2),node_range(maxres,2,maxnseg)
	integer node_nseg(maxres),node_size(maxres),printdepth
	real node_tau(maxres),compact
	character chainid
	character*4 code
	character*80 protname
c
	integer i,j,k,m,n
	integer node_parent(maxres)
	character*32 node_label(maxres)
	real node_tert(maxres)
	integer node_nbp(maxres)
	integer na,nb,na1,na2,nb1,nb2,ix
	character node_type(maxres),tab,struc(maxres)
	character*2 arm(0:25),dom(0:25),cd
	data arm/' ',' a',' b',' c',' d',' e',' f',' g',' h',' i',' j',
     $		     ' l',' m',' n',' o',' p',' q',' r',' s',' t',' u',
     $		     ' v',' w',' x',' y',' z'/
	data dom/'  ',' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',
     $		 '11','12','13','14','15','16','17','18','19','20',
     $		 '21','22','23','24','25'/
	character*5 resno(0:maxres)
	character*10 class
	character*11 shape
	integer nstrucdom,nseg,iseg,jseg,a1,a2,d(maxres),v(maxres)
	character c1,c2
c
	tab=char(9)
c
	call labels(code,chainid,ndom,node_child,node_label)
c
c	set ndom to zero if no residues
c
	if(node_size(1).eq.0) then
	 	write(*,*) 'WARNING: empty root ',code,chainid
		ndom=0
	end if
c
c	15,32: header !
c
	write(32,700) code
	write(15,700) code
	write(32,710) protname
	if(chainid.ne.' ') then
		write(32,720) chainid
		write(15,720) chainid
	end if
	nstrucdom=0
	do j=1,ndom
	   if(node_type(j).eq.'*'.or.node_type(j).eq.'=') 
     $          nstrucdom=nstrucdom+1
	end do

	m=0
	n=0
	ix=1
	do j=1,ndom
	  if(node_type(j).eq.'*'.or.node_type(j).eq.'=') then
		m=m+1
c!		if(node_size(1).gt.80) then
			if(node_size(j).le.40) then
				n=n+1
				node_type(j)='='
			end if
			ix=ix+1
c
c			32: first line
c
c			split range segments with multiple chains !
c
			if(chainid.eq.'*') then
			  nseg=node_nseg(j)
			  do i=1,nseg
			    c1=resno(node_range(j,1,i))(1:1)
			    c2=resno(node_range(j,2,i))(1:1)
			    if(c1.ne.c2) then
				iseg=i
				a1=node_range(j,1,iseg)
				a2=node_range(j,2,iseg)
				do k=a1+1,a2
				  c2=resno(k)(1:1)
				  if(c2.ne.c1) then
					node_range(j,2,iseg)=k-1
					node_nseg(j)=node_nseg(j)+1
					jseg=node_nseg(j)
					node_range(j,1,jseg)=k
					node_range(j,2,jseg)=a2
					c1=c2
					iseg=jseg
				  end if
				end do
			    end if
			 end do
			end if
c!		end if
	  end if
	end do
c
c	reorder domains to start sequentially !
c
	do ix=1,ndom
		d(ix)=ix
		v(ix)=node_range(ix,1,1)
	end do
	call j_index_Qsort(1,ndom,ndom,d,v)
c
	write(32,730) nstrucdom,m-n,n,'seqtl resnos  (PDB resnos    )'
	m=0
	n=0
	do ix=1,ndom
		j=d(ix)
		if(node_type(j).eq.'*'.or.node_type(j).eq.'=') then
			m=m+1
			if(node_size(j).le.40) then
				n=n+1
				if(n.eq.26) write(*,*) 'can only count to 25!'
				k=mod(n,25)
				if(k.eq.0) k=25
				cd=arm(k)
				node_type(j)='='
			else
				if(m-n.eq.26) write(*,*) 'can only count to 25!'
				k=mod(m-n,25)
				if(k.eq.0) k=25
				cd=dom(k)
			end if
			call getclass(struc,node_range,node_nseg,node_size,j,
     $				na,nb,na1,na2,nb1,nb2,class)
			if(node_tert(j).gt.compact) then
				shape='compact    '
			else
				shape='nonglobular'
			end if
			write(32,740) cd,class,shape,node_size(j),node_nseg(j),
     $			  node_range(j,1,1),node_range(j,2,1),
     $			  resno(node_range(j,1,1)),resno(node_range(j,2,1))
			do i=2, node_nseg(j)
			  write(32,750) node_range(j,1,i),node_range(j,2,i),
     $			  resno(node_range(j,1,i)),resno(node_range(j,2,i))
			end do
		end if
	end do
c
c	32: terminator !
c
	write(32,760) '//'
c
c	caution: writetree1 overwrites all node_data !
c
	call writetree1(15,ndom,node_tau,node_child,node_range,node_nseg,
     $		node_parent,node_size,node_label,
     $		node_tert,node_nbp,node_type,resno,printdepth)
c
500	format(a4,a1)
510	format('read mol "',a8,'";')
520	format(20i4)
530	format(a50)
540	format(i4,f8.3,20i4)
550	format(a4,73x,a50)
610	format(5i4,f12.3,20i4)
620	format(a32,i4,2f8.3,20i4)
630	format(32x,4x,16x,20x,20i4)
640	format(a4,a1,1x,i2,a1,1x,a1,a4,a1,1x,a23)
650	format(a4,a1,1x,i2,a1,5x,a1,1x,l1,1x,a3,i4,1x,a5,' to ',a5)
660	format(a4,a1,1x,i2,a1,17x,a5,' to ',a5)
c
700	format('PDBID     ',a4)
710	format('COMPND    ',a60)
720	format('CHAIN     ',a1)
730	format('NDOM      ',i5,'   (',i5,' structural domains +',
     $   i5,' small units)'/
     $   '# no. class       shape       size  no. segments  ',a30)
740	format(1x,a2,3x,a10,2x,a11,i5,i9,7x,i5,' -',i5,2x,
     $   '(',a5,' to ',a5,')')
750	format(50x,i5,' -',i5,2x,'(',a5,' to ',a5,')')
760	format(a2)

c
	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine findstrucdoms(ndom,node_child,node_size,node_tau,
     $	node_tert,node_nbp,compact,taucutoff,lsplit,node_type,minsize)
	implicit none
	include 'sizes-puu.for'
	integer ndom,node_child(maxres,2)
	integer node_size(maxres),minsize,halfminsize
	real node_tau(maxres),compact,taucutoff
	logical lsplit
c
	integer i,j,k,m,n
	logical lstop
	real node_tert(maxres)
	integer node_nbp(maxres)
	integer node_stop(maxres),node_split(maxres),node_del(maxres)
	character node_type(maxres)
	logical lstop1
c
	halfminsize=minsize/2
	do i=1,ndom
		node_stop(i)=0
		node_split(i)=0
		node_type(i)='*'
		node_del(i)=0
	end do
	do i=1,ndom
		j=node_child(i,1)
		k=node_child(i,2)
		if(node_del(i).ne.0) j=0
		call testsplit(i,j,k,lstop,node_split,node_stop,node_size,
     $			node_tau,node_tert,node_nbp,compact,taucutoff,lsplit,
     $			minsize)
c
c		1nrd hack
c
		if((j.gt.0).and.(k.gt.0).and.(lstop)) then
		  if(node_size(k).lt.node_size(j)) then
			m=k
			n=j
		  else
			m=j
			n=k
		  end if
		  lstop1=.true.
		  if((node_size(m).lt.halfminsize).and.(node_size(n).gt.minsize)
     $			.and.(node_tert(m).lt.compact)
     $			.and.(node_tert(n).gt.compact))
     $			call testsplit(n,node_child(n,1),node_child(n,2),
     $				lstop1,node_split,node_stop,node_size,
     $				node_tau,node_tert,node_nbp,compact,taucutoff,
     $				lsplit,minsize)
		  if(.not.lstop1) then
c			write(*,*) ' 1nrd-hack override'
			node_split(j)=6
			node_split(k)=6
			lstop=.false.
		  end if
		end if
c
c		prune tree
c
		if(lstop) then
c			undo all offspring of terminal node i
			call delete(j,node_child,node_del,node_type)
			call delete(k,node_child,node_del,node_type)
		else
c			children j,k replace parent i
			node_del(i)=999
			node_type(i)='+'
		end if
	end do
c	write(*,*) 'findstrucdoms'
c	do i=1,ndom
c	  write(*,*) i,node_split(i),node_stop(i),(node_child(i,j),j=1,2),
c     $		node_del(i),node_type(i),node_tau(i),node_tert(i)
c	end do

	return
	end
c
c-------------------------------------------------------------------------------
c
	subroutine getclass(struc,node_range,node_nseg,node_size,
     $		inode,na,nb,na1,na2,nb1,nb2,class)
	implicit none
	include 'sizes-puu.for'
	character struc(maxres)
	integer node_range(maxres,2,maxnseg),node_nseg(maxres),inode
	integer na,nb,na1,na2,nb1,nb2,node_size(maxres)
c
	integer i,j,n,m
	character c
	character*10 class
	real pa,pb,pa1,pa2,pb1,pb2
c
	class='mixed     '
	na=0
	nb=0
	na1=0
	na2=0
	nb1=0
	nb2=0
	m=nint(float(node_size(inode))/ 2)
	n=0
	do i=1,node_nseg(inode)
		do j=node_range(inode,1,i),node_range(inode,2,i)
			c=struc(j)
			n=n+1
			if((c.eq.'H').or.(c.eq.'I').or.(c.eq.'G')) then
				na=na+1
				if(n.lt.m) then
					na1=na1+1
				else
					na2=na2+1
				end if
			else if((c.eq.'E').or.(c.eq.'B')) then
				nb=nb+1
				if(n.lt.m) then
					nb1=nb1+1
				else
					nb2=nb2+1
				end if
			end if
		end do
	end do
c
c	"alpha": %alpha > 40 & %beta < 5
c	"beta":  %alpha < 10 & %beta > 30
c	"alpha/beta": %alpha > 15 & %beta > 15
c	"mixed": %alpha+%beta < 20
c	"alpha+beta": one half "alpha", other half "beta"
c
	pa=float(na)/max(1,node_size(inode))
	pb=float(nb)/max(1,node_size(inode))
	pa1=float(na1)/max(1,m)
	pb1=float(nb1)/max(1,m)
	pa2=float(na2)/max(1,node_size(inode)-m)
	pb2=float(nb2)/max(1,node_size(inode)-m)
	if((pa1.gt.0.40).and.(pb1.lt.0.15).and.
     $			(pa2.lt.0.15).and.(pb2.gt.0.30))then
		class='alpha+beta'
	else if ((pa2.gt.0.40).and.(pb2.lt.0.15).and.
     $			(pa1.lt.0.15).and.(pb1.gt.0.30))then
		class='alpha+beta'
	else if((pa.gt.0.40).and.(pb.lt.0.15)) then
		class='all-alpha '
	else if((pa.lt.0.15).and.(pb.gt.0.30)) then
		class='all-beta  '
	else if((pa.gt.0.15).and.(pb.gt.0.15)) then
		class='alpha/beta'
	else if(pa+pb.lt.0.20) then
		class='mixed     '
	end if
c	write(*,500) pa,pb,pa1,pb1,pa2,pb2,class
500	format(6f8.3,a10)

	return
	end





c
c-------------------------------------------------------------------------------
c
	subroutine writetree1(iunit,ndom,node_tau,node_child,node_range,
     $		node_nseg,node_parent,node_size,
     $		node_label,node_tert,node_nbp,node_type,resno,printdepth)
	implicit none
	include 'sizes-puu.for'
	integer iunit,ndom,node_child(maxres,2)
	integer node_range(maxres,2,maxnseg)
	integer node_parent(maxres),node_size(maxres)
	real node_tau(maxres)
	integer node_nseg(maxres)
	character node_type(maxres)
	character*32 node_label(maxres)
	real node_tert(maxres)
	integer node_nbp(maxres),printdepth
	character*5 resno(0:maxres)
c
	integer i,j,l,maxl,dep(maxres),list(maxres),tsil(maxres)
	character c
	integer k1,k2
c
c	select nodes downto printdepth; renumber tree
c
	j=0
	do i=1,ndom
		tsil(i)=0
		c=node_type(i)
		if(c.eq.'+') dep(i)=1
		if(c.eq.'*'.or.c.eq.'=') dep(i)=0
		if(c.eq.'-') then
			dep(i)=dep(node_parent(i))-1
			if(node_type(node_parent(i)).eq.'=') dep(i)=dep(i)-1
		end if
		if(dep(i).ge.printdepth) then
			j=j+1
			list(j)=i
			tsil(i)=j
		end if
	end do
	ndom=j
	do i=1,ndom
		j=list(i)
		if(node_child(j,1).gt.0) then
			node_child(i,1)=tsil(node_child(j,1))
			node_child(i,2)=tsil(node_child(j,2))
		else
			node_child(i,1)=0
			node_child(i,2)=0
		end if
		node_type(i)=node_type(j)
		node_label(i)=node_label(j)
		node_tau(i)=node_tau(j)
		node_tert(i)=node_tert(j)
		node_nbp(i)=node_nbp(j)
		node_size(i)=node_size(j)
		node_nseg(i)=node_nseg(j)
		do l=1,node_nseg(i)
			node_range(i,1,l)=node_range(j,1,l)
			node_range(i,2,l)=node_range(j,2,l)
		end do
	end do
c
	maxl=0
	do i=1,ndom
		maxl=max(maxl,node_nseg(i))
	end do


c The following is a change made to make puu3 compilation work jong@ebi.ac.uk
        write(iunit,630) ndom
        do i=1,ndom
                l=node_nseg(i)
                if(node_tert(i).gt.9999) node_tert(i)=9999.999
                write(iunit,620) i,
     $ node_child(i,1),node_child(i,2),node_type(i),node_label(i),
     $ node_tau(i),node_tert(i),
     $ node_nbp(i),node_size(i),l
c         	write four number on one line, times l lines
          	do j=1,l
			k1=node_range(i,1,j)
			k2=node_range(i,2,j)
      			write(iunit,621) k1,k2,
     $ resno(k1),resno(k2)
          	end do
        end do
        write(iunit,640) '//'

c The above is newly written down one from Liisa's

c
600	format('>>>> ',a4,a1,i10)
610	format(5i4,1x,a1,a32,f8.3,f8.3,20i4)
620     format(i5,1x,2i4,1x,a1,a32,2f9.3,i8,2i7)
621	format(2i10,2(5x,a5))
630	format('NUNITS    ',i5,
     $ /'# no. children type+label                     tau-squared ',
     $ 'compactness H-bonds size segments residues',8x,
     $ '(PDB residues   )')
640	format(a2)
c
	return
	end

c
c-------------------------------------------------------------------------------
c
	subroutine testsplit(i,j,k,lstop,node_split,node_stop,node_size,
     $		node_tau,node_tert,node_nbp,compact,taucutoff,
     $		lsplit,minsize)
	implicit none
	include 'sizes-puu.for'
	integer i,j,k,node_nbp(maxres),minsize
	logical lstop,lsplit
	integer node_split(maxres),node_stop(maxres),node_size(maxres)
	real node_tau(maxres),node_tert(maxres),compact,taucutoff
c
		lstop=.true.
		if((j.eq.0).or.(k.eq.0)) then
			lstop=.true.
			node_stop(i)=1
		else if(node_size(i).lt.minsize) then
			lstop=.true.
			node_stop(i)=-1
		else if(node_tau(i).gt.taucutoff) then
			lstop=.false.
			node_split(j)=1
			node_split(k)=1
c
c		lsplit=.true. ignores hbonds
c		lsplit=.false. precludes splitting sheets
c
		else if(node_nbp(i).gt.0.and..not.lsplit) then
			lstop=.true.
			node_stop(i)=2
		else if(node_tert(i).lt.compact) then
			lstop=.true.
			node_stop(i)=3
		else if((node_tert(j).le.compact).and.(node_tert(k).le.compact))
     $		  then
			lstop=.true.
			node_stop(i)=4
		else if((node_tert(j).gt.compact).and.(node_tert(k).gt.compact))
     $		  then
			lstop=.false.
			node_split(j)=2
			node_split(k)=2
		end if

	return
	end
c-------------------------------------------------------------------------------
c
	subroutine delete(knode,node_child,node_del,node_type)
	implicit none
	include 'sizes-puu.for'
	integer knode, node_child(maxres,2),node_del(maxres)
	character node_type(maxres)
c
	integer nstack,stack(maxres),i,j,inode
c
c	sink whole line from knode to depth 999
c
	nstack=0
	call putstack(knode,nstack,stack)
	do while (nstack.gt.0)
		inode=stack(nstack)
		node_type(inode)='-'
		node_del(inode)=999
		nstack=nstack-1
		i=node_child(inode,1)
		j=node_child(inode,2)
		call putstack(i,nstack,stack)
		call putstack(j,nstack,stack)
	end do

        return
        end
c
c------------------------------------------------------------------------------
c
	subroutine putstack(i,nstack,stack)
	implicit none
	include 'sizes-puu.for'
	integer nstack,stack(maxres),i

	if(i.gt.0) then
		nstack=nstack+1
		stack(nstack)=i
	end if

        return
        end
c
c------------------------------------------------------------------------------
c
