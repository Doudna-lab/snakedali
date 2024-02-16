c
c utilities used in many programs
c 
c no include files here


c
c-------------------------------------------------------------------------------
c
        function constructfilnam(cd,path,ext)
        implicit none
        character*(*) cd
        character*80 path
        character*(*) ext
        character*80 constructfilnam
c
        integer l1,l2,l3,i
c
        constructfilnam=''
        l1=len_trim(cd)
        if(l1.lt.1) return
        i=0                             ! need to parse path for l2
        do while(path(i+1:i+1).le.' ')  ! remove leading blanks
                i=i+1
        end do
        do while(path(i+1:i+1).ne.' ')  ! stop at first right-side blank
                i=i+1
        end do
        l2=i
        l3=len_trim(ext)
        !write(*,*) 'constructfilnam: |',cd,'|',path,ext,l1,l2,l3
        do i=1,80
                constructfilnam(i:i)=' '
        end do
        constructfilnam(1:l2)=path
        constructfilnam(l2+1:l2+l1)=cd
        if(cd(l1:l1).eq.' ') constructfilnam(l2+l1:l2+l1)='_'
        constructfilnam(l2+l1+1:l2+l1+l3)=ext
        !write(*,*) constructfilnam(1:80)

        return
        end
c
c-----------------------------------------------------------------------
c
        subroutine getdalidata(dalidatpath,cd1,iunit,nres,seq,stru,
     $         compnd,ca,resno)
c
c       reads DAT file
c
        implicit none
 	integer maxres,maxdomseg
	parameter(maxres=17000,maxdomseg=320) !!! better be consistent with main !!!
        integer iunit,nres,nsse,acc(maxres)
        character*80 dalidatpath
        character*80 filnam, compnd,constructdatfilnam
        character seq(maxres),stru(maxres),node_type(maxres)
        character secstr(maxres)
        real ca(3,maxres)
        character*5 resno(0:maxres),cd1
        integer l,m,i,j,segstart(maxres),segend(maxres),ndom,nh,ne,idom
        integer node_child(2,maxres),domns(maxres)
        integer domsize(maxres),domseglist(maxdomseg,maxres)
c
        filnam=constructdatfilnam(cd1,dalidatpath,'.dat')
        open(iunit,file=filnam,status='old',err=99)
c       parse all information in DAT-file, return subset
        read(iunit,500) nres, nsse, nh, ne,(secstr(i),i=1,nsse)
        ! segment.start.end
        if (nsse.gt.0) read(iunit,510) (segstart(i),segend(i),i=1,nsse)
        ! ca.x,y,z
10      read(iunit,520,err=10) ((ca(j,i),j=1,3),i=1,nres)
        ! sse-domains
        read(iunit,500) ndom
        do idom=1,ndom
           read(iunit,530) i,node_type(i),(node_child(j,i),j=1,2),
     $          domns(i),(domseglist(j,i),j=1,
     $          min(maxdomseg,domns(i)))
        end do
        ! residue-ranges of domains
        read(iunit,500) ndom
        do idom=1,ndom
          read(iunit,530) i,node_type(i),(node_child(j,i),j=1,2),
     $          domsize(i),domns(i),(domseglist(j,i),
     $          j=1,min(maxdomseg,domns(i)))
        end do
        ! keywords
        read(iunit,540) (stru(i),i=1,nres)
        read(iunit,540) (seq(i),i=1,nres)
        read(iunit,550) compnd
        do i=1,nres
                read(iunit,560) resno(i),acc(i)
        end do
        
        close(iunit)
        return
        
99      write(*,*) 'ERROR opening filnam=',filnam
c	stop "FATAL ERROR in getdata: Can't open DAT file "

500     format(10x,4i5,2x,200a1)
510     format(10x,2i10)
520     format(10f8.1)
530     format(i4,1x,a1,1x,1000i4)
540     format(11x,6000a1) ! <maxres>a1
550     format(11x,a80)
c560     format(13x,a5,i5)
560 	format(12x,a5,1x,i5)

        return
        end
c
c----------------------------------------------------------------------
c
        function constructdatfilnam(cd,path,ext)
        implicit none
        character*5 cd
        character*80 path
        character*(*) ext
        character*80 constructdatfilnam
c
        integer l1,l2,l3,i
c
        if(cd(5:5).eq.' ') cd(5:5)='_'
        i=0
	l1=len(cd)
	do while(i.lt.l1.and.cd(i+1:i+1).ne.' ') ! remove blank chainid
		i=i+1
	end do
	l1=i
c
        i=0                             ! need to parse path for l2
        do while(path(i+1:i+1).eq.' ')  ! remove leading blanks
                i=i+1
        end do
        do while(path(i+1:i+1).ne.' ')  ! stop at first right-side blank
                i=i+1
        end do
        l2=i
        l3=len(ext)
c        write(*,*) 'constructdatfilnam: |',cd,'|',path,ext,l1,l2,l3
        do i=1,80
                constructdatfilnam(i:i)=' '
        end do
        constructdatfilnam(1:l2)=path
        constructdatfilnam(l2+1:l2+l1)=cd
        constructdatfilnam(l2+l1+1:l2+l1+l3)=ext
        if(cd(5:5).eq.'_') cd(5:5)=' '

        return
        end
c
c-----------------------------------------------------------------------
c
	subroutine getdssp(filnam,chainid,iunit,maxres,header,compnd,
     $		source,author,nres,nchain,nss,nssintra,nssinter,area,
     $		resno,seq,struc,acc,dsspinfo,dsspresno,fsspresno,
     $		tco,kappa,alfa,phi,psi,ca,ierr)
c
c	reads information from a DSSP file
c	extracts residues with chainid (* is wildcard)
c	chain breaks (! residues) are ignored
c
c	use a frontend to transport wanted subset of information
c
c	input:	filnam
c		chainid
c		iunit
c		maxres 
c
c	output:	header, compnd, source, author -- echo PDB cards
c		nres
c		nchain
c		nss, nssintra, nssinter -- SS-bridges
c		area -- total accessible surface area
c		resno -- PDB residue number string with appended insertion
c			 indicator, e.g. 1ACX 63A
c		seq -- one-letter sequence
c		struc -- one-letter secondary structure
c		acc -- accessible surface area per residue
c		dsspinfo -- structure string with renumbered beta-sheet resnos
c				* marks H-bonds to outside given chain
c		dsspresno, fsspresno -- local work arrays
c		tco,kappa,alfa,phi,psi -- backbone conformational angles
c		ca -- CA coordinates
c		ierr = 0 -- normal
c		     = -1 -- file not found
c		     = -2 -- too many residues, chain truncated to maxres
c
	implicit none
	integer iunit,maxres,nres,nchain,nss,nssintra,nssinter,acc(maxres)
	character*80 filnam,header,compnd,source,author
	character chainid,seq(maxres),struc(maxres)
	real tco(maxres),kappa(maxres),alfa(maxres),phi(maxres),
     $		psi(maxres),ca(3,maxres)
	character*5 resno(0:maxres)
	character*33 dsspinfo(maxres)
	integer area,ierr,dsspresno(0:maxres),fsspresno(maxres)
c
	integer i,ntot,bp1,bp2,linecount
	character*136 line
	logical lheader,lcompnd,lsource,lauthor
	character se,c
	character*4 rno,str4
c
c	initialize
c
	ierr=-1
	do i=1,80
		header(i:i)=' '
		compnd(i:i)=' '
		source(i:i)=' '
		author(i:i)=' '
	end do
	nres=0
	ntot=99999
	nchain=0
	nss=0
	nssintra=0
	nssinter=0
	area=0.0
	lheader=.false.
	lcompnd=.false.
	lsource=.false.
	lauthor=.false.
c
	open(iunit,file=filnam,status='old',err=999)
	ierr=0
c
c	read header information until start of data ('  # ' line)
c
100	read(iunit,500,end=999) line
	if((line(1:6).eq.'HEADER').and.(.not.lheader)) then
		header=line(1:80)
		lheader=.true.
	end if
	if((line(1:6).eq.'COMPND').and.(.not.lcompnd)) then
		compnd=line(1:80)
		lcompnd=.true.
	end if
	if((line(1:6).eq.'SOURCE').and.(.not.lsource)) then
		source=line(1:80)
		lsource=.true.
	end if
	if((line(1:6).eq.'AUTHOR').and.(.not.lauthor)) then
		author=line(1:80)
		lauthor=.true.
	end if
	if(line(19:42).eq.'TOTAL NUMBER OF RESIDUES') then
		read(line(1:),*,err=100) ntot,nchain,nss,nssintra,nssinter
	end if
c	skip nhb, nhbp
	if(line(1:4).ne.'  # ') goto 100
c
c	read residue-information; do not read beyond ntot lines
c
	linecount=0
200	read(iunit,500,end=299) line
	if(line(14:14).ne.'!') linecount=linecount+1
	if(linecount.gt.ntot) then
		write(*,*) ' Skipping lines after ',ntot
		goto 299
	end if
c	skip breaks
	se=line(14:14)
	if(se.eq.'!') goto 200
c	check maxres
	read(line(1:5),*) i
	if((i.ge.maxres).or.(nres.eq.maxres)) then
		ierr=-2
		write(*,*) ' Too many residues in getdssp -- skipping ',
     $			i,line(1:38)
		goto 299
	end if
c	check chainid
	c=line(12:12)
	if(chainid.ne.'*') then
	  if(chainid.ne.c) then
		fsspresno(i)=0
		goto 200
	  end if
	end if
c	new residue
	nres=nres+1
	seq(nres)=se
	struc(nres)=line(17:17)
	resno(nres)(1:1)=c
	read(line(35:38),*) acc(nres)
c	append residue-insertion code + rightadjust
	rno=line(7:10)
	c=line(11:11)
	if(c.ne.' ') then
		do i=2,4
			rno(i-1:i-1)=rno(i:i)
		end do
		rno(4:4)=c
	end if
	resno(nres)(2:5)=rno
	read(line,510) tco(nres),kappa(nres),alfa(nres),phi(nres),
     $		psi(nres),ca(1,nres),ca(2,nres),ca(3,nres)
c	sheet info
	read(line(1:5),*) i
	dsspresno(nres)=i
	fsspresno(i)=nres
	dsspinfo(nres)=line(6:38)
c
	goto 200
c
c	renumber sheet records; -99 if H-bond partner is in other chain
c
299	do i=1,nres
		read(dsspinfo(i)(21:24),*) bp1
		read(dsspinfo(i)(25:28),*) bp2
c		write(*,*) 'dssp:',i,bp1,bp2,dsspresno(i)
		if(bp1.gt.0) then
			bp1=fsspresno(bp1)
			if(bp1.gt.0) then
				dsspinfo(i)(21:24)=str4(bp1)
			else
				dsspinfo(i)(21:24)=' -99'
			end if
		end if
		if(bp2.gt.0) then
			bp2=fsspresno(bp2)
			if(bp2.gt.0) then
				dsspinfo(i)(25:28)=str4(bp2)
			else
				dsspinfo(i)(25:28)=' -99'
			end if
		end if
	end do

500	format(a136)
c510     format(75x,f8.3,4f6.1,3f7.1) ! dsspcmbi format is wider
510     format(83x,f8.3,4f6.1,3f7.1)

999	close(iunit)

	return
	end
c
c-----------------------------------------------------------------------
c
	function onelettercode(aaa)
	implicit none
	character*3 aaa
	character onelettercode,a
c
c	converts residuename to one-letter-code
c
c	COMPATIBLE WITH classic DSSP
c	'-' is to be rejected !
c       accept MSE for compatibility with dsspCMBI
c
	character*52 aasymbol
	character*156 aminoacid
	data aasymbol
     $          /'ARNDCEQGHILKMFPSTWYVBZXXXXXXXXXXXXXXXX--CCCCIPPPW-XM'/
c
	integer i,k
c
	aminoacid(1:60)=
     $	'ALAARGASNASPCYSGLUGLNGLYHISILELEULYSMETPHEPROSERTHRTRPTYRVAL'
	aminoacid(61:120)=
     $	'ASXGLXACDALBALIABUAROBASBETHSEHYPHYLORNPCASARTAUTHYUNKACEFOR'
	aminoacid(121:156) = 'CYHCSHCSSCYXILUPRZPR0CPRTRYHOHSEMMSE'
c
	a='-'
	do i=1,154,3
		if(aminoacid(i:i+2).eq.aaa) then
			k=(i+2)/3
			a=aasymbol(k:k)
			goto 19
		end if
	end do
19	onelettercode=a
	return
	end

c
c-----------------------------------------------------------------------
c
	subroutine getcoor(iunit,maxres,maxatm,filnam,chainid,
     $		header,compnd,source,author,
     $		xyz,natom,atmnam,atmres,nres,resno,resnam,ierr,lhetatm)
c
c	COMPATIBLE WITH DSSPcmbi
c
c	if find an ENDMDL record, quit reading (NMR structures) -- 17-Aug-1993
c
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
c	output:	header, compnd, source, author -- echo PDB cards
c		nres,natom
c		resno -- PDB residue number string with appended insertion
c			 indicator, e.g. 1ACX 63A
c		resnam -- 'ALA'
c		atmnam -- 'CA '
c		atmres -- residue index of atom
c		bvalue -- real
c		xyz -- atomic coordinates
c		ierr = 0 -- normal
c		     = -1 -- file not found
c		     = -2 -- too many residues, chain truncated to maxres
c
c	acceptable residue: onelettercode.ne.'-'
c	new residue: different resno from previous one
c	ignore atoms which start 'H' or 'D'
c	ignore residues with incomplete backbone (N,CA,C,O)
c
	integer iunit,maxres,maxatm
	character*80 filnam,header,compnd,source,author
	character chainid
	real xyz(3,maxatm)
	integer natom,atmres(maxatm),nres,ierr
	character*3 atmnam(maxatm),resnam(maxres)
	character*5 resno(0:maxres)
	logical lhetatm
c
	character onelettercode
c
	character*80 line
	integer i,nig
	character*3 aaa
	character*6 oldresname,resname
	character c,d,altloc
	character*4 rno
	logical ln,lc,lca,lo,lsplit,lsource,lauthor,lcompnd,lheader
	logical lchain

c
c	initialize
c
	ierr=-1
	do i=1,80
		header(i:i)=' '
		compnd(i:i)=' '
		source(i:i)=' '
		author(i:i)=' '
	end do
	oldresname='!@#$'
	nig=0
	nres=0
	natom=0
	resno(0)='     '
	lheader=.false.
	lcompnd=.false.
	lchain=.false.
	lsource=.false.
	lauthor=.false.
	open(iunit,file=filnam,status='old',err=999)
	ierr=0
c
c	read data
c
100	read(iunit,500,end=199) line
	if(line(1:6).eq.'ENDMDL') then
		write(*,*) 'detected ENDMDL record in NMR entry - skip rest '
		goto 199
	end if
	if((line(1:6).eq.'HEADER').and.(.not.lheader)) then
		header=line(1:80)
		lheader=.true.
	end if
COMPND    MOL_ID: 1;
COMPND   2 MOLECULE: PUTATIVE UNCHARACTERIZED PROTEIN ECU08_1170;
COMPND   3 CHAIN: A;
COMPND   4 SYNONYM: TRM112 ACTIVATOR PROTEIN;
COMPND   5 ENGINEERED: YES;
COMPND   6 MOL_ID: 2;
COMPND   7 MOLECULE: N6 ADENINE SPECIFIC DNA METHYLASE;
COMPND   8 CHAIN: B;
COMPND   9 SYNONYM: METHYLTRANSFERASE SUPERFAMILY, MTQ2 CATALYTIC SUBUNIT;
COMPND  10 ENGINEERED: YES
	if((line(1:6).eq.'COMPND').and.(.not.lcompnd).and.
     $          (.not.lchain)) then
		compnd=line(1:80)
		lcompnd=.true.
	end if
	if((line(1:6).eq.'COMPND').and.(line(12:17).eq.'CHAIN:').and.
     $		(.not.lchain)) then
		c=line(19:19)
		if(c.eq.chainid) lchain=.true.
	end if
	if((line(1:6).eq.'SOURCE').and.(.not.lsource)) then
		source=line(1:80)
		lsource=.true.
	end if
	if((line(1:6).eq.'AUTHOR').and.(.not.lauthor)) then
		author=line(1:80)
		lauthor=.true.
	end if
	if(line(1:6).eq.'ATOM  '.or.(line(1:6).eq.'HETATM'
     $          .and.lhetatm)) then
		c=line(22:22)
		if(chainid.ne.'*') then
			if(chainid.ne.c) goto 109
		end if
		resname=line(22:27)
		altloc=line(17:17)
c		accept altloc=1 to make 2lbd compatible with dsspcmbi
		if(altloc.ne.' '.and.altloc.ne.'A'.and.altloc.ne.'1') then
			write(*,*) 'Alternate location indication is not ',
     $				'blank or A, ignore atom:',line(7:27)
			goto 109
		end if
		aaa=line(18:20)
		c=onelettercode(aaa)
		if(c.eq.'-'.and.line(1:6).ne.'HETATM') then
			write(*,*) ' ignoring nonstandard residue ',aaa
			goto 109
		end if
		d=line(14:14)
		if((d.eq.'H').or.(d.eq.'D')) then
			nig=nig+1
			goto 109
		end if
		if(resname.ne.oldresname) then
			if(nres.gt.0) then
			  if(ln.or.lo.or.lc.or.lca) then
	if(.not.lhetatm)write(*,520) resnam(nres),oldresname,nres
				i=natom+1
				do while((atmres(i-1).eq.nres).and.(i.gt.0))
					i=i-1
					nig=nig+1
				end do
				natom=i
				nres=nres-1
			  end if
			end if
			ln=.true.
			lo=.true.
			lc=.true.
			lca=.true.
			if(nres.eq.maxres) then
				ierr=-2
				write(*,*) ' Too many residues ',nres
				goto 199
			end if
			nres=nres+1
			resnam(nres)=aaa
			rno=line(23:26)
			c=line(27:27)
			if(c.ne.' ') then
				do i=2,4
					rno(i-1:i-1)=rno(i:i)
				end do
				rno(4:4)=c
			end if
			resno(nres)(2:5)=rno
			resno(nres)(1:1)=line(22:22)
			oldresname=resname
		end if
		if(natom.eq.maxatm) then
			ierr=-2
			write(*,*) ' Too many atoms ',natom
			goto 199
		end if
		natom=natom+1
		atmnam(natom)=line(14:16)
		if(atmnam(natom).eq.'N  ') ln=.false.
		if(atmnam(natom).eq.'CA ') lca=.false.
		if(atmnam(natom).eq.'C  ') lc=.false.
		if(atmnam(natom).eq.'O  ') lo=.false.
		atmres(natom)=nres

		read(line(31:54),*,err=109) xyz(1,natom),xyz(2,natom),
     $			xyz(3,natom)
109	end if
c
	goto 100
c
c	need to check backbone of last residue especially
c
199	if(nres.gt.0) then
	  if(ln.or.lo.or.lc.or.lca) then
		if(.not.lhetatm)write(*,520) resnam(nres),resname,nres
		i=natom+1
		do while((atmres(i-1).eq.nres).and.(i.gt.0))
			i=i-1
			nig=nig+1
		end do
		natom=i
		nres=nres-1
	  end if
	end if
	close(iunit)
	write(*,510) nres,natom,nig

500	format(a80)
510	format(/'Total number of residues read in:     ',i10,
     $       /'Total number of atoms read in:        ',i10,
     $       /'Total number of ignored atom records: ',i10)
520	format(/'Ignoring residue with incomplete backbone ',a3,1x,a6,i5)

999	return
	end

c
c-----------------------------------------------------------------------
c

        subroutine j_index_Qsort(i,j,n,d,v)
c
c       non-recursive Quicksort
c       sorts integer array v which has n elements from smallest to largest
c       d contains indices to v
c
        implicit none
        integer i,j,n,d(n),v(n)
c
        integer p,top,bottom,stack(100), nstack, partition

        stack(1)=j
        stack(2)=i
        nstack=2
        do while(nstack.ne.0)
                top=stack(nstack)
                bottom=stack(nstack-1)
                nstack=nstack-2
                do while(top.lt.bottom)
                        p=partition(top, bottom, n, d, v)
                        if((p-top).gt.(bottom-p)) then
                                stack(nstack+1)=p-1
                                stack(nstack+2)=top
                                top=p+1
                                nstack=nstack+2
                                if(nstack.gt.100) return 
c     $  stop ' Stack overflow in Qsort '
                        else
                                stack(nstack+1)=bottom
                                stack(nstack+2)=p+1
                                bottom=p-1
                                nstack=nstack+2
                                if(nstack.gt.100) return
c    $  stop ' Stack overflow in Qsort '
                        end if
                end do
        end do

        return
        end

c ======================================================================
        function partition(i,j,n,d,v)
        implicit none
        integer i,j,n,partition, d(n)
        integer v(n)
c
        integer upper, lower, save

        upper=i
        lower=j
        save=d(i)

        do while (upper.ne.lower)
                do while((upper.lt.lower).and.(v(save).le.v(d(lower))))
                        lower=lower-1
                end do
                if(upper.ne.lower) d(upper)=d(lower)
                do while((upper.lt.lower).and.(v(save).ge.v(d(upper))))
                        upper=upper+1
                end do
                if(upper.ne.lower) d(lower)=d(upper)
        end do
        d(upper)=save
        partition=upper

        return
        end
c
c-----------------------------------------------------------------------
c
	function onetothree(aa)
	implicit none
	character aa
	character*3 onetothree

	onetothree='UNK'
	if(aa.eq.'A') then
		onetothree='ALA'
	else if(aa.eq.'P') then
		onetothree='PRO'
	else if(aa.eq.'S') then
		onetothree='SER'
	else if(aa.eq.'C') then
		onetothree='CYS'
	else if(aa.eq.'T') then
		onetothree='THR'
	else if(aa.eq.'V') then
		onetothree='VAL'
	else if(aa.eq.'I') then
		onetothree='ILE'
	else if(aa.eq.'L') then
		onetothree='LEU'
	else if(aa.eq.'D') then
		onetothree='ASP'
	else if(aa.eq.'N') then
		onetothree='ASN'
	else if(aa.eq.'H') then
		onetothree='HIS'
	else if(aa.eq.'F') then
		onetothree='PHE'
	else if(aa.eq.'Y') then
		onetothree='TYR'
	else if(aa.eq.'W') then
		onetothree='TRP'
	else if(aa.eq.'M') then
		onetothree='MET'
	else if(aa.eq.'E') then
		onetothree='GLU'
	else if(aa.eq.'Q') then
		onetothree='GLN'
	else if(aa.eq.'K') then
		onetothree='LYS'
	else if(aa.eq.'R') then
		onetothree='ARG'
	else if(aa.eq.'G') then
		onetothree='GLY'
	else if((aa.ge.'a').and.(aa.le.'z')) then
		onetothree='CYS'
	end if

	return
	end
c
c----------------------------------------------------------------------
c
	function str4(i)
c
c	converts positive integer to 4-character string
c
	implicit none
	integer i
	character*4 str4,s
	character t(0:9)
	data t/'0','1','2','3','4','5','6','7','8','9'/
c
	integer k,i1,i2,i3,i4
c
	do k=1,3
		s(k:k)=' '
	end do
	s(4:4)='0'
	i1=mod(i,10)
	i2=(mod(i-i1,100))/10
	i3=(mod(i-i2-i1,1000))/100
	i4=(mod(i-i3-i2-i1,10000))/1000
	if(i.ge.0) s(4:4)=t(i1)
	if(i.ge.10) s(3:3)=t(i2)
	if(i.ge.100) s(2:2)=t(i3)
	if(i4.le.9) then
		if(i.ge.1000) s(1:1)=t(i4)
	else
		s(1:1)='*'
	end if
	str4=s

	return
	end
c
c-----------------------------------------------------------------------
c
        subroutine u2xyzrotate(u,
     $          rotateaxis,x_degrees,y_degrees,z_degrees)
        implicit none
        ! input
        real u(3,3)
        ! result
        character rotateaxis(3)
        real x_degrees,y_degrees,z_degrees
        ! internal
        real pi,tolerance
        parameter(pi=3.14159265,tolerance=1e-6)
        real cosa,sina,cosb,sinb,cosg,sing
c        real u11,u12,u13,u21,u22,u23,u31,u32,u33
        real ssq,s,bestssq,besta,bestb,bestg
        integer i,j,k
        real a(4),b(2),g(4)
        integer na,nb,ng

        ! case |u13|<1
        if(abs(u(1,3)).lt.1) then
                rotateaxis(1)='x'
                rotateaxis(2)='y'
                rotateaxis(3)='z'
                na=4
                a(1)=asin(u(2,3)/sqrt(1-u(1,3)*u(1,3)))
                a(2)=pi-a(1)
                a(3)=asin(-u(2,3)/sqrt(1-u(1,3)*u(1,3)))
                a(4)=pi-a(3)
                nb=2
                b(1)=asin(u(1,3))
                b(2)=pi-b(1)
                ng=4
                g(1)=acos(u(1,1)/sqrt(1-u(1,3)*u(1,3)))
                g(2)=-g(1)
                g(3)=acos(-u(1,1)/sqrt(1-u(1,3)*u(1,3)))
                g(4)=-g(3)
        else
                rotateaxis(1)='y'
                rotateaxis(2)='z'
                rotateaxis(3)='x'
                na=4
                a(1)=asin(-u(1,2))
                a(2)=pi-a(1)
                a(3)=asin(u(1,2))
                a(4)=pi-a(3)
                nb=2
                b(1)=pi/2
                b(2)=-pi/2
                ng=2
                g(1)=asin(u(2,1))
                g(2)=pi-g(1)
        end if
        ! test sign combinations
        bestssq=99999.9
        besta=0.0
        bestb=0.0
        bestg=0.0
        do i=1,na
          do j=1,nb
            do k=1,ng
                s=ssq(u,rotateaxis,a(i),b(j),g(k))
                if(s.lt.bestssq)then
                                bestssq=s
                                besta=a(i)
                                bestb=b(j)
                                bestg=g(k)
        if(s.lt.1e-6) goto 100 ! quit on first good fit
                end if
        if(s.lt.1e-6) write(*,*) '# fit:',besta,bestb,bestg,bestssq
            end do
          end do
        end do
        ! output result
100     x_degrees=besta/pi*180.0
        if(x_degrees.gt.180) x_degrees=x_degrees-360
        y_degrees=bestb/pi*180.0
        if(y_degrees.gt.180) y_degrees=y_degrees-360
        z_degrees=bestg/pi*180.0
        if(z_degrees.gt.180) z_degrees=z_degrees-360

        write(*,*) '# final result:',besta,bestb,bestg,bestssq
        if(bestssq.gt.tolerance) write(*,*) '# warning: bad fit!'
c        write(*,*) '# in degrees:',x_degrees,y_degrees,z_degrees

        return
        end

        function ssq(u,rotateaxis,a,b,g)
        implicit none
        real u(3,3),a,b,g,ssq
        real v(3,3),err(3,3)
        character rotateaxis(3)

        real sina,cosa,sinb,cosb,sing,cosg,s
        integer i,j
        logical lverbose
        parameter(lverbose=.false.)

        sina=sin(a)
        cosa=cos(a)
        sinb=sin(b)
        cosb=cos(b)
        sing=sin(g)
        cosg=cos(g)

        if(rotateaxis(1).eq.'y') then
                ! construct v=[y_rotate(b)][z_rotate(g),[x_rotate(a)]
                v(1,1)=cosb*cosg
                v(1,2)=-cosa*cosb*sing+sina*sinb
                v(1,3)=sina*cosb*sing+cosa*sinb
                v(2,1)=sing
                v(2,2)=cosa*cosg
                v(2,3)=-sina*cosg
                v(3,1)=-sinb*cosg
                v(3,2)=cosa*sinb*sing+sina*cosb
                v(3,3)=-sina*sinb*sing+cosa*cosb
        else if(rotateaxis(1).eq.'x') then
                ! construct v=[x_rotate(a)][y_rotate(b),[z_rotate(g)]
                v(1,1)=cosb*cosg
                v(1,2)=-cosb*sing
                v(1,3)=sinb
                v(2,1)=sina*sinb*cosg+cosa*sing
                v(2,2)=-sina*sinb*sing+cosa*cosg
                v(2,3)=-sina*cosb
                v(3,1)=-cosa*sinb*cosg+sina*sing
                v(3,2)=cosa*sinb*sing+sina*cosg
                v(3,3)=cosa*cosb
        else
                stop '# unknown rotateaxis in ssq'
        end if

        ssq=0.0
        do i=1,3
                do j=1,3
                        err(i,j)=(u(i,j)-v(i,j))*(u(i,j)-v(i,j))
                        ssq=ssq+err(i,j)
                end do
        end do
        ! print difference
        if(lverbose) then
                write(*,*) '# a,b,g:',a,b,g,ssq
                do i=1,3
                        write(*,500) (err(i,j),j=1,3),(v(i,j),j=1,3)
                end do
        end if

500     format('#ssq',3f8.5,8x,3f8.5)

        return
        end
c
c-----------------------------------------------------------------------
c

