        integer maxstack2,maxdom,maxres0,maxres1,maxres2,exdim
        parameter (maxres1=4000,maxres2=4000,maxres0=maxres2/10)
        parameter (exdim=1000000)
        integer maxstack,maxseg,maxatm,maxres,maxgrid,maxrot
        parameter (maxstack=10000,maxstack2=10000)
        parameter (maxseg=200, maxdom=maxseg*2)
        ! maxprot is number of PDB chains
        integer maxprot,maxpair, maxnoccs
        integer weight(1000)
        integer scoretable(160,400)
        integer infinit,nul
        parameter (infinit=1e7,nul=-99)
        parameter (maxprot=1000000)
        parameter (maxatm=9999)
        parameter (maxres=maxres2)
        parameter (maxgrid=20)
        parameter (maxrot=500)
        parameter (maxpair=50000)
        parameter (maxnoccs=300000)
        integer rankcount,callcount,testcount

        integer boxdim
        parameter(boxdim=501)

        common weight,scoretable,rankcount,callcount,testcount

        logical lskipflag
        common lskipflag
