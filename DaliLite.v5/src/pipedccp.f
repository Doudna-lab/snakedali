c This module/program is part of DaliLite (c) L. Holm 1999
c
        program refinput
c
c       derivative of pipe: input from DCCP file
c
        implicit none
        include 'parsizes.for'
        character*5 cd1,cd2
        character*10 line
        real score
        integer nblock,l1(maxres),l2(maxres),r1(maxres),r2(maxres)
        integer i
        integer ide,lali
        real rmsd
        logical lfirst
c
c
c
        lfirst=.true.
10      read(*,600,err=10,end=19) line,score,rmsd,lali,ide,nblock,
     $          cd1,cd2
        if(line(1:9).eq.'DCCP   1 '.or.line(2:10).eq.'DCCP   1 ') then
100             read(*,605,end=19) line
                if(line(1:9).ne.'alignment'.and.line(2:10).ne.
     $                  'alignment')goto 100
                read(*,610,err=10) (l1(i),r1(i),i=1,nblock)
                read(*,610,err=10) (l2(i),r2(i),i=1,nblock)
        end if
c
c       output
c
        if(.not.lfirst) write(*,500) 'END  '
        lfirst=.false.
        if(cd1(5:5).eq.'_') cd1(5:5)=' '
        write(*,500) cd1

        if(cd2(5:5).eq.'_') cd2(5:5)=' '
        write(*,500) cd2,'*'
        write(*,*) nblock
        write(*,*) (l1(i),r1(i),i=1,nblock)
        write(*,*) (l2(i),r2(i),i=1,nblock)
        goto 10
19      close(98)
        write(*,500) 'END  '
c        write(*,500) 'END  '

500     format(a5,a1,a1)
600     format(1x,a9,f8.1,f4.1,i4,16x,i4,3x,i4,16x,a5,1x,a5)
601     format(1x,a9,f8.1,f4.1,i4,f8.1,8x,i4,3x,i4,16x,a5,1x,a5)
605     format(1x,a9)
610     format(8(i4,2x,i4))
!DCCP   1    222.7-1.0  55     3.3          -1     55                9rnt  1bsaA
        end
c
c------------------------------------------------------------------------------
c

