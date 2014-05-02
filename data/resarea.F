      subroutine RESAREA(IATMR1,IATMR2,NRES,ATMNAM,NATMS,ACCSS,RESNAM,
     +                   RESNUM,RESACC,NUMACC,NACCESS,RACCSS,MAXRES,
     +                   ACCTYP,OUTHET)
C
C CALCULATE ABSOLUTE AND PERCENTAGE ACCESSIBILITIES FOR:
C                     WHOLE RESIDUE
C                     NON POLAR SIDECHAIN ATOMS (INCLUDING CA)
C                     POLAR SIDECHAIN ATOMS
C                     TOTAL SIDECHAIN ATOMS
C                     TOTAL MAINCHAIN ATOMS (EXCLUDING CA)
C
C SUMMATION FOR SIDE CHAIN includeS CA: GLY then HAS A SIDECHAIN AND THERE
C IS NO SPECIAL CASE.
C
      parameter (NRT=20)
      integer IATMR1(NRES), IATMR2(NRES), I
      real RACCSS(10,MAXRES),ACCSS(NATMS)
      real ATOTAL(NRT,2), AMAIN(NRT,2), ASIDE(NRT,2)
      real APOLSIDE(NRT,2), ANPOLSIDE(NRT,2)
      character*4 ATMNAM(NATMS)
      character*3 RESNAM(NRES), RESACC(MAXRES)
      character*5 NUMACC(MAXRES), RESNUM(NRES)
      character*(*) ACCTYP*(*)
      logical OUTHET
C
C --- New standard accessibilities by Simon Hubbard (Dec. 1989) ------------
C
C TOTAL CONTACT AREA IN FULLY EXTENDED open CHAIN FORM
C
      data (ATOTAL(I,1),I=1,NRT)/ 33.27, 72.58, 40.91, 39.38, 41.90,
     +                            51.60, 48.53, 23.75, 54.98, 55.52,
     +                            56.53, 61.45, 61.86, 60.81, 44.58,
     +                            33.77, 41.85, 75.22, 62.17, 47.71/
C
C TOTAL ACCESSIBLE SURFACE AREA IN FULLY EXTENDED open CHAIN FORM
C
      data (ATOTAL(i,2),i=1,nrt)/109.19,243.03,144.93,141.31,136.20,
     +                           180.60,174.46, 81.20,183.87,177.22,
     +                           180.35,204.10,197.05,199.73,140.92,
     +                           118.24,142.04,248.09,212.29,153.33/
C
C TOTAL NON POLAR AREA OF SIDECHAIN IN EXTENDED FORM (INCLUDING CA)
C
      data (ANPOLSIDE(I,1),I=1,NRT)/ 23.16, 25.92, 14.56, 15.95, 31.96,
     +                               17.00, 19.75, 10.84, 30.43, 45.71,
     +                               46.61, 38.72, 51.94, 51.54, 39.28,
     +                               15.76, 24.92, 58.15, 42.35, 37.89/
C
C TOTAL NON-POLAR ACCESSIBLE SURFACE AREA OF SIDECHAIN IN EXTENDED
C FORM (INCLUDING CA)
C
      data (anpolside(i,2),i=1,nrt)/ 70.81, 79.70, 45.01, 49.50, 98.3,
     +                               52.53, 61.20, 33.14, 96.38,139.76,
     +                              142.52,118.39,159.23,164.49,120.12,
     +                               48.19, 76.21,185.88,134.89,115.87/
C
C TOTAL POLAR AREA OF SIDE CHAIN IN EXTENDED CONFORMATION
C
      data (APOLSIDE(I,1),I=1,NRT)/   0.0, 36.74, 16.38, 13.44,  0.0,
     +                              24.68, 18.86,   0.0, 15.18,  0.0,
     +                                0.0, 12.81,   0.0,   0.0,  0.0,
     +                               7.86,  7.01,  6.93, 10.56,  0.0/
C
C TOTAL POLAR ACCESSIBLE SURFACE AREA OF SIDECHAIN IN EXTENDED CONFORM
C
      DATA (apolside(i,2),i=1,nrt)/ 0.00,125.54, 61.88, 53.78,   0.0,
     +                             90.25, 75.42,   0.0, 51.88,   0.0,
     +                               0.0, 47.88,   0.0,  0.00,  0.00,
     +                             31.44, 28.04, 23.69, 42.22,  0.00/
C
C TOTAL SIDE CHAIN AREA (INCLUDING CA)
C
      data (ASIDE(I,1),I=1,NRT)/23.16,62.66,30.93,29.40,31.96,
     +                          41.67,38.61,10.84,45.61,45.71,
     +                          46.61,51.53,51.94,51.54,39.28,
     +                          23.62,31.93,65.08,52.91,37.89/
C
C TOTAL SIDE CHAIN ACCESSIBLE SURFACE AREA (INCLUDING CA)
C
      DATA (aside(i,2),i=1,nrt)/ 70.81,205.23,106.90,103.28, 98.3,
     +                          142.78,136.63, 33.14,148.25,139.76,
     +                          142.52,166.28,159.23,164.49,120.12,
     +                           79.63,104.24,209.56,177.11,115.87/
C
C MAIN CHAIN AREAS IN EXTENDED FORM (NOT INCLUDING CA)
C
      data (AMAIN(I,1),I=1,NRT)/10.11, 9.92, 9.98, 9.98, 9.94,
     +                           9.92, 9.92,12.91, 9.37, 9.81,
     +                           9.92, 9.92, 9.92, 9.28, 5.29,
     +                          10.15, 9.92,10.14, 9.26, 9.82/
C
C MAIN CHAIN ACCESSIBLE SURFACE AREA IN EXTENDED FORM (NOT INCLUDING CA)
C
      DATA (amain(i,2),i=1,nrt)/38.38,37.83,38.03,38.03,37.90,
     +                          37.83,37.83,48.06,35.62,37.45,
     +                          37.83,37.83,37.83,35.24,20.80,
     +                          38.61,37.80,38.53,35.18,37.46/
C
C    +----------------------------------------------------------------------
C
      if (ACCTYP .eq. 'CONTACT') then
        IA=1
      else
        IA=2
      end if
C
      NACCESS=0
C
      do IRES=1, NRES
C
C is the current residue one of the 20 standard residue types (ASX,
C and GLX are used as ASN, and GLN respectively): UNK used to be treated
C the same as GLY, not now
C
        ITYP=ICODE(RESNAM(IRES))
C
C if yes, process and write a line to output file, otherwise ignore it
C
        if (ITYP .eq. 21) then
         if (OUTHET) then
c ------- for HETATM only this:
          NACCESS=NACCESS+1
          NUMACC(NACCESS)=RESNUM(IRES)
          RESACC(NACCESS)=RESNAM(IRES)
          do K=1, 10
            RACCSS(K,NACCESS)=0.0
          end do
         end if
        else
          NACCESS=NACCESS+1
          NUMACC(NACCESS)=RESNUM(IRES)
          RESACC(NACCESS)=RESNAM(IRES)
          do K=1, 10
            RACCSS(K,NACCESS)=0.0
          end do
C
C add an area of each atom to the appropriate bin:
C
          do IAT=IATMR1(IRES), IATMR2(IRES)
            if (IMNCH(ATMNAM(IAT)) .gt. 0) then
C
C it is a main chain atom:
C
                RACCSS(9,NACCESS)=RACCSS(9,NACCESS)+ACCSS(IAT)
                RACCSS(1,NACCESS)=RACCSS(1,NACCESS)+ACCSS(IAT)
            else
              if (IPOLSDCH(ATMNAM(IAT)) .gt. 0) then
C
C it is a polar side chain atom:
C
                RACCSS(7,NACCESS)=RACCSS(7,NACCESS)+ACCSS(IAT)
                RACCSS(5,NACCESS)=RACCSS(5,NACCESS)+ACCSS(IAT)
                RACCSS(1,NACCESS)=RACCSS(1,NACCESS)+ACCSS(IAT)
              else
C
C it must be a non-polar side chain atom then:
C
                RACCSS(7,NACCESS)=RACCSS(7,NACCESS)+ACCSS(IAT)
                RACCSS(3,NACCESS)=RACCSS(3,NACCESS)+ACCSS(IAT)
                RACCSS(1,NACCESS)=RACCSS(1,NACCESS)+ACCSS(IAT)
              end if
            end if
          end do
C
C get the percentages of contact areas for this residue:
C
          RACCSS( 2,NACCESS)=PERCENTAGE(RACCSS(1,NACCESS),
     +                       ATOTAL(ITYP,IA))
          RACCSS( 4,NACCESS)=PERCENTAGE(RACCSS(3,NACCESS),
     +                       ANPOLSIDE(ITYP,IA))
          RACCSS( 6,NACCESS)=PERCENTAGE(RACCSS(5,NACCESS),
     +                       APOLSIDE(ITYP,IA))
          RACCSS( 8,NACCESS)=PERCENTAGE(RACCSS(7,NACCESS),
     +                       ASIDE(ITYP,IA))
          RACCSS(10,NACCESS)=PERCENTAGE(RACCSS(9,NACCESS),
     +                       AMAIN(ITYP,IA))
        end if
      end do
C
      return
      end
