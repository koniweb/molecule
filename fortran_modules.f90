MODULE fortran_modules

  IMPLICIT NONE

CONTAINS 
  
  subroutine define_bonds(cutoff,cutmin,natoms,ndim,coords,alist,ncalc,periodicity,bondmax,vec,bonding,ERROR)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)                           :: cutoff
    DOUBLE PRECISION, INTENT(IN)                           :: cutmin
    DOUBLE PRECISION, DIMENSION(0:natoms,ndim), INTENT(IN) :: coords
    INTEGER, DIMENSION(0:ncalc)                            :: alist
    LOGICAL,INTENT(IN)                                     :: periodicity
    DOUBLE PRECISION, DIMENSION(ndim,ndim), INTENT(IN)     :: vec
    INTEGER(KIND=8), DIMENSION(bondmax,2+ndim),INTENT(INOUT)    :: bonding
    LOGICAL, INTENT(INOUT)                                 :: ERROR
    ! local data
    INTEGER :: natoms,ndim,bondmax,ncalc
    ! counter
    DOUBLE PRECISION :: distsq,cutoffsq,cutminsq
    INTEGER :: i,j,nper,x,y,z,dim
    INTEGER :: bondcnt

    ! set error
    ERROR=.FALSE.

    ! loop over atoms
    IF (periodicity .eqv. .TRUE.) THEN
       nper=-1
    ELSE
       nper=0
    END IF

    ! squared
    cutoffsq=cutoff**2
    cutminsq=cutmin**2
    
    ! print vec
    !write (*,*) vec(1,:)

    ! loop over periodic boxes

    bondcnt=1
    ATOMLOOP: DO i=0,natoms
       IF ((mod(i,100) .eq. 0) .OR. (i==natoms)) THEN
          write(0,*) "...bonding for atom ",i+1, " of ",natoms+1," calculated"
       END IF
       ! only atomrange
       IF ( ANY( alist==i)) THEN
          DO j=0,natoms
          
             DO x=-nper,nper
                DO y=-nper,nper
                   DO z=-nper,nper
                      
                      
                      ! distance sq calculation
                      distsq=0.0D0
                      DO dim=1,3
                         distsq=distsq+(coords(i,dim)-(coords(j,dim) +x*vec(dim,1) +y*vec(dim,2) +z*vec(dim,3)) )**2
                      END DO
                         
                      ! cutoffsq check
                      IF (distsq < cutoffsq .AND. distsq > cutminsq) THEN
                         bonding(bondcnt,1)=i
                         bonding(bondcnt,2)=j
                         bonding(bondcnt,3)=x
                         bonding(bondcnt,4)=y
                         bonding(bondcnt,5)=z
                         bondcnt=bondcnt+1
                         !write(0,*) distsq,cutoff
                         !write (0,*) bondcnt
                      END IF
                      
                      ! CHECK if bondlist to small
                      IF (bondcnt>= bondmax) THEN
                         WRITE(0,*) "ERROR WITH FORTRAN DEFINE BONDS"
                         ERROR=.TRUE.
                         EXIT ATOMLOOP
                      END IF
                         
                   END DO
                END DO
             END DO
          
          END DO
       END IF
    END DO ATOMLOOP
    
    
  END subroutine define_bonds
  
END MODULE fortran_modules
