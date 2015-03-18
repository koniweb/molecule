MODULE fortran_modules

  IMPLICIT NONE

CONTAINS 
  
  subroutine define_bonds(cutoff,natoms,ndim,coords,arange,periodicity,bondmax,vec,bonding)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)                           :: cutoff
    DOUBLE PRECISION, DIMENSION(0:natoms,ndim), INTENT(IN) :: coords
    INTEGER, DIMENSION(0:1)                                :: arange
    LOGICAL,INTENT(IN)                                     :: periodicity
    DOUBLE PRECISION, DIMENSION(ndim,ndim), INTENT(IN)     :: vec
    INTEGER(KIND=8), DIMENSION(bondmax,2+ndim),INTENT(INOUT)    :: bonding
    ! local data
    INTEGER :: natoms,ndim,bondmax
    ! counter
    DOUBLE PRECISION :: distsq,cutoffsq
    INTEGER :: i,j,nper,x,y,z,dim
    INTEGER :: bondcnt=1

    ! set maximum range
    IF(arange(1)<arange(0)) THEN 
       arange(1)=huge(arange(1))
    END IF
    
    ! loop over atoms
    IF (periodicity .eqv. .TRUE.) THEN
       nper=-1
    ELSE
       nper=0
    END IF

    ! cutoffsq
    cutoffsq=cutoff**2
    
    ! print vec
    !write (*,*) vec(1,:)

    ! loop over periodic boxes
    DO x=-nper,nper
       DO y=-nper,nper
          DO z=-nper,nper
             
             DO i=0,natoms
                DO j=0,natoms

                   ! distance sq calculation
                   IF (i /= j) THEN 
                      distsq=0.0D0
                      DO dim=1,3
                         distsq=(coords(i,dim)-coords(j,dim) +x*vec(dim,1) +y*vec(dim,2) +x*vec(dim,3) )**2
                      END DO
                      
                      ! cutoffsq check
                      IF (distsq < cutoff) THEN
                         bonding(bondcnt,1)=i
                         bonding(bondcnt,2)=j
                         bonding(bondcnt,3)=x
                         bonding(bondcnt,4)=y
                         bonding(bondcnt,5)=z
                         bondcnt=bondcnt+1
                      END IF
                      
                   END IF

                END DO
             END DO
             
          END DO
       END DO
    END DO
    
    
  END subroutine define_bonds
  
END MODULE fortran_modules
