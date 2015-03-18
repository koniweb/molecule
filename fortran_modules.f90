MODULE fortran_modules

  IMPLICIT NONE

CONTAINS 
  
  subroutine define_bonds(cutoff,cutmin,natoms,ndim,coords,arange,periodicity,bondmax,vec,bonding)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)                           :: cutoff
    DOUBLE PRECISION, INTENT(IN)                           :: cutmin
    DOUBLE PRECISION, DIMENSION(0:natoms,ndim), INTENT(IN) :: coords
    INTEGER, DIMENSION(0:1)                                :: arange
    LOGICAL,INTENT(IN)                                     :: periodicity
    DOUBLE PRECISION, DIMENSION(ndim,ndim), INTENT(IN)     :: vec
    INTEGER(KIND=8), DIMENSION(bondmax,2+ndim),INTENT(INOUT)    :: bonding
    ! local data
    INTEGER :: natoms,ndim,bondmax
    ! counter
    DOUBLE PRECISION :: distsq,cutoffsq,cutminsq
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

    ! squared
    cutoffsq=cutoff**2
    cutminsq=cutmin**2
    
    ! print vec
    !write (*,*) vec(1,:)

    ! loop over periodic boxes
    DO i=0,natoms
       IF (mod(i,100) .eq. 0) write(0,*) "...bonding for atom ",i, " of ",natoms," calculated"
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
                      STOP
                   END IF
                      
                END DO
             END DO
          END DO

       END DO
    END DO
    
    
  END subroutine define_bonds
  
END MODULE fortran_modules
