! This routine calls the DAVIDSON package to find the ground state of a large sparse Hermitian matrix
! without constructing the full Hamiltonian matrix, but giving only the instructions
! to calculate:  |PsiOut> = Ham * |PsiIn>    [see subroutine AMUL]
! see: http://www.staff.science.uu.nl/~sleij101/JD_software/jd.html
     CALL DAVIDSON(NumTot, Kmax, Energy, GS1)
  ENDIF

  CALL Magnet(GS1, MagnetX, MagnetZ)

! Prints out the lowest three energy levels
  PRINT *,'First three energy levels:   ', Energy(1:3)

! Prints out the magnetizations along X and Z
  IF (PBC) THEN
     PRINT *,'GS magnetization along X:   ', MagnetX(1)
     PRINT *,'GS magnetization along Z:   ', MagnetZ(1)
  ELSE
     WRITE (6,111)  MagnetX
     WRITE (6,112)  MagnetZ
111  FORMAT ('GS magnetization profile along X:   ', 100(ES15.8,3x))
112  FORMAT ('GS magnetization profile along Z:   ', 100(ES15.8,3x))
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Magnet(Psi, MagX, MagZ)
    USE SystemParameters

    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(ell)    :: MagX, MagZ
    DOUBLE COMPLEX,   DIMENSION(NumTot) :: Psi
    DOUBLE COMPLEX   :: Mx_sum
    DOUBLE PRECISION :: MagnetZ
    INTEGER          :: iSite, ii, Exc, Mag_ii

    MagX = 0.d0
    MagZ = 0.d0
    DO iSite = 1,ell
       Mx_sum = 0.d0
       DO ii = 1,NumTot
          IF (iSpin(ii,iSite) == 1)   MagZ(iSite) = MagZ(iSite) + Abs(Psi(ii))**2
          IF (iSpin(ii,iSite) == 0)   MagZ(iSite) = MagZ(iSite) - Abs(Psi(ii))**2

          IF (iSpin(ii,iSite) == 1)   Exc = ii -2**(iSite-1)
          IF (iSpin(ii,iSite) == 0)   Exc = ii +2**(iSite-1)
          Mx_sum = Mx_sum + dConjg(Psi(ii)) * Psi(Exc)
       ENDDO
       IF (Abs(aImag(Mx_sum)) > 1.d-10)  STOP 'Non real magnetization'
       MagX(iSite) = dReal(Mx_sum)
    ENDDO

    DO ii = 1,NumTot
       Mag_ii = 0
       DO iSite = 1,ell
          IF (iSpin(ii,iSite) == 1)   Mag_ii = Mag_ii + 1
          IF (iSpin(ii,iSite) == 0)   Mag_ii = Mag_ii - 1
       ENDDO
       MagnetZ = MagnetZ + Abs(1.d0*Mag_ii) * Abs(Psi(ii))**2
    ENDDO
    PRINT *,'Symmetry-broken magnetization along Z:   ', MagnetZ/ell

  END SUBROUTINE Magnet

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Davidson(n, Kmax, Eigenvalue, EigenVector)
    USE SystemParameters

    IMPLICIT NONE
    LOGICAL           UseGuess,Wanted
    INTEGER           kmax,jmax,jmin,maxstep,method,m,l,maxnmv,order,testspace,j,lwork,istate,ii,n,kmaxuser
    DOUBLE PRECISION  tol,lock,targetEn,Norm,emin,etemp
    DOUBLE PRECISION, DIMENSION(Kmax)             :: EigenValue
    DOUBLE COMPLEX,   DIMENSION(n)                :: EigenVector
    DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: alpha,beta,tmp,residu
    DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: eivec,zwork


    !!  INIZIALIZATION OF PARAMETERS  !!
    Useguess = .false.
    KMaxUser = KMax
    targetEn = -5.d0*ell
    tol = 1.d-9    ! Tolerance of the eigensolutions: $\Vert \beta H_{SB} x - \alpha x \vert$
    maxnmv = 100    ! Maximum number of matvecs in cgstab or gmres (very problem dependent; typically choose in [5-100])
    wanted = .true. ! If wanted=.true. then computes the converged eigenvectors
    order = -1      ! Selection criterion for Ritz values:  0 (nearest to target);  -1 (smallest real part)
    IF(order == 0)  testspace = 3 ! put 3 if a reasonable value for target is known, else take 2
    IF(order /= 0)  testspace = 2

    IF (3*KmaxUser <= 20) jmax=20          ! Maximum size of the search space:
    IF (3*KmaxUser >  20) jmax=3*KmaxUser
    jmin=2*KmaxUser                        ! Minimum size of the search space
    maxstep = 1000                         ! Maximum number of Jacobi-Davidson iterations
    lock = 1.d-12                          ! Tracking parameter
    method = 2                             ! Method for the linear equation solver  1: gmres(m)  2: cgstab(l)
    m = 30                                 ! Maximum dimension of searchspace for gmres(m):
    l= 2                                   ! Degree of gmres-polynomial in bi-cgstab(l):
    IF (method == 1) lwork =  4 +  m  + 5*jmax + 3*KmaxUser  ! Size of workspace
    IF (method == 2) lwork = 10 + 6*l + 5*jmax + 3*KmaxUser  !KmaxUser is used since Kmax = 1 gives problems ...!
    !!  END OF INIZIALIZATION  !!

    ALLOCATE (alpha(jmax), beta(jmax), eivec(n,Kmax))
    Alpha=0.d0
    Beta=0.d0
    EiVec=0.d0
    ALLOCATE (tmp(n), residu(n), zwork(n,lwork))
    tmp=0.d0
    residu=0.d0
    zwork=0.d0

    CALL JDQZ(ALPHA, BETA, EIVEC, wanted, n, targetEn, tol, Kmax, jmax, jmin, method, m, l, maxnmv, maxstep, &
              lock, order, testspace, zwork, lwork, UseGuess )

    !     Computes the norms of the residuals:
    DO j = 1, Kmax
       CALL AMUL  ( n, eivec(1,j), residu )
       CALL ZSCAL ( n, beta(j), residu, 1 )
       CALL BMUL  ( n, eivec(1,j), tmp )
       CALL ZAXPY ( n, -alpha(j), tmp, 1, residu, 1 )
    ENDDO
    DEALLOCATE (zwork,tmp,residu)
    Eigenvalue(1:Kmax) = dReal(alpha(1:Kmax)/beta(1:Kmax))

    !     Calculates the smallest eigenvalue (ground state)
    emin=eigenvalue(1)
    istate = 1
    DO ii=2,Kmax
       IF (eigenvalue(ii) < emin) THEN
          emin=eigenvalue(ii)
          istate = ii
       ENDIF
    ENDDO
    IF (istate /= 1) THEN
       etemp=eigenvalue(1)
       eigenvalue(1)=eigenvalue(istate)
       eigenvalue(istate)=etemp
    ENDIF
    DEALLOCATE (alpha,beta)

  !  print *,'istate',istate
  !  Chooses the eigenvector corresponding to the selected eigenvalue
    EigenVector = eivec(:,istate)
    Norm = Sum(dConjg(EigenVector)*(EigenVector))
    EigenVector = EigenVector/(Norm**0.5d0)
    DEALLOCATE (eivec)

  END SUBROUTINE Davidson

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
