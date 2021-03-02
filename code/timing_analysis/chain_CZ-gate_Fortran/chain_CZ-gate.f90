MODULE utilities
	IMPLICIT NONE
    CONTAINS

	FUNCTION tensor_product(a,b) result(c)

        COMPLEX*16 :: a(:,:), b(:,:)
        COMPLEX*16, DIMENSION (:,:), ALLOCATABLE :: c
        INTEGER :: NA, NB, ii, jj

        NA = size(A,1)
        NB = size(B,1)
        ALLOCATE( C(NA*NB,NA*NB) )

        DO ii=1, NA
          DO jj=1, NA
            C((ii-1)*NB+1:ii*NB ,(jj-1)*NB+1:jj*NB ) = A(ii,jj)*B
          END DO
        END DO

        RETURN
    END FUNCTION tensor_product

	FUNCTION tensor_product_vect(A,B) result(C)

        COMPLEX*16, DIMENSION (:), ALLOCATABLE :: A, B, C
        INTEGER :: NA, NB, ii, jj

        NA = size(A,1)
        NB = size(B,1)
        ALLOCATE( C(NA*NB) )

        DO ii=1, NA
            C((ii-1)*NB+1:ii*NB) = A(ii)*B
        END DO

        RETURN
    END FUNCTION tensor_product_vect

	FUNCTION outer_product(psi1,psi2) RESULT(rho)

        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi1,psi2
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: rho
		INTEGER :: N, M

        N = SIZE(psi1,1)
        M = SIZE(psi2,1)
        ALLOCATE(rho(N,M))

        rho = MATMUL( RESHAPE(psi1,(/N,1/)), RESHAPE(CONJG(psi2),(/1,M/)) )

        RETURN
    END FUNCTION outer_product

	SUBROUTINE eig_comp(mat, eigenvalues)

        COMPLEX*16, DIMENSION(:,:) :: mat
        REAL*8, DIMENSION(SIZE(mat,1)) :: eigenvalues
        REAL*8, DIMENSION(:), ALLOCATABLE :: rwork
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: work
        CHARACTER(1) :: jobz, uplo
        INTEGER :: dims, lda, lwork, info

        dims = SIZE(mat,1)
        lda = dims
        lwork = 2 * dims - 1
        jobz = "V"
        uplo = "U"
        ALLOCATE(rwork(3 * dims - 2))
        ALLOCATE(work(2* dims - 1))

        CALL ZHEEV(jobz, uplo, dims, mat, lda, eigenvalues, work, lwork, rwork, info)
    END SUBROUTINE eig_comp

	SUBROUTINE inv(mat)
		COMPLEX*16, DIMENSION(:,:) :: mat
		COMPLEX*16, DIMENSION(:), ALLOCATABLE :: work
		INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
		INTEGER :: n, info

		ALLOCATE(work(SIZE(mat,1)))
		ALLOCATE(ipiv(SIZE(mat,1)))

		n = size(mat,1)

		call ZGETRF(n, n, mat, n, ipiv, info)
		call ZGETRI(n, mat, n, ipiv, work, n, info)

	END SUBROUTINE inv

	FUNCTION unit_evol(H, psi, tau) RESULT(psi_fin)

		COMPLEX*16, DIMENSION(:,:) :: H
		COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H_copy, exp_diag, U
		COMPLEX*16, DIMENSION(:) :: psi
		COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi_fin
		REAL*8, DIMENSION(:), ALLOCATABLE :: eigenvalues
		REAL*8 :: tau
		INTEGER :: ii

		ALLOCATE(eigenvalues(SIZE(H,1)))

       	CALL eig_comp(H,eigenvalues)

		H_copy = H

        CALL inv(H)

		ALLOCATE(exp_diag(SIZE(eigenvalues,1),SIZE(eigenvalues,1)))

		exp_diag = 0.0

		DO ii = 1, SIZE(eigenvalues,1)
        	exp_diag(ii,ii) = EXP(CMPLX(0.0,-1.0)*eigenvalues(ii)*tau) !
		END DO

        U = MATMUL(H_copy,MATMUL(exp_diag,H))
		psi_fin = MATMUL(U,psi)

    END FUNCTION unit_evol

	FUNCTION basis(dim,state) RESULT(psi)

		INTEGER dim, state
		COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi

		ALLOCATE(psi(dim))

		psi = 0.0
		psi(state) = 1.0

		RETURN
	END FUNCTION basis

	FUNCTION hamiltonian(Omega,Delta) RESULT(H)

		COMPLEX*16 :: Omega, Delta
		COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi00, psi01, psi0r, psi10, psi11, psi1r, psir0, psir1, psirr
		COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H0, H01, H10, H2, H, H01_bis, H10_bis, H_a, H_b, H_c, H_d

		psi00 = tensor_product_vect(basis(3,1),basis(3,1))
		psi01 = tensor_product_vect(basis(3,1),basis(3,2))
		psi0r = tensor_product_vect(basis(3,1),basis(3,3))
		psi10 = tensor_product_vect(basis(3,2),basis(3,1))
		psi11 = tensor_product_vect(basis(3,2),basis(3,2))
		psi1r = tensor_product_vect(basis(3,2),basis(3,3))
		psir0 = tensor_product_vect(basis(3,3),basis(3,1))
		psir1 = tensor_product_vect(basis(3,3),basis(3,2))
		psirr = tensor_product_vect(basis(3,3),basis(3,2))

		H0  = 0 * outer_product(psi00,psi00)

		H01 = 1.0/2.0 * ( Omega * outer_product(psi01,psi0r) + CONJG(Omega) * outer_product(psi0r,psi01) )
		H01_bis = -1.0 * Delta * outer_product(psi0r,psi0r)

		H10 = 10./2.0 * ( Omega * outer_product(psi10,psir0) + CONJG(Omega) * outer_product(psir0,psi10) )
		H10_bis = -1.0 * Delta * outer_product(psir0,psir0)

		H_a = Omega * ( outer_product(psi11,psir1) + outer_product(psi11,psi1r) )
		H_b = CONJG(Omega) * ( outer_product(psir1,psi11) + outer_product(psi1r,psi11) )
		H_c = outer_product(psir1,psir1) + outer_product(psir1,psi1r)
		H_d = outer_product(psi1r,psir1) + outer_product(psi1r,psi1r)
		H2  = 1.0/2.0 * ( H_a + H_b ) - Delta/2.0 * ( H_c  + H_d )

		H = H0 + H01 + H01_bis + H10 + H10_bis + H2

	END FUNCTION hamiltonian

	FUNCTION chain_init(N,state_first,state_last) RESULT(psi)

		INTEGER :: N, state_first, state_last, ii
		COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi

		psi = basis(3,state_first)

		DO ii = 1, N-2
			psi = tensor_product_vect(psi,basis(3,1))
		END DO

		psi = tensor_product_vect(psi,basis(3,state_last))

	END FUNCTION chain_init

	FUNCTION chain_CZ_gate_time(psi,Omega,Delta,tau) RESULT(time)

		COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi, psi_fin
		COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: id, Hp1, Hp2, H1, H2, H11, H22
		COMPLEX*16 :: Omega, Delta, exp_xi
		REAL*8 :: tau
		REAL*8 :: start, finish, time
		INTEGER :: N, ii, jj

		CALL CPU_TIME(start)

		N = INT(LOG10(REAL(SIZE(psi,1)))/LOG10(3.0))
		exp_xi = CMPLX(-0.7242632186311062,0.6895235965056673)

		ALLOCATE(id(3,3))
		id = 0.0

		DO ii = 1, 3
        	id(ii,ii) = 1.0
		END DO

		Hp1 = hamiltonian(Omega,Delta)
    	Hp2 = hamiltonian(Omega * exp_xi, Delta)

		DO ii = 1, N-1

			IF (ii == 1) THEN
				H1 = Hp1
				H2 = Hp2
			ELSE
				H1 = id
				H2 = id
			END IF

			DO jj = 1, N-2

				IF (jj == ii-1) THEN
					H1 = tensor_product(H1,Hp1)
					H2 = tensor_product(H2,Hp2)
				ELSE
					H1 = tensor_product(H1,id)
					H2 = tensor_product(H2,id)
				END IF

			END DO

			psi = unit_evol(H1, psi, tau)
           	psi = unit_evol(H2, psi, tau)

		END DO

		psi_fin = psi
		CALL CPU_TIME(finish)
		time = finish - start

	END FUNCTION chain_CZ_gate_time

	FUNCTION mean(a) RESULT(m)

		REAL*8, DIMENSION(:), ALLOCATABLE :: a
		REAL*8 :: m

		m = SUM(a)/SIZE(a,1)

	END FUNCTION mean

	FUNCTION std(a) RESULT(s)

		REAL*8, DIMENSION(:), ALLOCATABLE :: a
		REAL*8 :: s, var, m
		INTEGER :: ii

		var = 0.0
		m = SUM(a)/SIZE(a,1)

		DO ii = 1, SIZE(a,1)
			var = var + (a(ii) - m)**2
		END DO

		var = var / (SIZE(a,1) - 1)
		s = SQRT(var)

	END FUNCTION std

END MODULE utilities

PROGRAM chain_CZ_gate

	USE utilities

   	IMPLICIT NONE

	REAL*8, DIMENSION(:), ALLOCATABLE:: times
	COMPLEX*16, DIMENSION(:), ALLOCATABLE :: psi
	INTEGER, DIMENSION(:), ALLOCATABLE :: N
	COMPLEX*16 :: Omega, Delta
	REAL*8 :: tau
	INTEGER :: ii, jj, niter

	Omega = 1.0
	Delta = 0.377371
	tau = 4.29268

	niter = 10
	ALLOCATE(times(niter))

	N = (/ 2, 3, 4, 5, 6/)

	OPEN(1, file = 'fortran.txt', status = "replace")
	DO ii = 1, SIZE(N,1)
		DO jj = 1, niter
			PRINT *, "N: ", N(ii), "iter: ", jj
			psi = chain_init(N(ii),2,2)
			times(jj) = chain_CZ_gate_time(psi,Omega,Delta,tau)
		END DO
		WRITE (1,*) N(ii), mean(times), std(times)/SQRT(REAL(N(ii)))
	END DO
	CLOSE(1)

END PROGRAM chain_CZ_gate
