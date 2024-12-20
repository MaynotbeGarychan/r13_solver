        subroutine matrix_inverse(in_mat,n,c)
        !============================================================
        ! Inverse matrix
        !-----------------------------------------------------------
        ! input ...
        ! in_mat(n,n) - array of coefficients for matrix A
        ! n      - dimension
        ! output ...
        ! c(n,n) - inverse matrix of A
        ! comments ...
        !===========================================================
            implicit none 
            integer n
            double precision a(n,n), c(n,n), in_mat(n,n)
            double precision L(n,n), U(n,n), b(n), d(n), x(n)
            double precision coeff
            integer i, j, k
            
            ! step 0: initialization for matrices L and U and b
            ! Fortran 90/95 aloows such operations on matrices
            L=0.0
            U=0.0
            b=0.0
            do i=1,n
                do j=1,n
                    a(i,j)=in_mat(i,j)
                enddo
            enddo
            
            ! step 1: forward elimination
            do k=1, n-1
            do i=k+1,n
                coeff=a(i,k)/a(k,k)
                L(i,k) = coeff
                do j=k+1,n
                    a(i,j) = a(i,j)-coeff*a(k,j)
                end do
            end do
            end do
            
            ! Step 2: prepare L and U matrices 
            ! L matrix is a matrix of the elimination coefficient
            ! + the diagonal elements are 1.0
            do i=1,n
            L(i,i) = 1.0
            end do
            ! U matrix is the upper triangular part of A
            do j=1,n
            do i=1,j
                U(i,j) = a(i,j)
            end do
            end do
            
            ! Step 3: compute columns of the inverse matrix C
            do k=1,n
            b(k)=1.0
            d(1) = b(1)
            ! Step 3a: Solve Ld=b using the forward substitution
            do i=2,n
                d(i)=b(i)
                do j=1,i-1
                d(i) = d(i) - L(i,j)*d(j)
                end do
            end do
            ! Step 3b: Solve Ux=d using the back substitution
            x(n)=d(n)/U(n,n)
            do i = n-1,1,-1
                x(i) = d(i)
                do j=n,i+1,-1
                x(i)=x(i)-U(i,j)*x(j)
                end do
                x(i) = x(i)/u(i,i)
            end do
            ! Step 3c: fill the solutions x(n) into column k of C
            do i=1,n
                c(i,k) = x(i)
            end do
            b(k)=0.0
            end do

        end subroutine matrix_inverse

        subroutine matrix_identity(num,idm)
        !============================================================
        ! Return an identity matrix
        !------------------------------------------------------------
        ! input: 
        ! idm(num,num)  - a given empty matrix
        ! num           - dimension of the matrix
        ! output:
        ! idm(num,num)  - an identity matrix
        !============================================================
            implicit none
            integer i,j,num
            double precision idm(num,num)
    
            do i=1,num
              do j=1,num
                idm(i,j)=0.
              enddo
            enddo
    
            do i=1,num
              idm(i,i)=1.
            enddo
            
        end subroutine matrix_identity

        subroutine matrix_zero(d1,d2,mat)
        !============================================================
        ! Return an zero matrix
        !------------------------------------------------------------
        ! input: 
        ! d1,d2           - dim of the matrix
        ! output:
        ! mat(d1,d2)      - a matrix contain zeros
        !------------------------------------------------------------
        implicit none
        integer d1,d2
        integer i,j
        double precision mat(d1,d2)
        
        do i=1,d1
            do j=1,d2
                mat(i,j)=0.
            enddo
        enddo

        end subroutine matrix_zero

        subroutine matrix_subtract(a,b,d1,d2,c)
        !============================================================
        ! c = a - b , pos-to-pos substract
        !------------------------------------------------------------
        ! input: 
        ! a(d1,d2),b(d1,d2) - input matrix
        ! d1,d2           - dim of the matrix
        ! output:
        ! c(d1,d2)      - output matrix
        !============================================================
        implicit none
        integer d1,d2
        double precision a(d1,d2),b(d1,d2),c(d1,d2)
        integer i,j

        do i=1,d1
            do j=1,d2
                c(i,j)=a(i,j)-b(i,j)
            enddo
        enddo
        
        end subroutine matrix_subtract

        subroutine matrix_add(a,b,d1,d2,c)
        !============================================================
        ! c = a + b , pos-to-pos add
        !------------------------------------------------------------
        ! input: 
        ! a(d1,d2),b(d1,d2) - input matrix
        ! d1,d2           - dim of the matrix
        ! output:
        ! c(d1,d2)      - output matrix
        !============================================================
            implicit none
            integer d1,d2
            double precision a(d1,d2),b(d1,d2),c(d1,d2)
            integer i,j
    
            do i=1,d1
                do j=1,d2
                    c(i,j)=a(i,j)+b(i,j)
                enddo
            enddo

        end subroutine matrix_add

        subroutine matrix_multipy_coeff(a,coeff,d1,d2,b)
        !============================================================
        ! c = a(i,j)*coeff 
        !------------------------------------------------------------
        ! input: 
        ! a(d1,d2)        - input matrix
        ! coeff           - coefficient for multiplication
        ! d1,d2           - dim of the matrix
        ! output:
        ! b(d1,d2)      - output matrix
        !============================================================
            implicit none
            integer d1,d2
            double precision a(d1,d2),b(d1,d2)
            integer i,j
            double precision coeff
    
            do i=1,d1
                do j=1,d2
                    b(i,j)=a(i,j)*coeff
                enddo
            enddo

        end subroutine matrix_multipy_coeff


        subroutine matrix_inner_product(a,b,d1,d2,d3,c)
        !============================================================
        ! a:b = c
        !------------------------------------------------------------
        ! input: 
        ! a(d1,d2),b(d2,d3)  - input matrix
        ! d1,d2,d3         - dim of the matrix
        ! output:
        ! c(d1,d3)      - output matrix
        !============================================================
        implicit none
        integer i,j
        integer d1,d2,d3
        double precision a(d1,d2),b(d2,d3)
        double precision c(d1,d3)
        
        do i=1,d1
            do j=1,d3
                c(i,j)=dot_product(a(i,:),b(:,j))
            enddo
        enddo
        
        end subroutine matrix_inner_product

        subroutine matrix_transpose(a,d1,d2,a_t)
        !============================================================
        ! transpose the matrix
        !------------------------------------------------------------
        ! input: 
        ! a(d1,d2)  - input matrix
        ! d1,d2     - dim of the matrix
        ! output:
        ! a_t(d2,d1)  - matrix after transpose
        !============================================================
        implicit none
        integer d1,d2
        double precision a(d1,d2),a_t(d2,d1)
        integer i,j

        do i=1,d1
            do j=1,d2
                a_t(j,i)=a(i,j)
            enddo
        enddo

        end subroutine matrix_transpose

        subroutine matrix33_det(a,det)
        !============================================================
        ! calculate the determinant of the matrix 3*3
        !------------------------------------------------------------
        ! input: 
        ! a(3,3)  - input matrix
        ! output:
        ! det  - determinant of the matrix
        !============================================================
        implicit none
        double precision a(3,3)
        double precision det

        det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-
     1      a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+
     2      a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

        end subroutine matrix33_det





