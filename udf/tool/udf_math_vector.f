        subroutine vector_cross_product(a,b,n)
        !============================================================
        ! calculate the cross product of two vectors
        !------------------------------------------------------------
        ! input: 
        ! a(3),b(3)      - input vectors
        ! output:
        ! n(3)           - output of cross vectors
        !------------------------------------------------------------
        implicit none
        double precision a(3)
        double precision b(3)
        double precision n(3)

        n(1)=a(2)*b(3)-a(3)*b(2)
        n(2)=a(3)*b(1)-a(1)*b(3)
        n(3)=a(1)*b(2)-a(2)*b(1)

        end subroutine vector_cross_product

        function vector_norm(v,n) result(val)
        !============================================================
        ! calculate the norm of the vectors
        !------------------------------------------------------------
        ! input: 
        ! v(n)      - input vectors
        ! n         - dim of the vector
        ! output:
        ! val       - norm of the vectors
        !------------------------------------------------------------
        implicit none
        integer n
        integer i
        double precision tol,val
        double precision v(n)

        tol=0.
        do i=1,n
            tol=tol+v(i)*v(i)
        enddo
        val=sqrt(tol)

        end function vector_norm

        subroutine vector_normalize(v,n,v11)
        !============================================================
        ! calculate the normalized vector
        !------------------------------------------------------------
        ! input: 
        ! v(n)      - input vectors
        ! n         - dim of the vector
        ! output:
        ! v11(n)    - normalized input vectors
        !------------------------------------------------------------
        implicit none
        integer i
        integer n
        double precision v(n),v11(n)
        double precision val
        double precision vector_norm

        val=vector_norm(v,n)
        do i=1,n
            v11(i)=v(i)/val
        enddo

        end subroutine vector_normalize

        subroutine vector_dyadic_product(a,b,n,mat)
        !============================================================
        ! calculate the dyadic product of vector
        !------------------------------------------------------------
        ! input: 
        ! a(n),b(n)      - input vectors
        ! n              - dim of the vector
        ! output:
        ! mat(n,n)       - dyadic product
        !------------------------------------------------------------
        implicit none
        integer n
        double precision a(n),b(n)
        double precision mat(n,n)
        integer i,j

        do i=1,n
            do j=1,n
                mat(i,j)=a(i)*b(j)
            enddo
        enddo

        end subroutine vector_dyadic_product