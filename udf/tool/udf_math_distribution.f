      function distribution_nml(mag,mean,sd,x) result(val)
        !============================================================
        ! Return a value from a normal distribution (untested)
        !------------------------------------------------------------
        ! input: 
        ! mag   - magnitute of this distribution
        ! mean,sd   - mean and std of this distribution
        ! x         - current x
        ! output:
        ! val   - return value
        !------------------------------------------------------------

        implicit none
        double precision mag,mean,sd
        double precision x
        double precision val
        double precision a,b
        double precision, parameter :: pi=acos(-1.0d0)
        double precision, parameter :: e=2.718281828459045

        a=mag/(sd*sqrt(2*pi))
        b=(((x-mean)/sd)**2)*(-0.5)
        val=a*(e**b)

      end function distribution_nml