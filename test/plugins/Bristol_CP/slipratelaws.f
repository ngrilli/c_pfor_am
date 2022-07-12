c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c
c Aug. 12th, 2021 - 1st working version
c
	module slipratelaws
      implicit none
      contains



      subroutine sliprate(tau,sstate,temp,sXdist,ph_no,gdot,dgdot_dtau)
      use globalvars, only: modelno, sliprate_param, sliphard_param,
     &creepno, creep_param, numstvar, G, KB
      implicit none
      real(8) tau, tauc, x, temp, sXdist, gdot, dgdot_dtau
      real(8)	sstate(numstvar), gdotc, dgdotc_dtau
      integer ph_no
c     model constants for slip/creep
      real(8) gdot0, m, crit
c     model constants for slip with backstress (Code Aster model)
      real(8) K, n, pl
c     model constatns for GND using slip gradients
      real(8) alpha, b, tauc0
c     model constatns for slip/creep with backstress (Power Law)
      real(8) C
c     Irradiation hardening law constants
      real(8) lambda
c     Creep model-1 constants
      real(8) gdot0c, mc
c     Creep model-2 constants
      real(8) A, D
c     Creep model-3 constants
      real(8) CD, V
c     Creep model-4 constants
      real(8) Qc

c     Slip/Creep  with no kinematic hardening (x = 0), s it has no effect
      if (modelno.eq.1d+0) then

c          write(6,*) 'sliprate'
c          write(6,*) sstate

      tauc = sstate(1)


      gdot0 = sliprate_param(ph_no,1)

      m = sliprate_param(ph_no,2)


      gdot = gdot0 * ((dabs(tau)/tauc)**(1.0d+0/m))
     &*dsign(1.0d+0,tau)

      dgdot_dtau = gdot0/m * ((dabs(tau)/tauc)
     &**((1.0d+0/m) - 1.0d+0)) / tauc




c     Slip and Creep with kinematic hardening
      elseif (modelno.eq.2d+0) then





      tauc = sstate(1)
      x = sstate(2)


      gdot0 = sliprate_param(ph_no,1)

      m = sliprate_param(ph_no,2)




      gdot=gdot0*((dabs(tau-x)/tauc)**(1.0d+0/m))*
     &dsign(1.0d+0,tau-x)

      dgdot_dtau=gdot0/m*((dabs(tau-x)/tauc)
     &**(1.0d+0/m-1.0d+0))/tauc



c     Slip with  kinematic hardening (Code Aster)
      elseif (modelno.eq.3d+0) then

c          write(6,*) 'sliprate'
c          write(6,*) sstate

      tauc = sstate(1)
      x = sstate(2)

      K = sliprate_param(ph_no,1)

      n = sliprate_param(ph_no,2)

      pl = (dabs(tau-x)-tauc)/K


      if (pl.le.0.0) then

        gdot = 0.0d+0

        dgdot_dtau = 0.0d+0

      else

        gdot=(pl**n)*dsign(1.0d+0,tau-x)

        dgdot_dtau = pl**(n-1.0d+0)/K/n


      endif



c     Slip/Creep considering GNDs with slip gradients
      elseif (modelno.eq.4d+0) then

c          write(6,*) 'sliprate'
c          write(6,*) sstate

c         x is not a state but it will be derived from states
c         calculate backstress here!
      x = sstate(4)


			gdot0 = sliprate_param(ph_no,1)

			m = sliprate_param(ph_no,2)


c         Geometric factor for slip resistance calculation
			alpha = sliphard_param(ph_no,4)



c         Burgers vector for slip resistance calculation
			b = sliphard_param(ph_no,3)

c         Initial slip resistance
			tauc0 = sliphard_param(ph_no,1)


c         Calculate slip resistance
      tauc = tauc0 + alpha * G(ph_no) * b *
     & dsqrt(sstate(1) + sstate(2) +sstate(3))

c          write(6,*) 'tauc', tauc



c          write(6,*) 'x', x


c          if (dabs(tau).gt.dabs(x)) then
      gdot=gdot0*((dabs(tau-x)/tauc)**(1.0d+0/m))*
     &dsign(1.0d+0,tau-x)
c
      dgdot_dtau=gdot0/m*((dabs(tau-x)/tauc)**
     &((1.0d+0/m)-1.0d+0))/tauc
c          else
c              gdot= 0.0
c              dgdot_dtau= 0.0
c          endif


			!if (isnan(gdot)) then
			!    write(6,*) 'tau', tau
			!    write(6,*) 'tauc', tauc
			!    write(6,*) 'x', x
			!    write(6,*) 'sstate(1)', sstate(1)
			!    write(6,*) 'sstate(2)', sstate(2)
			!    write(6,*) 'sstate(3)', sstate(3)
			!endif

c          write(6,*) 'gdot', gdot





c     Slip/Creep  with kinematic hardening
      elseif (modelno.eq.5d+0) then

c          write(6,*) 'sliprate'
c          write(6,*) sstate


			x = sstate(2)


c         Geometric factor for slip resistance calculation
			alpha = sliphard_param(ph_no,6)



c         Burgers vector for slip resistance calculation
			b = sliphard_param(ph_no,3)

c         Initial slip resistance
			tauc0 = sliphard_param(ph_no,1)


c         Scaling factor for size dependent hardening
			C = sliphard_param(ph_no,9)


			gdot0 = sliprate_param(ph_no,1)

			m = sliprate_param(ph_no,2)




c         Calculate slip resistance
			tauc = tauc0 + alpha*G(ph_no)*b*(dsqrt(sstate(1))+C/sXdist)


      gdot=gdot0*((dabs(tau-x)/tauc)**(1.0d+0/m))*
     &dsign(1.0d+0,tau-x)

      dgdot_dtau=gdot0/m*((dabs(tau-x)/tauc)
     &**((1.0d+0/m)-1.0d+0))/tauc




c     Irradiation hardening model
      elseif (modelno.eq.6d+0) then

c          write(6,*) 'sliprate'
c          write(6,*) sstate

c         x is not a state but it will be derived from states
c         calculate backstress here!

      x = sstate(4)
c         Initial slip rate
      gdot0 = sliprate_param(ph_no,1)
c         Power law exponent
      m = sliprate_param(ph_no,2)
c         Statistical coefficient
      lambda = sliprate_param(ph_no,3)


c         Geometric factor for slip resistance calculation
      alpha = sliphard_param(ph_no,4)


c         Burgers vector for slip resistance calculation
      b = sliphard_param(ph_no,3)

c         Initial slip resistance
      tauc0 = sliphard_param(ph_no,1)


c         Calculate slip resistance
      tauc = tauc0+lambda*G(ph_no)*b*
     & dsqrt(sstate(1)+sstate(2)+sstate(3)+sstate(5))




      gdot=gdot0*((dabs(tau-x)/tauc)**(1.0d+0/m))*
     & dsign(1.0d+0,tau-x)
c
      dgdot_dtau=gdot0/m*((dabs(tau-x)/tauc)**
     & ((1.0d+0/m)-1.0d+0))/tauc






c     Other slip rate laws come here!!!


      endif




c     CREEP LAWS!!!
      gdotc = 0.0d+0
      dgdotc_dtau = 0.0d+0
c     Creep law of Dylan and Abdullah
      if (creepno.eq.1d+0) then


        gdot0c = creep_param(ph_no,1)

        mc = creep_param(ph_no,2)



        gdotc = gdot0c*((dabs(tau-x)/tauc)**(1.0d+0/mc))
     & *dsign(1.0d+0,tau-x)


        dgdotc_dtau = gdot0c/mc*((dabs(tau-x)/tauc)
     & **(1.0d+0/mc-1.0d+0))/tauc



c     Creep law of Takeuchi and Argon
      elseif (creepno.eq.2d+0) then


        A = creep_param(ph_no,1)

        b = creep_param(ph_no,2)

        D = creep_param(ph_no,3)

        n = creep_param(ph_no,4)


        gdotc = A*G(ph_no)*b*D/KB/temp*((dabs(tau-x)/G(ph_no))**n)
     & *dsign(1.0d+0,tau-x)*1.0d+3


        dgdotc_dtau = A*b*D/KB/temp/n*((dabs(tau-x)/G(ph_no))
     & **(n-1.0d+0))*dsign(1.0d+0,tau-x)*1.0d+3





c     Sinh law
      elseif (creepno.eq.3d+0) then


        CD = creep_param(ph_no,1)

        n = creep_param(ph_no,2)

        V = creep_param(ph_no,3)


        gdotc = CD*(dsinh(dabs(tau-x)*V/KB/temp))**n
     & *dsign(1.0d+0,tau-x)

        dgdotc_dtau = CD*n*(dsinh(dabs(tau-x)*V/KB/temp))
     & **(n-1.0d+0)*dcosh(dabs(tau-x)*V/KB/temp)*V/KB/temp



c    Exponential law
      elseif (creepno.eq.4d+0) then



        B = creep_param(ph_no,1)

        Qc = creep_param(ph_no,2)

        V = creep_param(ph_no,3)


        gdotc = B*dexp(-Qc/KB/temp)*
     & dexp(dabs(tau-x)*V/KB/temp)*dsign(1.0d+0,tau-x)


        dgdotc_dtau = B*dexp(-Qc/KB/temp)*
     & dexp(dabs(tau-x)*V/KB/temp)*V/KB/temp






      endif





c     Add creep rates to slip rate
      gdot = gdot + gdotc

      dgdot_dtau = dgdot_dtau + dgdotc_dtau






      return
      end subroutine sliprate









      end module slipratelaws
