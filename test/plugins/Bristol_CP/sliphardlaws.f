c Chris Allen
c Edward Horton
c Eralp Demir
c Hugh Dorward
c Michael Salvini
c
c Aug. 12th, 2021 - 1st working version
c
	module sliphardlaws
      implicit none
      contains
c
c
c
      subroutine sliphard(sstate,gsum,gint,gdot,temp,sstate0,sXdist,
     & ph_no,sstatedot)
      use globalvars, only: modelno,sliphard_param,Rgas,numstvar,
     &numslip, G
      implicit none
			real(8) gsum, gint, gdot, temp, sXdist
      real(8)	sstatedot(numstvar+1), sstate(numstvar), sstate0(numstvar)
      integer ph_no
c     generic variables
      real(8) tauc, x
c     variables for Voce hardening law
      real(8) h0, ss, n, hb
c     variables for isotropic and kinematic hardening with AR thermal softening
      real(8) tauc0, A, Q, d, h, hD
c     variables for isotropic and kinematic hardening (Code Aster model)
      real(8) b, C
c     variable used in GND model with slip gradients
      real(8) yc, K, rhoSSD, rhotot, rhoLOOP
c     variables for Voce hardening law with slip distance
      real(8) alpha, mKM
c     Irradiation hardening parameters
      real(8) AL, rhos
c
c
      sstatedot = 0.0d+0
c
c
c     Three parameter hardening law
      if (modelno.eq.1d+0) then
c
c
c
c
c
      tauc = sstate(1)
c
      h0 = sliphard_param(ph_no,2)
      ss = sliphard_param(ph_no,3)
      n = sliphard_param(ph_no,4)
c
      hb = h0*((1.0d+0-(tauc/ss))**n)

c         self-hardening
      sstatedot(1) = hb*dabs(gdot)
c
c
c          write(6,*) 'inside sliphar'
c          write(6,*) 'selfhrd'
c         Other slip hardening laws come here!!!
c
c          write(6,*) 'sstatedot(1)', sstatedot(1)
c
c           write(6,*) 'sstate(1)', sstate(1)
c
      elseif (modelno.eq.2d+0) then
c
c
c
c
      tauc = sstate(1)
      x = sstate(2)
c
      tauc0 = sliphard_param(ph_no,1)
      h0 = sliphard_param(ph_no,2)
      n = sliphard_param(ph_no,3)
      Q = sliphard_param(ph_no,4)
      d = sliphard_param(ph_no,5)
      A = sliphard_param(ph_no,6)
      h = sliphard_param(ph_no,7)
      hD = sliphard_param(ph_no,8)
c
c
c         self-hardening
      sstatedot(1) = h0*dabs(gdot)*(1.0d+0 +
     & (h0*gsum/tauc0/n))**(n-1.0d+0)
c
c         Thermal softening
      sstatedot(3) = -A * (tauc**d) * dexp(-Q/Rgas/temp)
c
c         kinematic hardening
      sstatedot(2) = h*gdot - hD*x*dabs(gdot)
c
c
      elseif (modelno.eq.3d+0) then
c
c
c
c
      x = sstate(2)
c
      b = sliphard_param(ph_no,2)
c          Q = sliphard_param(3)
      C = sliphard_param(ph_no,4)
      D = sliphard_param(ph_no,5)
c
c
c         self-hardening
c		sstatedot(1) = Q * (1. - dexp(-b*gint))
      sstatedot(1) = (1.0d+0 - dexp(-b*gint))


c         kinematic hardening
      sstatedot(2) = C*gdot - D*x*dabs(gdot)





      elseif (modelno.eq.4d+0) then




          rhoSSD = sstate(1)
c          write(6,*) 'rhoSSD'
c          write(6,*) rhoSSD



c         No rates are defined for GND and backstress state variables

          rhotot = sstate(1) + dabs(sstate(2)) + dabs(sstate(3))

          b = sliphard_param(ph_no,3)
          K = sliphard_param(ph_no,5)
          yc = sliphard_param(ph_no,6)




c        SSD evolution
         sstatedot(1)=dabs(gdot)/b*(K*dsqrt(rhotot) -
     & 2.0d+0*yc*rhoSSD)

c          write(6,*) 'sstatedot(1)', sstatedot(1)


c         GND evolution is calculated at the end of the former time increment

c         GND rate is set to zero for implicit integration (states = 2 & 3)
c         Backstress rate is set to zero for implicit integration (state=4)




c     Three parameter hardening law with slip distance effect (size)
      elseif (modelno.eq.5d+0) then
c
c
          rhoSSD = sstate(1)
          x = sstate(2)
c
c
c
c
c
          b = sliphard_param(ph_no,3)
          K = sliphard_param(ph_no,4)
          yc = sliphard_param(ph_no,5)
          alpha = sliphard_param(ph_no,6)
          h = sliphard_param(ph_no,7)
          hD = sliphard_param(ph_no,8)
          C = sliphard_param(ph_no,9)
c
c
c
c         Slip distance will be sent as an input to this subroutine
c         Slip is not actually the statedot, statedot is used just to avoid new
c         definition of input variables
c         Xdist is nonlocal slip distance, that has a different value for each IP
c
c
c         The RHS term represents the Taylor hardening term:
c         alpha*G*b/L scaled with a constant C
      mKM = K*(dsqrt(rhoSSD) + C/sXdist)
     & - 2.0d+0*yc*rhoSSD
c
c         self-hardening
      sstatedot(1) = mKM * dabs(gdot) / b
c
c
c         kinematic hardening
      sstatedot(2) = h*gdot - hD*x*dabs(gdot)
c
c
c          write(6,*) 'inside sliphar'
c          write(6,*) 'selfhrd'
c         Other slip hardening laws come here!!!
c
c          write(6,*) 'sstatedot(1)', sstatedot(1)
c
c           write(6,*) 'sstate(1)', sstate(1)
c
      elseif (modelno.eq.6d+0) then
c
c         Statistically stored dislocation density
          rhoSSD = sstate(1)
c         Frank dislocation loop density
          rhoLOOP = sstate(5)
c
c
c
c         No rates are defined for GND and backstress state variables
c
          rhotot = sstate(1) + dabs(sstate(2)) + dabs(sstate(3))
c
          b = sliphard_param(ph_no,3)
          K = sliphard_param(ph_no,5)
          yc = sliphard_param(ph_no,6)
c
c Irradiation induced frank dislocation loop annihilation area
          AL = sliphard_param(ph_no,10)
c frank loop dislocation saturation density
          rhos = sliphard_param(ph_no,11)
c
c        SSD evolution
         sstatedot(1)=dabs(gdot)/b*(K*dsqrt(rhotot)-
     & 2.0d+0*yc*rhoSSD)
c        Irradiation induced frank dislocation loop density evolution
c
          sstatedot(5)=AL*(rhos-rhoLOOP)*rhoLOOP*dabs(gdot)
c
c
c         GND evolution is calculated at the end of the former time increment
c         GND rate is set to zero for implicit integration (states = 2 & 3)
c         Backstress rate is set to zero for implicit integration (state=4)

c
      endif
c
c
c
c
c
c
c
c
c
      return
      end subroutine sliphard
c
c
c
c
c
c
c
c
      end module sliphardlaws
