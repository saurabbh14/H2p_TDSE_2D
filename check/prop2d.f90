module wavepacket_vars
complex*16, allocatable, dimension(:,:):: psi, kprop, psi0
double precision, allocatable, dimension(:):: idensR, idensx, idensR1, idensx1
double precision, allocatable, dimension(:):: idenspR, idenspx

end module

!with vector field correction
!left and right x-absorbers are seperate
!also the field after ionization is turned OFF
subroutine propagation_2D(ewf,chi0, E2, A)

use data_grid
use pot_param
use data_au
use FFTW3
use omp_lib
use wavepacket_vars

 implicit none
! include "/usr/include/fftw3.f"

 integer I, J, K
 integer L, M, LL, void
 integer I_cpm, I_cpmR
 integer*8 planF, planB, planFd, planBd,planBx, planFx 
 integer*8 planFR, planBR, planFx1, planBx1

 character(150) filename
 integer cof_choice

 double precision dt2, time
 double precision E2(Nt), E1, E21, E22
 double precision A(Nt),E2_New(Nt)
 double precision v12, vcol_prop(NR)
 double precision norm
 double precision norm_outx, norm_outR, norm_outx1, norm_outR1
 double precision norm_outxR, norm_outRx, norm_outxR1, norm_outRx1
 double precision norm_addx, norm_addR, norm_add
 double precision norm_addx1, norm_addR1

 double precision evR, evx, epx, epr
 double precision x_out, R_out
 double precision momt_x,momt_R, momt_x1
 
 double precision, intent(in):: ewf(Nx, NR, Nstates),chi0(nr,vstates)
 double precision:: a1e_m(NR), g(Nx)
 complex*16, allocatable, dimension(:,:):: psi_cont
 complex*16, allocatable, dimension(:,:):: psi_out_x_cont, psi_out_x_cont1

 double precision:: cof(Nx), cofR(NR), cofx(Nx)
 double precision:: cofxl(Nx), cofxr(Nx), cpmxl

 complex*16:: corrf, normp, normn
! complex*16, allocatable, dimension(:,:):: psi, kprop, psi0
 complex*16, allocatable, dimension(:,:):: kprop_vel
! double precision, allocatable, dimension(:):: idensR, idensx, idensR1, idensx1
! double precision, allocatable, dimension(:):: idenspR, idenspx
 complex*16, allocatable, dimension(:,:):: psi_out_x, psi_out_xR
 complex*16, allocatable, dimension(:,:):: psi_out_x_Rpx
 complex*16, allocatable, dimension(:,:):: psi_out_x1, psi_out_xR1
 complex*16, allocatable, dimension(:,:):: psi_g1, psi_dum
 double precision, allocatable, dimension(:,:):: JES_en_out_x
 double precision, allocatable, dimension(:,:):: JES_en_out_xR, JES_en_out_Rx
 complex*16, allocatable, dimension(:,:):: psi_out_x_intervals
 double precision, allocatable, dimension(:,:):: rho_psi_out_x_intervals
! complex*16, allocatable, dimension(:,:):: psi_out_xl, psi_out_xr
 complex*16, allocatable, dimension(:):: psi_dumx, psi_dumR
 complex*16, allocatable, dimension(:,:):: psi_out_R, psi_out_Rx
 complex*16, allocatable, dimension(:,:):: psi_out_R1, psi_out_Rx1
 complex*16, allocatable, dimension(:,:):: psi_out_x_momt
 ! complex*16, allocatable, dimension(:,:):: psi_out_x101, psi_out_x151
complex(kind=kind(0.d0)):: psi_ex(NR)

 open(96,file='cof_2d.out',status='unknown')
 open(97,file='cofx_2d.out',status='unknown')
 open(98,file='cofxr_2d.out',status='unknown')
 open(99,file='cofR_2d.out',status='unknown')
 open(100,file='psi0_2d.out',status='unknown')
 open(198,file='gaussian_window.out',status='unknown')
 open(200,file='dens_R.out',status='replace')
 open(201,file='dens_x.out',status='replace')
 open(202,file='Pdens_R.out',status='replace')
 open(203,file='Pdens_x.out',status='replace')
 open(204,file='dens_outx.out',status='replace')
 open(205,file='Pdens_outx.out',status='replace')
 open(206,file='dens_outR.out',status='replace')
 open(207,file='Pdens_outR.out',status='replace')

 open(500,file='pop1_2d.out',status='replace')
 open(501,file='pop2_2d.out',status='replace')
 open(502,file='pl_pop1_2d.out',status='replace')
 open(503,file='neg_pop2_2d.out',status='replace')
 open(504,file='densR1_2d.out',status='replace')
 open(505,file='densR2_2d.out',status='replace')
 open(600,file='ampl1_2d.out',status='replace')
 open(601,file='ampl2_2d.out',status='replace')
 open(800,file='pz_2d.out',status='replace')
 open(800,file='R_2d.out',status='replace')
 open(801,file='x_2d.out',status='replace')
 open(802,file='norm_2d.out',status='replace')
! open(803,file='corr_2d.out',status='replace')
 open(850,file='ionization_yield.out',status='replace')
 open(851,file='dissociation_yield.out',status='replace')

 open(909,file='field_2d.out',status='replace')
! open(906,file='vector_field_2d.out',status='replace')
! open(907,file='field1_2d.out',status='replace')
! open(908,file='field2_2d.out',status='replace') 

 allocate(psi(NR,Nx), idensR(NR), idensx(Nx))
! allocate(psi_cont(NR,Nx),psi_out_x_cont(NR,Nx),psi_out_x_cont1(NR,Nx))
 allocate(kprop(NR,Nx), psi0(nr,nx))
 allocate(idenspr(nr), idenspx(nx), idensx1(Nx), idensR1(NR))
 allocate(psi_out_x(NR,Nx),psi_out_xR(NR,Nx))
 allocate(psi_out_x_Rpx(NR,Nx))
 allocate(psi_g1(NR,Nx),psi_dum(NR,Nx))
 allocate(psi_out_x_intervals(NR,Nx),rho_psi_out_x_intervals(NR,Nx))
 allocate(psi_out_x1(NR,Nx),psi_out_xR1(NR,Nx))
 allocate(psi_dumx(Nx), psi_dumR(NR))
! allocate(psi_out_xl(NR,Nx),psi_out_xr(NR,Nx))
 allocate(psi_out_R(NR,Nx),psi_out_Rx(NR,Nx))
 allocate(psi_out_R1(NR,Nx),psi_out_Rx1(NR,Nx),kprop_vel(NR,Nx))
 allocate(psi_out_x_momt(NR,Nx), JES_en_out_x(NR/2,Nx/2))
 allocate(JES_en_out_xR(NR/2,Nx/2),JES_en_out_Rx(NR/2,Nx/2))
! allocate(psi_out_Rl1(NR,Nx),psi_out_Rr1(NR,Nx))
 !allocate(psi_out_x101(NR,Nx),psi_out_x151(NR,Nx))
 

 print*
 print*,'Tuning FFTW...'

 void=fftw_init_threads( )
 if (void==0) then
    write(*,*) 'Error in fftw_init_threads, quitting'
    stop
 endif

 call fftw_plan_with_nthreads(omp_get_max_threads())
 call dfftw_plan_dft_2d(planF, NR, Nx, psi, psi, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_2d(planB, NR, Nx, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)
 call dfftw_plan_dft_2d(planFd, NR, Nx, psi_dum, psi_dum, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_2d(planBd, NR, Nx, psi_dum, psi_dum, FFTW_BACKWARD,FFTW_MEASURE)

 call dfftw_plan_dft_1d(planFx, Nx, psi_dumx, psi_dumx, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_1d(planBx, Nx, psi_dumx, psi_dumx, FFTW_BACKWARD,FFTW_MEASURE)
 call dfftw_plan_dft_1d(planFR, NR, psi_dumR, psi_dumR, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_1d(planBR, NR, psi_dumR, psi_dumR, FFTW_BACKWARD,FFTW_MEASURE)

 print*,'Done.'
 print*


 do I = 1, NR! pR**2 /2* red_mass
   v12=1.0d0/R(I)
   !vcol_prop(I)=exp(-im*dt*v12)
   vcol_prop(I)=exp(-im*0.5d0*dt*v12)
   do J = 1, Nx ! px**2 / 2* m_eff 
    kprop(I,J) = exp(-im *dt *((PR(I)**2/(2.d0*m_red))+(Px(J)**2/(2.d0 *m_eff))))
   end do
 end do
 
do I=1,Nx
  g(I)=(exp(-1.0*(x(I)-(25.0d0))**2))+exp(-1.0*(x(I)+(25.0d0))**2)
  write(198,*) x(I)*au2a, g(I)
enddo


 do I = 1, NR
  do J = 1, Nx
   !psi(I,J) =ewf(J,I,1) * exp(kappa*(R(I)-RI)**2)
   psi(I,J) =ewf(J,I,1) * chi0(I,1) !exp(kappa*(R(I)-RI)**2)
  end do
 end do 

 !cpm=40.0d0!/au2a
 !cpmx=105.0d0!/au2a !detector
 !cpmxl=20.d0/au2a
 !cpmR=2.0d0/au2a
 I_cpm = minloc(abs(x(NR/2:NR)-cpm),1) 
 I_cpmR = minloc(abs(R(:)-cpmR),1) 

 cof_choice = 1
 call cofs(cof_choice,cof, cofxl,cofxr, cofR)
 cof_choice = 2
 call cofs(cof_choice, cofx,cofxl, cofxr,cofR)
! call cofs(cpm,cpmR,cof, cofxl,cofxr, cofR)
! call cofs(cpmx,cpmR,cofx,cofxl, cofxr,cofR)
! call cofs(cpmxl,cpmR,cofxl,cofR)

  do I = 1, Nx
    write(96,*) sngl(x(I)), sngl(cof(i))
    write(97,*) x(I), cofx(I)
    write(98,*) x(I), cofxr(I)
  end do

  do j = 1, NR
   write(99,*) sngl(R(J)), sngl(cofR(j))
  end do
 close(96, status='keep')
 close(97, status='keep')
 close(98, status='keep')
 close(99, status='keep')

 call integ_2d(psi, norm)
 print*,'norm1 =', sngl(norm)

 psi = psi / sqrt(norm)

 call integ_2d(psi, norm)
 print*,'norm2 =', sngl(norm)

 
 do I = 1, NR/2
   do J = 1, Nx, 4
    write(100,*) sngl(R(I) *au2a), sngl(x(J) *au2a),sngl(abs(psi(I,J))**2)!density?
   end do
    write(100,*)
 end do

 psi0 = psi

!______________________________________________________________________
!
!                   Propagation Loop
!______________________________________________________________________


 print*
 print*,'2D propagation...'
 print*

! !psi_out initiation___!
! do J=1,nx
!    psi_out_x(:,J)=psi(:,J)*(1.0d0-cof(J))
! enddo
! psi_dum=psi_out_x
! call dfftw_execute(planFd)
! psi_out_x=psi_dum
!
! do I=1,NR
!    psi_out_R(I,:)=psi(I,:)*(1.0d0-cofR(I))
! enddo
! psi_dum=psi_out_R
! call dfftw_execute(planFd)
! psi_out_R=psi_dum
! ! Done________________!
! call pulse(E2,A)
 LL=10000
 psi_out_x = (0.d0,0.d0)
 psi_out_R = (0.d0,0.d0)
 psi_out_xR = (0.d0,0.d0)
 psi_out_Rx = (0.d0,0.d0)
 !starting timeloop____!
 timeloop: do K = 1, Nt

    time = K * dt

    if(mod(K,10).eq.0) then
      print*,'time:', sngl(time *au2fs)
    end if

   evR = 0.d0
   evx = 0.d0
   epR = 0.d0
   epx = 0.d0
   
!    E21 = E01 *exp(-fwhm1 * (time - t_start1)**2)* (cos(omega1 * (time-t_start1)+phi1))! &
      ! & - ((2.d0 * fwhm) / omega1) * (time - t_start) * sin(omega1 * time+phi1))
    
!    E22 = E02 *exp(-fwhm2 * (time - t_start2)**2)* (cos(omega2 * (time-t_start2)+phi2))! &
      ! & - ((2.d0 * fwhm) / omega2) * (time - t_start) * sin(omega2 * time+phi2))
    
!    E2(K)=E21+E22  

! !vector field
!  A(k)=(-1.0)*sum(E2(1:K))*dt
!  E2_New(k)=-(A(k)-A(k-1))/dt
!  write(906,*) time*au2fs, A(k), E2_New(k)
! !........................

!!###############################################
!!##### length gauge propagation ################
!!###############################################
  call length_gauge_propagation(planF, planB, E2(K))
!      do I =1,NR 
!        psi(I,:) = psi(I,:) * exp(-0.5d0*im * dt * (pot(I,:)+x(:)*E2(k)))
!      enddo
!      call dfftw_execute(planF)
!      psi = psi * kprop
!      psi = psi / sqrt(dble(NR * Nx))
!      call density(psi,idenspR,idenspx)
!      call dfftw_execute(planB)
!      psi = psi / sqrt(dble(NR * Nx))
!      do I =1,NR 
!        psi(I,:) = psi(I,:) * exp(-0.5d0*im * dt * (pot(I,:)+x(:)*E2(k)))
!      enddo
!!################################################
!!##### velocity gauge propagation ###############
!!################################################
!      psi = psi * exp(-0.5d0*im * dt * pot)
!
!      call dfftw_execute(planF)
!!      psi = psi * kprop
!      psi = psi / sqrt(dble(NR * Nx))
! 
!      do J=1,Nx
!         psi(:,J)=psi(:,J)* exp(-im*dt*((pR(:)**2 /(2*m_red)) + &
!                 & (px(J) + A(K))**2 /(2*m_eff)))
!      enddo
!
!!      call density(psi,idenspR,idenspx)
!!      do I = 1, NR
!!       epR = epR + dble(pR(i) *idenspR(I))
!!      end do
!!       epR = epR * dR/ norm
!!      do J = 1, Nx
!!       epx = epx + dble(px(j) *idenspx(J))
!!      end do
!!       epx = epx * dx/ norm
!       
!      call dfftw_execute(planB)
!      psi = psi / sqrt(dble(NR * Nx))
!
!      psi = psi * exp(-0.5d0* im * dt * pot)
! 
! ............................................

!+++++++ Gaussian window +++++++++++++++++++++++  
  do I =1, NR
   psi_g1(I,:)=g(:)*psi(I,:)*exp(-im*x(:)*A(k))
  enddo
   psi_dum=psi_g1
   call dfftw_execute(planFd)
   psi_g1=psi_dum/sqrt(dble(Nx*NR))
   
  call density(psi_g1, idensR, idensx)
   momt_x=sum(dble(px(:)*idensx(:)))*dx !/sum(abs(psi_g1(:)**2))
   write(799,*) time*au2fs, A(K), momt_x, sum(abs(psi_g1(:,:))**2)*dx*dR
!+++++++++++++++++++++++++++++++++++++++++++++++


  call density(psi,idensR,idensx)
  call integ_2D(psi, norm)

   do I = 1, NR
    evR = evR + dble(R(I) *idensR(I))
   end do
    evR = evR * dR

   do J = 1, Nx
    evx = evx + dble(x(J) *idensx(J))
   end do
    evx = evx * dx

    evR = evR / norm
    evx = evx / norm

! call overlap_2D(psi0, psi, corrf)
! call pop_analysis(psi, time, ewf, K)
 call localization(psi, time, ewf, K)

    write(800,*) sngl(time *au2fs), sngl(evR)!, sngl(epR)
    write(801,*) sngl(time *au2fs), sngl(evx)!, sngl(epx)
    write(802,*) sngl(time *au2fs), sngl(norm) 
!    write(803,*) sngl(time *au2fs), sngl(abs(corrf))
    write(909,*) sngl(time *au2fs), sngl(E2(k))
!    write(907,*) sngl(time *au2fs), sngl(E21)
!    write(908,*) sngl(time *au2fs), sngl(E22)
  
  
!  if(mod(K,1000).eq.0) then
!   LL=LL+1
!  do I= 1, NR
!    R=R0+(I-1)*dR
!    do J=1, Nx
!    x=x0+(J-1)*dx
!    write(LL,*) R*au2a, x*au2a, abs(psi(I,J))**2
!    enddo
!    write(LL,*)
!  enddo
!  close(LL)
!  endif   
  
 if(mod(K,50).eq.0) then
  do I = 1, NR, 4
    write(200,*) sngl(time *au2fs), sngl(R(I)),sngl(idensR(I))
  end do
   write(200,*)
  
   do I = NR/2+1, NR, 4
    write(202,*) sngl(time *au2fs), sngl(Pr(i)),sngl(idenspR(I))
   end do
   do I = 1, NR/2, 4
    write(202,*) sngl(time *au2fs), sngl(pr(i)),sngl(idenspR(I))
   end do
   write(202,*)
 end if
 if(mod(K,20).eq.0) then
  do J = 1, Nx, 4
    write(201,*) sngl(time *au2fs), sngl(x(J)), sngl(idensx(J))
  end do
    write(201,*)

  do j = Nx/2 +1, Nx, 4
    write(203,*) sngl(time *au2fs), sngl(Px(j)),sngl(idenspx(j))
  end do
  do j = 1, Nx/2, 4
    write(203,*) sngl(time *au2fs), sngl(px(j)),sngl(idenspx(j))
  end do
   write(203,*)
 end if
 
  
  
 !...............CONTINUUM....................................
 !projectiong continuum states       !ref : https://iopscience.iop.org/article/10.1088/1361-6455/ab8c21/pdf
!  psi_cont=0.0d0 
!  do M=1,Nstates
!    a1e_M=0.0d0
!    do I=1,NR
!     do J=1,Nx
!      a1e_m(I)=a1e_M(I)+ewf(J,I,M)*psi(I,J)
!     enddo
!     a1e_M(I)=a1e_M(I)*dx
!     do J=1,Nx
!      psi_cont(I,J)=psi_cont(I,J)+a1e_M(I)*ewf(J,I,M)
!     enddo
!    enddo
!  enddo
!  psi_cont(:,:)=psi(:,:)-psi_cont(:,:)  
!  psi_dum=psi_cont
!  call dfftw_execute(PlanFd)                  
!  psi_cont=psi_dum
!  psi_cont=psi_cont/sqrt(dble(NR*Nx))
!  call density(psi_cont, idensR, idensx)
!  call integ_2D(psi_cont, norm_outx)
!  do J=1,Nx
!    momt_x=momt_x+dble(px(J)*idensx(J))
!  enddo
!  momt_x=momt_x*dx 
!  write(3111,*) time*au2fs, a(k), momt_x, norm_outx 
 
! ! Electronic coordinate wavefunction prpopagation____!
!   do J=1,Nx
!    psi_dumx(:)=psi_out_xl(:,J)
!    call dfftw_execute(planBx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xl(:,J)=psi_dumx(:)
!
!    psi_out_xl(:,J)=vcol_prop(:)*psi_out_xl(:,J) !coloumb repulsion propagator
!!  
!    psi_dumx(:)=psi_out_xl(:,J)
!    call dfftw_execute(planFx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xl(:,J)=psi_dumx(:)
!!
!    psi_dumx(:)=psi_out_xr(:,J)
!    call dfftw_execute(planBx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xr(:,J)=psi_dumx(:)
!
!    psi_out_xr(:,J)=vcol_prop(:)*psi_out_xr(:,J) !coloumb repulsion propagator
!!  
!    psi_dumx(:)=psi_out_xr(:,J)
!    call dfftw_execute(planFx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xr(:,J)=psi_dumx(:)
!
!    psi_out_xl(:,j) = psi_out_xl(:,j) * exp(-im * dt *  &
!                   &  px(J) *A(k))
!               !& + (1.d0 + 1.d0/(2.d0*m_red +1.d0)) *x *E2))
!    psi_out_xr(:,j) = psi_out_xr(:,j) * exp(-im * dt *  &
!                   &  px(J) *A(k))
! !                  &  (1.d0 + 1.d0/(2.d0*m_red +1.d0)) *x *E2)
!    psi_out_x_cont(i,j) = psi_out_x_cont(i,j) * exp(-im * dt *  &
!                   &  px(J) *A(K))
! !                  &  (1.d0 + 1.d0/(2.d0*m_red +1.d0)) *x *E2)
!   enddo
!   
!   psi_out_xl=psi_out_xl*kprop
!   psi_out_xr=psi_out_xr*kprop
!   psi_out_x_cont=psi_out_x_cont*kprop

! ! Nuclear coordinate wavefunction propagation________!
!   do I=1,NR
!    do J=1,Nx
!   ! x=x0+(J-1)*dx
!!     psi_out_Rl(i,j) = psi_out_Rl(i,j) * exp(-im * dt *  &
!                   &  px(J) *A(k))
!              !& + (1.d0 + 1.d0/(2.d0*m_red +1.d0)) *x *E2))
!     psi_out_Rr(i,j) = psi_out_Rr(i,j) * exp(-im * dt *  &
!                   &  px(J) *A(k))
!     psi_out_R2(i,j) = psi_out_R2(i,j) * exp(-im * dt *  &
!                   &  px(J) *A(k))
!    enddo
!   enddo
!   psi_out_Rl=psi_out_Rl*kprop
!   psi_out_Rr=psi_out_Rr*kprop
!   psi_out_R2=psi_out_R2*kprop

!############# Basic continuum propagation ############################
!%%%%%% for primary length gauge propagation %%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!****** First propapgating the dissociated parts **********************
  ! Only dissociated
  do I = 1, NR
    psi_out_R(I,:) = psi_out_R(I,:)* exp(-im*0.5d0*dt*pot(NR-I_cpmR,:))
    psi_dumx(:) = psi_out_R(I,:)
    call dfftw_execute(planFx)
    psi_out_R(I,:) = psi_dumx(:)/sqrt(dble(Nx))
  enddo
  do J =1,Nx
  psi_out_R(:,J) = psi_out_R(:,J) &
          & * exp(-im * dt * (pR(:)**2 /(2*m_red)+ (px(J)+A(K))**2 /(2*m_eff)))
  enddo
  do I = 1, NR
    psi_dumx(:) = psi_out_R(I,:)
    call dfftw_execute(planBx)
    psi_out_R(I,:) = psi_dumx(:)/sqrt(dble(Nx))
    psi_out_R(I,:) = psi_out_R(I,:)* exp(-im*0.5d0*dt*pot(NR-I_cpmR,:))
  enddo
! Note these are the velocity gauge wfs
  do J = 1, Nx/2
  psi_out_xR(:,J) = psi_out_xR(:,J) * exp(-im*dt*pot(NR-I_cpmR,I_cpm)) &
          & * exp(-im * dt * (pR(:)**2 /(2*m_red)+ (px(J)+A(K))**2 /(2*m_eff)))
  psi_out_Rx(:,J) = psi_out_Rx(:,J) * exp(-im*dt*pot(NR-I_cpmR,I_cpm)) &
          & * exp(-im * dt * (pR(:)**2 /(2*m_red)+ (px(J)+A(K))**2 /(2*m_eff)))
  enddo
  do J = Nx/2 +1, Nx
  psi_out_xR(:,J) = psi_out_xR(:,J) * exp(-im*dt*pot(NR-I_cpmR,Nx-I_cpm)) &
          & * exp(-im * dt * (pR(:)**2 /(2*m_red)+ (px(J)+A(K))**2 /(2*m_eff)))
  psi_out_Rx(:,J) = psi_out_Rx(:,J) * exp(-im*dt*pot(NR-I_cpmR,Nx-I_cpm)) &
          & * exp(-im * dt * (pR(:)**2 /(2*m_red)+ (px(J)+A(K))**2 /(2*m_eff)))
  enddo

!****** now the propagation of ionized parts ***************************
 ! Note 
 ! in the nucl-nucl repulsion potential
 do J = 1, Nx/2 
 ! psi_out_x(:,J) = psi_out_x(:,J) * vcol_prop(:)!* &
  psi_out_x(:,J) = psi_out_x(:,J) * exp(-im*0.5*dt*pot(:,I_cpm))!* &
           ! & exp(-0.5d0*im * dt * (x(:)*E2(k)))
  !%% R-grid FFT
  psi_dumR(:) = psi_out_x(:,J)   
  call dfftw_execute(planFR)
  psi_out_x(:,J) = psi_dumR(:)/sqrt(dble(NR))
 enddo
 do J = Nx/2 +1, Nx 
 ! psi_out_x(:,J) = psi_out_x(:,J) * vcol_prop(:)!* &
  psi_out_x(:,J) = psi_out_x(:,J) * exp(-im*0.5*dt*pot(:,Nx-I_cpm))!* &
           ! & exp(-0.5d0*im * dt * (x(:)*E2(k)))
  !%% R-grid FFT
  psi_dumR(:) = psi_out_x(:,J)   
  call dfftw_execute(planFR)
  psi_out_x(:,J) = psi_dumR(:)/sqrt(dble(NR))
 enddo
  !%%
  !%% kinetic operator
  psi_out_x_Rpx=(0.0d0,0.d0)
 do J=1,Nx
  psi_out_x(:,J) = psi_out_x(:,J) &
          & * exp(-im * dt * (pR(:)**2 /(2*m_red)+ (px(J)+A(K))**2 /(2*m_eff)))
 enddo
  !%%
  !%% R-grid FFT
 do J = 1, Nx/2
  psi_dumR(:) = psi_out_x(:,J)
  call dfftw_execute(planBR)
  psi_out_x(:,J) = psi_dumR(:)/sqrt(dble(NR))! * Nx))
  
  !%%
  !%% Potential operator
 ! psi_out_x(:,J) = psi_out_x(:,J) * vcol_prop(:) !* &
  psi_out_x(:,J) = psi_out_x(:,J) * exp(-im*0.5*dt*pot(:,I_cpm))!* &
            !& exp(-0.5d0*im * dt * (x(:)*E2(k)))
  psi_dum(:,J) =psi_out_x(:,J)*exp(im*A(K)*x(J))
 enddo
 do J = Nx/2 +1, Nx
  psi_dumR(:) = psi_out_x(:,J)
  call dfftw_execute(planBR)
  psi_out_x(:,J) = psi_dumR(:)/sqrt(dble(NR))! * Nx))
  
  !%%
  !%% Potential operator
 ! psi_out_x(:,J) = psi_out_x(:,J) * vcol_prop(:) !* &
  psi_out_x(:,J) = psi_out_x(:,J) * exp(-im*0.5*dt*pot(:,Nx-I_cpm))!* &
            !& exp(-0.5d0*im * dt * (x(:)*E2(k)))
  psi_out_x_Rpx(:,J) =psi_out_x(:,J)*exp(im*A(K)*x(J))
 enddo

 call density(psi_out_x_Rpx,idensR,idensx)
 if(mod(K,100).eq.0) then
  do I = 1, NR, 2
    write(206,*) sngl(time *au2fs), sngl(R(I)),sngl(idensR(I))
  end do
  write(206,*)
  
  do I = NR/2+1, NR, 2
    write(207,*) sngl(time *au2fs), sngl(Pr(i)),sngl(idenspR(I)), &
            & sum(abs(psi_out_xR(I,:))**2)*dx, sum(abs(psi_out_Rx(I,:))**2)*dx
  end do
  do I = 1, NR/2, 2
    write(207,*) sngl(time *au2fs), sngl(pr(i)),sngl(idenspR(I)), &
            & sum(abs(psi_out_xR(I,:))**2)*dx, sum(abs(psi_out_Rx(I,:))**2)*dx
  end do
  write(207,*)
  
  do J = 1, Nx, 4
    write(204,*) sngl(time *au2fs), sngl(x(J)), sngl(idensx(J))
  end do
  write(204,*)

  do j = Nx/2 +1, Nx, 4
    write(205,*) sngl(time *au2fs), sngl(Px(j)),sngl(idenspx(j)), &
            & sum(abs(psi_out_xR(:,J))**2)*dR, sum(abs(psi_out_Rx(:,J))**2)*dR
  end do
  do j = 1,Nx/2, 4
  write(205,*) sngl(time *au2fs), sngl(px(j)),sngl(idenspx(j)), &
            & sum(abs(psi_out_xR(:,J))**2)*dR, sum(abs(psi_out_Rx(:,J))**2)*dR
  end do
  write(205,*)
 end if

 if (mod(K,2000) .eq. 0 .and. k .ge. 4000) then
  if (time*au2fs .gt. 99.99d0) then
    write(filename, '(a,f6.2,a)') 'psi_outx_R-px_time_',time*au2fs,'fs.out'
    open(2000,file=filename, status='unknown') 
  else
    write(filename, '(a,f5.2,a)') 'psi_outx_R-px_time_',time*au2fs,'fs.out'
    open(2000,file=filename, status='unknown') 
  endif

  do I = 1, NR, 4
     do J = Nx/2 +1, Nx, 4
        write(2000,*) R(I), px(J), abs(psi_out_x_Rpx(I,J))**2
     enddo
     do J = 1, Nx/2, 4 
        write(2000,*) R(I), px(J), abs(psi_out_x_Rpx(I,J))**2
     enddo
     write(2000,*)
  enddo
!  enddo
  close(2000)
  psi_dum=(0.d0,0.d0)
!  psi_dumx=(0.d0,0.d0)
 endif


!%%%%%%%%%%%%%%%%%%%% applying cutoff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

! do J = 1, Nx

!    psi_dumx(:)=psi_out_xl(:,J)
!    call dfftw_execute(planBx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xl(:,J)=psi_dumx(:)
!
!    psi_dumx(:)=psi_out_xr(:,J)
!    call dfftw_execute(planBx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xr(:,J)=psi_dumx(:)
!
   do I = 1, NR

!######## commented fancy stuff, dividing left and right absrober ##############
!   psi_out_xl1(i,j)=psi(i,j)*(1.0d0-cofxl(j))   !ionized wavefunction
!   psi_out_xr1(i,j)=psi(i,j)*(1.0d0-cofxr(j))   !ionized wavefunction
!   psi_out_Rl1(i,j)=psi_out_xl(i,j)*(1.0d0-cofR(i)) !ionized and then dissociated wf; cut out of psi_out_x  !Let's try not
!   psi_out_Rr1(i,j)=psi_out_xr(i,j)*(1.0d0-cofR(i)) !propogating psi_out_x

!   psi_out_x_cont1(i,j)=psi(i,j)*(1.0d0-cof(j))  !using absorber for countinuum state projection

!   psi_out_xl(i,j)=psi_out_xl(i,j) *cofR(i) !dissociated ionized wf cut out of psi_out_x 
!   psi_out_xr(i,j)=psi_out_xr(i,j) *cofR(i) !dissociated ionized wf cut out of psi_out_x 

!###### Let's turn back the clock and use simple method #########################
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !% Dissociative ionization %%%%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   psi_out_xR1(I,:) = psi_out_x(I,:)*(1.0d0-cofR(I)) ! dissociative ionization
   psi_out_x(I,:) = psi_out_x(I,:)*cofR(I)   ! cutting off the dissociated wavefunction
 !  psi_out_x_intervals(I,:) = psi_out_x_intervals(I,:)*cofR(I)   ! cutting off the dissociated wavefunction
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !%% localization and dissociation
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   psi_out_Rx1(I,:)=psi_out_R(I,:)*(1.0d0-cof(:)) ! ionization after dissociation
   psi_out_R(I,:)=psi_out_R(I,:)*cof(:) ! cutting off ionized part from dissociated wavefunction
   psi_out_x1(i,:) = psi(i,:)*(1.d0-cof(:))  ! ionization
   psi_out_R1(i,:)=psi(i,:)*(1.0d0-cofR(i)) ! dissociating localized wavefunction
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !%% Total cut off %%%%%%%%%%%%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   psi(i,:) = psi(i,:) * cof(:)*cofR(i)   !cut by cofR is important otherwise reflections
   end do

! Converting escapees to velocity guage... you go guys ;)
 do J= 1,Nx
   psi_out_x1(:,J) = psi_out_x1(:,J) * exp(-im*x(J)*A(K))  
   psi_out_R1(:,J) = psi_out_R1(:,J) * exp(-im*x(J)*A(K))  
 enddo


!    psi_dumx(:)=psi_out_xl(:,J)
!    call dfftw_execute(planFx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xl(:,J)=psi_dumx(:)
!
!    psi_dumx(:)=psi_out_xr(:,J)
!    call dfftw_execute(planFx)
!    psi_dumx=psi_dumx/sqrt(dble(NR))
!    psi_out_xr(:,J)=psi_dumx(:)
!  end do

 ! Fourier transform of the cutoffs____________________!
!##### Back to basics ####################################
  do J=1,Nx
   psi_dumR(:)=psi_out_R1(:,J)
   call dfftw_execute(PlanFR)
   psi_dumR=psi_dumR / sqrt(dble(NR))
   psi_out_R1(:,J)= psi_dumR(:)
  enddo

  do I =1,NR
   psi_dumx(:)=psi_out_Rx1(I,:)
   call dfftw_execute(PlanFx)
   psi_out_Rx1(I,:)=psi_dumx(:)/sqrt(dble(Nx))
  enddo

  do J =1,Nx
   psi_dumR(:)=psi_out_xR1(:,J)
   call dfftw_execute(PlanFR)
   psi_out_xR1(:,J)=psi_dumR(:)/sqrt(dble(NR))
  enddo

! print*, "test"
  do I=1,NR
   psi_dumx(:)=psi_out_x1(I,:)
   call dfftw_execute(PlanFx)
   psi_dumx=psi_dumx / sqrt(dble(Nx))
   psi_out_x1(I,:)= psi_dumx(:)
  enddo
  

! print*, "test"
 ! Accumulation________________________________________!
!  psi_out_xl=psi_out_xl+psi_out_xl1
!  psi_out_xr=psi_out_xr+psi_out_xr1
!  psi_out_x_cont=psi_out_x_cont+psi_out_x_cont1

!  psi_out_x=psi_out_x+psi_out_x1 * exp(im*dt*x(:)*A(K)) ! in velocity gauge
!%%%%%% coneverting dissociating wf from length gauge to velocity gauge
  do I=1,NR
  psi_out_x(I,:)=psi_out_x(I,:)+psi_out_x1(I,:)! * exp(im*dt*x(:)*A(K)) ! in velocity gauge
  psi_out_xR(I,:)=psi_out_xR(I,:) + psi_out_xR1(I,:)   ! alredy velocity gauge
  psi_out_R(I,:)=psi_out_R(I,:) + psi_out_R1(I,:)
  psi_out_Rx(I,:)=psi_out_Rx(I,:) + psi_out_Rx1(I,:)   ! alredy velocity gauge
  enddo

  
  if(mod(K+2,2000) .eq. 0) then
     rho_psi_out_x_intervals = (0.0d0,0.0d0)
  else
  do J = 1, Nx
     psi_dumR(:) = psi_out_x1(:,J)   
     call dfftw_execute(planFR)
     psi_dumR(:) = psi_dumR(:)/sqrt(dble(NR))
     rho_psi_out_x_intervals(:,J) = rho_psi_out_x_intervals(:,J) + abs(psi_dumR(:))**2
  enddo
  endif

 if (mod(K,2000) .eq. 0 .and. k .ge. 4000) then
  if (time*au2fs .gt. 99.99d0) then
    write(filename, '(a,f6.2,a)') 'psi_out_x_time_',time*au2fs,'fs.out'
    open(2000,file=filename, status='unknown')
  else
    write(filename, '(a,f5.2,a)') 'psi_out_x_time_',time*au2fs,'fs.out'
    open(2000,file=filename, status='unknown')
  endif

  do I = NR/2 +1,NR, 4
     do J = Nx/2 +1, Nx, 4
        write(2000,*) pR(I), px(J), rho_psi_out_x_intervals(I,J) 
     enddo
     do J = 1, Nx/2, 4
        write(2000,*) pR(I), px(J), rho_psi_out_x_intervals(I,J)
     enddo
     write(2000,*)
  enddo
  do I = 1, NR/2, 4
     do J = Nx/2 +1, Nx, 4
        write(2000,*) pR(I), px(J), rho_psi_out_x_intervals(I,J)
     enddo
     do J=1, Nx/2, 4
        write(2000,*) pR(I), px(J), rho_psi_out_x_intervals(I,J)
     enddo
     write(2000,*)
  enddo
  close(2000)
  psi_dum=(0.d0,0.d0)
  psi_dumx=(0.d0,0.d0)
 endif

! Joint electron nuclear spectra ionization (psi_out_x)
 if (mod(K,2000) .eq. 0 .and. k .ge. 4000) then
  if (time*au2fs .gt. 99.99d0) then 
    write(filename, '(a,f6.2,a)') 'JENS_momt_ionized-wf_time_',time*au2fs,'fs.out'
    open(2230,file=filename, status='unknown')
    write(filename, '(a,f6.2,a)') 'JENS_en_ionized-wf_time_',time*au2fs,'fs.out'
    open(3330,file=filename, status='unknown')
    write(filename, '(a,f6.2,a)') 'JENS_x-momt_R-en_ionized-wf_time_',time*au2fs,'fs.out'
    open(3340,file=filename, status='unknown')
  else
    write(filename, '(a,f5.2,a)') 'JENS_momt_ionized-wf_time_',time*au2fs,'fs.out'
    open(2230,file=filename, status='unknown')
    write(filename, '(a,f5.2,a)') 'JENS_en_ionized-wf_time_',time*au2fs,'fs.out'
    open(3330,file=filename, status='unknown')
    write(filename, '(a,f5.2,a)') 'JENS_x-momt_R-en_ionized-wf_time_',time*au2fs,'fs.out'
    open(3340,file=filename, status='unknown')
  endif
  psi_out_x_momt = psi_out_x
  do J = 1, Nx
  psi_dumR(:) = psi_out_x(:,J)   
  call dfftw_execute(planFR)
  psi_dumR(:) = psi_dumR(:)/sqrt(dble(NR))
  psi_out_x_momt(:,J) = psi_dumR(:)
  enddo

!  call integ_2d(psi_out_x_momt, norm_outx)
!  psi_out_x_momt = psi_out_x_momt/sqrt(norm_outx)
!  call integ_2d(psi_out_xR_momt, norm_outxR)
!  psi_out_xR_momt = psi_out_xR_momt/sqrt(norm_outxR)
!  call integ_2d(psi_out_Rx_momt, norm_outRx)
!  psi_out_Rx_momt = psi_out_Rx_momt/sqrt(norm_outRx)
  do I = 1, NR/2
    do J = 1, NR/2 
      JES_en_out_x(I,J) =(m_eff*m_red) * ((abs(psi_out_x_momt(I,J))**2)/(pR(I)*px(J)) &
                 & + (abs(psi_out_x_momt(I,NR-J))**2)/(pR(I)*px(NR-J)))
      JES_en_out_xR(I,J) =(m_eff*m_red) * ((abs(psi_out_xR(I,J))**2)/(pR(I)*px(J)) &
                 & + (abs(psi_out_xR(I,NR-J))**2)/(pR(I)*px(NR-J)))
      JES_en_out_Rx(I,J) =(m_eff*m_red) * ((abs(psi_out_Rx(I,J))**2)/(pR(I)*px(J)) &
                 & + (abs(psi_out_Rx(I,NR-J))**2)/(pR(I)*px(NR-J)))
      write(3330,*) pR(I)*pR(I)/(2*m_red), px(J)*px(J)/2, &
                  & JES_en_out_x(I,J), JES_en_out_xR(I,J), JES_en_out_Rx(I,J)
    enddo
    write(3330,*)
  enddo
  do I=NR/2+1,NR, 4
   do J=Nx/2+1,Nx, 4
     write(2230,*) pR(I),px(J), abs(psi_out_x_momt(I,J))**2
   enddo
   do J=1, Nx/2, 4
     write(2230,*) pR(I), px(J), abs(psi_out_x_momt(I,J))**2
   enddo
   write(2230,*)
  enddo
  do I=1,NR/2, 4
   do J=Nx/2+1,Nx
     write(2230,*) pR(I),px(J), abs(psi_out_x_momt(I,J))**2
     write(3340,*) pR(I)*pR(I)/(2*m_red), px(J), &
                  & m_red*abs(psi_out_x_momt(I,J))**2 / pR(I), &
                  & m_red*abs(psi_out_xR(I,J))**2 / pR(I), &
                  & m_red*abs(psi_out_Rx(I,J))**2 / pR(I)
   enddo
   do J=1, Nx/2, 4
     write(2230,*) pR(I), px(J), abs(psi_out_x_momt(I,J))**2
     write(3340,*) pR(I)*pR(I)/(2*m_red), px(J), &
                  & m_red*abs(psi_out_x_momt(I,J))**2 / pR(I), &
                  & m_red*abs(psi_out_xR(I,J))**2 / pR(I), &
                  & m_red*abs(psi_out_Rx(I,J))**2 / pR(I)
   enddo
   write(2230,*)
   write(3340,*)
  enddo
  close(2230)
  close(3330)
  close(3340)
  psi_dumR = (0.d0,0.d0)
  endif

! do I=1,NR
!  do J=1,Nx
!  psi_out_R2(i,j)=psi_out_R21(i,j)+psi_out_R2(i,j)!*exp(-im*pot(445,j)*dt)
!!  psi_out_Rl(i,j)=psi_out_Rl1(i,j)+psi_out_Rl(i,j)*exp(-im*dt/25.0d0)        !potential at 13 A (25 au)
!!  psi_out_Rr(i,j)=psi_out_Rr1(i,j)+psi_out_Rr(i,j)*exp(-im*dt/25.0d0)    
!!  psi_out_R(i,j)=psi_out_R(i,j)+psi_out_R1(i,j)*exp(-im*pot(410,j)*dt) &
!!          & + psi_out_R2(i,j)*exp(-im*pot(410,j)*dt)
!  enddo
! enddo

 ! psi_cont=psi_cont+psi_out_x_cont
 
 ! Analysis___________________________________________!
 ! ionizaion probability
 call integ_2D(psi_out_x, norm_outx)
 call integ_2D(psi_out_x1, norm_outx1)
 call integ_2D(psi_out_xR, norm_outxR)
 call integ_2D(psi_out_xR1, norm_outxR1)
 call integ_2D(psi, norm)
 norm_addx=norm_addx+norm_outx1
 write(850,*) time*au2fs, norm_outx, norm_outx1, norm_addx, norm_addx1, &
         & norm_outxR, norm_outxR1, norm

 !dissociation probability
 call integ_2D(psi_out_R, norm_outR)
 call integ_2D(psi_out_R1, norm_outR1)
 call integ_2D(psi_out_Rx, norm_outRx)
 call integ_2D(psi_out_Rx1, norm_outRx1)
 norm_addR=norm_addR+norm_outR1
 write(851,*) time*au2fs, norm_outR, norm_outR1, norm_addR, norm_addR1, &
         & norm_outRx, norm_outRx1
 ! Electronic Momentum
!  !$OMP parallel do
!  do J=1,Nx
!   psi_dumR(:)=psi_out_x1(:,J)
!   call dfftw_execute(PlanFR)
!   psi_dumR=psi_dumR / sqrt(dble(NR))
!   psi_out_x1(:,J)= psi_dumR(:)
!  enddo
!  !$OMP end parallel do

  call density(psi_out_x1,idenspR,idenspx)
  if(mod(K,100).eq. 0) then
  do j = Nx/2+1, Nx
    write(1114,*) sngl(time *au2fs), sngl(Px(j)),sngl(idenspx(j))
  end do
  do j = 1,Nx/2
    write(1114,*) sngl(time *au2fs), sngl(px(j)),sngl(idenspx(j))
  end do
   write(1114,*)
  endif

!  psi_dum=psi_out_x
!  call dfftw_execute(PlanFd)
!  psi_dum=psi_dum/sqrt(dble(NR*Nx))
  call density(psi_out_R,idenspR,idenspx)
  if(mod(K,100).eq. 0) then
  do j = Nx/2+1, Nx
    write(1124,*) sngl(time *au2fs), sngl(Px(j)),sngl(idenspx(j))
  end do
  do j = 1,Nx/2
    write(1124,*) sngl(time *au2fs), sngl(px(j)),sngl(idenspx(j))
  end do
   write(1124,*)
  endif
! !Left Momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  momt_x=0.0d0
!  momt_x1=0.0d0
!  call density(psi_out_xl, idensR, idensx)
!  call integ_2D(psi_out_xl, norm_outx)
!  call density(psi_out_xl1, idensR1, idensx1)
!  call integ_2D(psi_out_xl1, norm_outxl)
!  do J=1,Nx
!    momt_x=momt_x+dble(px(J)*idensx(J))
!    momt_x1=momt_x1+dble(px(J)*idensx1(J))
!  enddo
!    momt_x=momt_x*dx
!    momt_x1=momt_x1*dx
!    write(1114,*) time*au2fs,a(K), momt_x, norm_outx, momt_x1, norm_outxl
!
! ! right momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  momt_x=0.0d0
!  momt_x1=0.0d0
!  call density(psi_out_xr, idensR, idensx)
!  call integ_2D(psi_out_xr, norm_outx)
!  call density(psi_out_xr1, idensR1, idensx1)
!  call integ_2D(psi_out_xr1, norm_outxr)
!  do J=1,Nx
!    momt_x=momt_x+dble(px(J)*idensx(J))
!    momt_x1=momt_x1+dble(px(J)*idensx1(J))
!  enddo
!    momt_x=momt_x*dx
!    momt_x1=momt_x1*dx
!    write(1124,*) time*au2fs, a(k), momt_x, norm_outx, momt_x1, norm_outxr
! 
! ! continuum momemtum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  momt_x=0.0d0 
!  call density(psi_cont, idensR, idensx)
!  call integ_2D(psi_cont, norm_outx)
!  do J=1,Nx
!    momt_x=momt_x+dble(px(J)*idensx(J))
!  enddo
!    momt_x=momt_x*dx
!    write(3112,*) time*au2fs, a(k), momt_x, norm_outx
!
!
!
! ! Nuclear Momentum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!  momt_R=0.0d0
!  call density(psi_out_R2, idensR, idensx)
!  call integ_2D(psi_out_R2, norm_outR)
!  do I=1,NR
!    momt_R=momt_R+dble(pR(I)*idensR(I))
!  enddo
!  momt_R=momt_R*dR
!!  momt_R=momt_R/norm_outR
!  write(1113,*) sngl(time*au2fs), sngl(momt_R), sngl(norm_outR)
!
!  momt_R=0.0d0
!  call density(psi_out_Rl, idensR, idensx)
!  call integ_2D(psi_out_Rl, norm_outR)
!  do I=1,NR
!    momt_R=momt_R+dble(pR(I)*idensR(I))
!  enddo
!  momt_R=momt_R*dR
!!  momt_R=momt_R/norm_outR
!  write(1213,*) sngl(time*au2fs), a(k),sngl(momt_R), sngl(norm_outR), &
!          & sum(px(:)*idensx(:))*dx
!
!  momt_R=0.0d0
!  call density(psi_out_Rr, idensR, idensx)
!  call integ_2D(psi_out_Rr, norm_outR)
!  do I=1,NR
!    momt_R=momt_R+dble(pR(I)*idensR(I))
!  enddo
!  momt_R=momt_R*dR
!!  momt_R=momt_R/norm_outR
!  write(1223,*) sngl(time*au2fs), a(k),sngl(momt_R), sngl(norm_outR), &
!          & sum(px(:)*idensx(:))*dx
!
!  momt_R=0.0d0
!  call density(psi_out_Rl1, idensR, idensx)
!  call integ_2D(psi_out_Rl1, norm_outR)
!  do I=1,NR
!    momt_R=momt_R+dble(pR(I)*idensR(I))
!  enddo
!  momt_R=momt_R*dR
!  write(1233,*) sngl(time*au2fs), a(k),sngl(momt_R), sngl(norm_outR), &
!          & sum(px(:)*idensx(:))*dx
!
!  momt_R=0.0d0
!  call density(psi_out_Rr1, idensR, idensx)
!  call integ_2D(psi_out_Rr1, norm_outR)
!  do I=1,NR
!    momt_R=momt_R+dble(pR(I)*idensR(I))
!  enddo
!  momt_R=momt_R*dR
!  write(1333,*) sngl(time*au2fs), a(k),sngl(momt_R), sngl(norm_outR), &
!          & sum(px(:)*idensx(:))*dx

end do timeloop

!Adding left and right part of H+ + H+ + e channel
!  psi_out_R=psi_out_Rl+psi_out_Rr

!  call integ_2D(psi_out_xl, norm_outx)
!  psi_out_xl=psi_out_xl/sqrt(norm_outx)
!  call integ_2D(psi_out_xr, norm_outx)
!  psi_out_xr=psi_out_xr/sqrt(norm_outx)
!  call integ_2D(psi_out_x15, norm_outx)
!  psi_out_x=psi_out_x15/sqrt(norm_outx)
!  call integ_2D(psi_out_Rl, norm_outR)
!  psi_out_Rl=psi_out_Rl/sqrt(norm_outR)
!  call integ_2D(psi_out_Rr, norm_outR)
!  psi_out_Rr=psi_out_Rr/sqrt(norm_outR)
!  call integ_2D(psi_out_R2, norm_outR)
!  psi_out_R2=psi_out_R2/sqrt(norm_outR)

  call integ_2D(psi_out_R, norm_outR)
  psi_out_R=psi_out_R/sqrt(norm_outR) ! x-position grid
  print*, "Norm psi_out_R (dissociated wf: H + H+)", norm_outR
  call integ_2D(psi_out_xR, norm_outxR)
  psi_out_xR=psi_out_xR/sqrt(norm_outxR) ! momentum grid
  print*, "Norm psi_out_xR (CEI wf: H+ + H+ + e-)", norm_outxR
  call integ_2D(psi_out_Rx, norm_outRx)
  psi_out_Rx=psi_out_Rx/sqrt(norm_outRx) ! momentum grid
  print*, "Norm psi_out_Rx (dissociated ionized wf: H+ + [H+ + e-])", norm_outRx
  call integ_2D(psi_out_x, norm_outx)
  psi_out_x=psi_out_x/sqrt(norm_outx) ! R-position grid
  print*, "Norm psi_out_x (ionized wf: H_2^2+ + e-)", norm_outx ! should be close to zero

 do J =1,Nx
  psi_dumR(:)=psi_out_x(:,J)
  call dfftw_execute(PlanFR)
  psi_out_x(:,J)=psi_dumR(:)/sqrt(dble(NR))
 enddo

 do I =1,NR
  psi_dumx(:)=psi_out_R(I,:)
  call dfftw_execute(PlanFx)
  psi_out_R(I,:)=psi_dumx(:)/sqrt(dble(Nx))
 enddo

! Momemtum distribution
!  call density(psi_out_xl, idensR, idensx)
!  do J=Nx/2+1, Nx
!   write(1111,*) px(J), idensx(J)
!  enddo
!  do J=1,Nx/2
!   write(1111,*) px(J), idensx(J)
!  enddo
!
!  call density(psi_out_xr, idensR, idensx)
!  do J=Nx/2+1, Nx
!   write(1121,*) sngl(px(J)), sngl(dble(idensx(J)))
!  enddo
!  do J=1,Nx/2
!   write(1121,*) sngl(px(J)), sngl(dble(idensx(J)))
!  enddo

!  call density(psi_out_x15, idensR, idensx)
!  do J=Nx/2+1, Nx
!   write(1131,*) sngl(px(J)), sngl(dble(idensx(J)))
!  enddo
!  do J=1,Nx/2
!   write(1131,*) sngl(px(J)), sngl(dble(idensx(J)))
!  enddo

!  psi_dum=psi_out_x
!  call dfftw_execute(PlanFd)
!  psi_out_x1=psi_dum
!  psi_out_x1=psi_out_x1/sqrt(dble(NR*Nx))

  open(1110,file='psi_out_x_px.out',status='unknown')
  do J=Nx/2+1, Nx
   write(1110,*) px(J), sum(abs(psi_out_x(:,J))**2)*dR, sum(abs(psi_out_xR(:,J))**2)*dR 
  enddo
  do J=1,Nx/2
   write(1110,*) px(J), sum(abs(psi_out_x(:,J))**2)*dR, sum(abs(psi_out_xR(:,J))**2)*dR 
  enddo
  close(1110)

  open(1111,file='psi_out_x_pR.out',status='unknown')
  do I = NR/2+1, NR
   write(1111,*) pR(I), sum(abs(psi_out_x(I,:))**2)*dx, sum(abs(psi_out_xR(I,:))**2)*dx
  enddo
  do I=1,NR/2
   write(1111,*) pR(I), sum(abs(psi_out_x(I,:))**2)*dx, sum(abs(psi_out_xR(I,:))**2)*dx
  enddo

! Joint electron nuclear spectra ionization (psi_out_x)
! Note: the ionized and the dissociated wavepacket (Coulomb explosion)
!%%%% momentum spectra
  open(2210,file='Joint_elec_nucl_momt_spectra_CEI.out',status='unknown')
  open(2240,file='Joint_elec_momt_nucl_en_spectra_CEI.out',status='unknown')
  do I=NR/2+1,NR
   do J=Nx/2+1,Nx
     write(2210,*) pR(I),px(J), abs(psi_out_x(I,J))**2, abs(psi_out_xR(I,J))**2
   enddo
   do J=1, Nx/2
     write(2210,*) pR(I),px(J), abs(psi_out_x(I,J))**2, abs(psi_out_xR(I,J))**2
   enddo
   write(2210,*)
  enddo
  do I=1,NR/2
   do J=Nx/2+1,Nx
     write(2210,*) pR(I),px(J), abs(psi_out_x(I,J))**2, abs(psi_out_xR(I,J))**2
     write(2240,*) pR(I)*pR(I)/m_red, px(J), m_red*abs(psi_out_x(I,J))**2 / pR(I)
   enddo
   do J=1, Nx/2
     write(2210,*) pR(I),px(J), abs(psi_out_x(I,J))**2, abs(psi_out_xR(I,J))**2
     write(2240,*) pR(I)*pR(I)/m_red, px(J), m_red*abs(psi_out_x(I,J))**2 / pR(I)
   enddo
   write(2210,*)
   write(2240,*)
  enddo
  close(2210)

  open(3310,file='Joint_elec_nucl_KER_spectra_CEI.out',status='unknown')
  do I = 1, NR/2
    do J = 1, NR/2 
      JES_en_out_x(I,J) =(m_eff*m_red) * ((abs(psi_out_x(I,J))**2)/(pR(I)*px(J)) &
                 & + (abs(psi_out_x(I,NR-J))**2)/(pR(I)*px(NR-J)))
      JES_en_out_xR(I,J) =(m_eff*m_red) * ((abs(psi_out_xR(I,J))**2)/(pR(I)*px(J)) &
                 & + (abs(psi_out_xR(I,NR-J))**2)/(pR(I)*px(NR-J)))
      write(3310,*) pR(I)*pR(I)/(2*m_red), px(J)*px(J)/2, &
                  & JES_en_out_x(I,J), JES_en_out_xR(I,J)
    enddo
    write(3310,*)
  enddo
  close(3310)
  
  open(1112,file='psi_out_R_px.out',status='unknown')
  do J=Nx/2+1, Nx
   write(1112,*) px(J), sum(abs(psi_out_R(:,J))**2)*dR, sum(abs(psi_out_Rx(:,J))**2)*dR 
  enddo
  do J=1,Nx/2
   write(1112,*) px(J), sum(abs(psi_out_R(:,J))**2)*dR, sum(abs(psi_out_Rx(:,J))**2)*dR 
  enddo
  close(1112)

  open(1113,file='psi_out_R_pR.out',status='unknown')
  do I = NR/2+1, NR
   write(1113,*) pR(I), sum(abs(psi_out_R(I,:))**2)*dx, sum(abs(psi_out_Rx(I,:))**2)*dx
  enddo
  do I=1,NR/2
   write(1113,*) pR(I), sum(abs(psi_out_R(I,:))**2)*dx, sum(abs(psi_out_Rx(I,:))**2)*dx
  enddo

! Joint electron nuclear spectra dissociation after ionization (psi_out_R)
  open(2211,file='Joint_elec_nucl_momt_spectra_diss.out',status='unknown')
  open(2241,file='Joint_elec_momt_nucl_en_spectra_diss.out',status='unknown')
  do I=NR/2+1,NR
   do J=Nx/2+1,Nx
     write(2211,*) pR(I),px(J), abs(psi_out_R(I,J))**2, abs(psi_out_Rx(I,J))**2
   enddo
   do J=1, Nx/2
     write(2211,*) pR(I),px(J), abs(psi_out_R(I,J))**2, abs(psi_out_Rx(I,J))**2
   enddo
   write(2211,*)
  enddo
  do I=1,NR/2
   do J=Nx/2+1,Nx
     write(2211,*) pR(I),px(J), abs(psi_out_R(I,J))**2, abs(psi_out_Rx(I,J))**2
     write(2241,*) pR(I)*pR(I)/m_red, px(J), m_red*abs(psi_out_R(I,J))**2 / pR(I), &
             & m_red*abs(psi_out_Rx(I,J))**2 / pR(I)
   enddo
   do J=1, Nx/2
     write(2211,*) pR(I),px(J), abs(psi_out_R(I,J))**2, abs(psi_out_Rx(I,J))**2
     write(2241,*) pR(I)*pR(I)/m_red, px(J), m_red*abs(psi_out_R(I,J))**2 / pR(I), &
             & m_red*abs(psi_out_Rx(I,J))**2 / pR(I)
   enddo
   write(2211,*)
   write(2241,*)
  enddo
  close(2211)
  close(2241)

  open(3311,file='Joint_elec_nucl_KER_spectra_diss.out',status='unknown')
  do I = 1, NR/2
    do J = 1, NR/2 
      JES_en_out_x(I,J) =(m_eff*m_red) * ((abs(psi_out_R(I,J))**2)/(pR(I)*px(J)) &
                 & + (abs(psi_out_R(I,NR-J))**2)/(pR(I)*px(NR-J)))
      JES_en_out_Rx(I,J) =(m_eff*m_red) * ((abs(psi_out_Rx(I,J))**2)/(pR(I)*px(J)) &
                 & + (abs(psi_out_Rx(I,NR-J))**2)/(pR(I)*px(NR-J)))
      write(3311,*) pR(I)*pR(I)/(2*m_red), px(J)*px(J)/2, &
                 & JES_en_out_x(I,J), JES_en_out_Rx(I,J)
    enddo
    write(3311,*)
  enddo
  close(3311)
!  
  momt_R=0.0d0 
  call density(psi_out_R, idensR, idensx)
  do I=1,NR
    momt_R=momt_R+dble(pR(I)*idensR(I))
  enddo
    momt_R=momt_R*dR
   ! momt_R=momt_R/norm_outR
  Print*, 'Nuclear momentum: H + H+ pathway'
  print*, 'average momentum=', sngl(momt_R), 'a.u.'
  print*, 'average velocity=', sngl((momt_R*au2a)/(mass*au2fs)), 'A°/fs'

  call density(psi_out_Rx, idensR, idensx)
  do I=1,NR
    momt_R=momt_R+dble(pR(I)*idensR(I))
  enddo
    momt_R=momt_R*dR
  Print*, 'Nuclear momentum: H+ + [H+ + e-] pathway'
  print*, 'average momentum=', sngl(momt_R), 'a.u.'
  print*, 'average velocity=', sngl((momt_R*au2a)/(mass*au2fs)), 'A°/fs'

  call density(psi_out_x, idensR, idensx)
  do I=1,NR
    momt_R=momt_R+dble(pR(I)*idensR(I))
  enddo
    momt_R=momt_R*dR
  Print*, 'Nuclear momentum: H_2^+ +  e- pathway'
  print*, 'average momentum=', sngl(momt_R), 'a.u.'

  call density(psi_out_xR, idensR, idensx)
  do I=1,NR
    momt_R=momt_R+dble(pR(I)*idensR(I))
  enddo
    momt_R=momt_R*dR
  Print*, 'Nuclear momentum: H+ + H+ + e- pathway'
  print*, 'average momentum=', sngl(momt_R), 'a.u.'

 deallocate(idensR, idensx, psi, kprop, psi0, idenspx, idenspR)
 deallocate(psi_out_x, psi_out_R, psi_out_R1, psi_out_x1, psi_out_x_momt) 
 deallocate(psi_out_Rx, psi_out_xR, psi_out_Rx1, psi_out_xR1)
 deallocate(JES_en_out_x, JES_en_out_xR, JES_en_out_Rx)
! deallocate(psi_out_Rl1, psi_out_Rr1)

! close(98, status='keep')
 close(100, status='keep')
 close(200, status='keep')
 close(201, status='keep')
 close(202, status='keep')
 close(203, status='keep')
 close(500, status='keep')
 close(501, status='keep')
 close(502, status='keep')
 close(503, status='keep')
! close(550, status='keep')
 close(600, status='keep')
 close(601, status='keep')
! close(800, status='keep')
! close(801, status='keep')
 close(802, status='keep')
! close(803, status='keep')
 close(909, status='keep')
! close(907, status='keep') 
! close(908, status='keep')


 call dfftw_destroy_plan(planF)
 call dfftw_destroy_plan(planB)
 call dfftw_destroy_plan(planFd)
 call dfftw_destroy_plan(planBd)
 call dfftw_destroy_plan(planFx)
 call dfftw_destroy_plan(planBx)
 call dfftw_destroy_plan(planFR)
 call dfftw_destroy_plan(planBR)

 

return
end subroutine

!_________________________________________________________

subroutine integ_2D(psi, norm)

use data_grid
 implicit none
 integer I, J

 double precision,intent(out):: norm
 complex*16,intent(in):: psi(NR, Nx)

 norm = 0.d0
 norm = sum(abs(psi(:,:))**2) * dx*dR
return
end subroutine


! ..................................................................

subroutine density(psi,idensR,idensx)

use data_grid
use pot_param,only:R0,x0

 implicit none
 integer:: I, J
 double precision,intent(out)::idensx(Nx), idensR(NR) 
 complex*16,intent(in):: psi(NR,Nx)
 
   idensR = 0.d0
   idensx = 0.d0

  do I = 1, NR
    idensR(I) = sum(abs(psi(I,:))**2)*dx
  end do


  do J = 1, Nx
    idensx(J) = sum(abs(psi(:,J))**2)*dR
  end do

return
end subroutine

!...................................................

subroutine overlap_2d(psi1, psi2, C)

use data_grid
 implicit none
 integer:: I, J
 complex*16:: C, F(NR)
 complex*16,intent(in):: psi1(NR,Nx), psi2(NR,Nx)

 C = (0.d0, 0.d0); F = (0.d0, 0.d0)

 do I = 1, NR
  do J = 1, Nx

   f(i) = f(i) + conjg(psi1(i,j)) * psi2(i,j)

  end do
 end do

  f = f * dx

  do i = 1, NR
   C = C + F(I)
  end do

  C = C * dR


return
end subroutine

!........................................................................

subroutine pop_analysis(psi, time, ewf, K)

use data_grid
use data_au
use pot_param

 implicit none

 integer:: I, J, K, N

 double precision,intent(in):: time, ewf(Nx, NR, Nstates)
 double precision:: B(Nstates),ax(2), cx(nr,2)
 complex(kind=kind(0.d0)),intent(in):: psi(NR,Nx)
 complex(kind=kind(0.d0)):: a(NR,Nstates),psi2(NR,Nx)
 complex(kind=kind(0.d0)):: sg(nr,nx), su(nr,nx)


  a = (0.d0,0.d0)
  b = 0.d0
  cx = (0.d0, 0.d0)
  ax = 0.d0

 psi2 = psi
 

  do N = 1, Nstates
   do I = 1, NR

    do J = 1, Nx
     a(i,n) = a(i,n) + ewf(j,i,n) * psi2(i,j)
    end do

   end do
  end do

  a = a * dx


  do N = 1, Nstates
   do I = 1, NR
     b(n) = b(n) + abs(a(i,n))**2
   end do
  end do

  b = b * dR


  write(500,*) sngl(time *au2fs), b(1)
  write(501,*) sngl(time *au2fs), sngl(b(2)), sngl(b(3)), sngl(b(4))


 if(mod(K,100).eq.0) then  ! writing out the probability amplitudes
   do I = 1, NR
    write(600,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(a(i,1))**2)
    write(601,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(a(i,2))**2)
   end do

   write(600,*)
   write(601,*)
 end if


   
return
end subroutine

! .......................................................................

subroutine osc_dipole(psi, d_t1, d_t2,grad)

use data_grid
use pot_param

 implicit none

 integer:: I, J
 double precision,intent(out):: d_t1, d_t2(NR)
 double precision,intent(in)::grad(nr,nx)
 complex*16,intent(in):: psi(NR,Nx)
 
 d_t1 = 0.d0
 d_t2 = 0.d0

 do i = 1, Nr
  do j = 1, Nx
   d_t2(i) = d_t2(i) + abs(psi(i,j))**2 * grad(i,j)    
  end do
   d_t2(i) = -d_t2(i) * dx
 end do 
 
 
 do i = 1, Nr
  d_t1 = d_t1 + d_t2(i)
 end do
 
  d_t1 = d_t1 * dR
  

return 
end subroutine

!........................................................................

subroutine localization(psi2, time, ewf, K)

use data_grid
use data_au
use pot_param

 implicit none

 integer:: I, J, K, N
 
 double precision time, ewf(Nx,NR,Nstates)
 double precision pl_loc(NR,Nx), neg_loc(NR,Nx)
 double precision b1(2)
 complex*16 amp_neg(NR), amp_pl(NR)
 complex*16 psi2(NR,Nx)


! do I = 1, NR
!  do J = 1, Nx
!   psi2(J,I) = psi(I,J)
!  enddo
! enddo
 
  
 pl_loc = 0.d0
 neg_loc = 0.d0
 
 do i = 1, NR
   pl_loc(I,:) = 1./sqrt(2.d0)* (ewf(:,I,1) + ewf(:,I,2))
   neg_loc(I,:) = 1./sqrt(2.d0)* (ewf(:,I,1) - ewf(:,I,2))
 end do


  amp_pl = 0.d0
  amp_neg = 0.d0
  b1 = 0.d0

   do I=1,NR 
    amp_pl(I) = sum(1.d0/sqrt(2.d0)* (ewf(:,I,1)+ewf(:,I,2)) * psi2(I,:))*dx
    amp_neg(I) = sum(1.d0/sqrt(2.d0)* (ewf(:,I,1)-ewf(:,I,2)) * psi2(I,:))*dx
   enddo

!  amp1 = amp1 * dx

   
  
   b1(1) = sum(abs(amp_pl(:))**2)*dR
   b1(2) = sum(abs(amp_neg(:))**2)*dR

  
  
 if(mod(K,100).eq.0) then  ! writing out the probability amplitudes
   do I = 1, NR
    write(504,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(amp_pl(I))**2)
    write(505,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(amp_neg(I))**2)
   end do

   write(504,*)
   write(505,*)
 end if

  write(502,*) sngl(time *au2fs), b1(1)
  write(503,*) sngl(time *au2fs), b1(2)
 


return
end subroutine


!............... Cut off Functions ................


subroutine cofs(cof_choice, cof,cofxl,cofxr, cofR)

use data_grid
use data_au
use pot_param

 implicit none
 integer:: i, j
 integer:: cof_choice 
 double precision:: cof(nx), cofr(nr)
 double precision:: cofxl(nx), cofxr(nx)

  cof = 0.d0
  cofR = 0.d0
!  cpm = 20.d0 / au2a
!  cpmr = 2.d0 / au2a

 select case (cof_choice)
!  call cutoff_cos_2d(cpm, cpmR, cof, cofR)
  case(1)
  call cutoff_ex_2d(cpm,cpmR,cof,cofR)
  case(2)
  call cutoff_ex_2d(cpmx,cpmR,cof,cofR)
 end select
  cofxr=cof
  do I=1,Nx
   cofxl(I)=cof(Nx-i+1)
  enddo

  do I = 1, Nx/2
   cof(i) = cof(nx-i+1)
  end do
 

!  do I = 1, Nx
!    x = x0 + (I - 1) * dx
!    write(98,*) I, sngl(x), sngl(cof(i))
!  end do

!  do j = 1, NR
!   R = R0 + (j - 1) * dR
!   write(99,*) J, sngl(R * au2a), sngl(cofR(j))
!  end do
! close(98, status='keep')
! close(99, status='keep')

 
return
end subroutine

!------------------------------------------------
subroutine cutoff_cos_2d(cpm, cpmR,cof, cofR)
use data_grid
use data_au
use pot_param, only: xend, Rend

implicit none

   integer :: I,J
   double precision:: cpmR
   double precision:: cof(Nx),cpm,cofR(NR)
   do I = 1, Nx
     if(abs(x(I)).lt.abs((xend - cpm))) then
     cof(i) = 1.d0
     else
     cof(i) = cos(((x(I) - xend + cpm) / -cpm) * (0.5d0 * pi))
     cof(i) = cof(i)**2
     end if
   end do


   do j = 1, NR
    if(R(J).lt.(Rend - cpmR)) then
    cofR(j) = 1.d0
    else
    cofR(j) = dcos(((R(J) - Rend + cpmR) / -cpmR) * (0.5d0 * pi))
    cofR(j) = cofR(j)**2
    end if
   end do

 return
 end subroutine
 !------------------------------------------------
 subroutine cutoff_ex_2d(cpm, cpmR,cof, cofR)
 use data_grid
 use data_au
 use pot_param, only: xend, Rend

  implicit none

    integer ::I, J
    double precision:: cof(Nx),cpm,c
    double precision:: cofR(NR), cpmR, cr
open(1,file='c.out')
  c=0.1500d0
  cr=1.00d0
  do I = 1, Nx
    cof(i) = 1.0d0/(1.0d0+exp(c*(x(I)-xend+cpm)))
 end do


 do j = 1, NR
   cofR(j)=1.0d0/(1.0d0+exp(cr*(R(J)-Rend+cpmR)))
   write(1,*) J, sngl(R(J)*au2a)
 end do
  write(1,*) J
  close(1)
 return
 end subroutine

 !----------------------------------------------------
 
 subroutine length_gauge_propagation(planF, planB, E)
 use data_grid
 use wavepacket_vars
 use data_au, only: im

 implicit none
 integer:: I, J
 integer*8:: planF, planB
 double precision:: E

!###############################################
!##### length gauge propagation ################
!###############################################
  do I =1,NR 
    psi(I,:) = psi(I,:) * exp(-0.5d0*im * dt * (pot(I,:)+x(:)*E))
  enddo
  call dfftw_execute(planF)
  psi = psi * kprop
  psi = psi / sqrt(dble(NR * Nx))
  call density(psi,idenspR,idenspx)
  call dfftw_execute(planB)
  psi = psi / sqrt(dble(NR * Nx))
  do I =1,NR 
    psi(I,:) = psi(I,:) * exp(-0.5d0*im * dt * (pot(I,:)+x(:)*E))
  enddo

end subroutine
