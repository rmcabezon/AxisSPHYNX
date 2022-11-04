PROGRAM AxisSPHYNX_MHD

      implicit none

      integer, parameter :: n1=131072
      integer, parameter :: n2=2*n1
      integer, parameter :: nnl=500
      integer, parameter :: lmax=10000000
      integer, parameter :: nc=8     ! number of boxes (normal+refective+periodic)
      integer, parameter :: ni0r=60
      integer, parameter :: vecMax=200
      integer, parameter :: nt=2*n2
      integer, parameter :: ntree=n2
      integer, parameter :: iodelay=nnl+50
      integer, parameter :: initial_adjust_iterations=15

      logical, parameter :: averVelr=.true.
      logical, parameter :: crossedMom=.true.
      logical, parameter :: timeaxis=.false.
      logical, parameter :: timeTemp=.false.
      logical, parameter :: balsara_lim=.true.
      logical, parameter :: iad0=.false.
      logical, parameter :: cleandiv=.true.
      logical, parameter :: profiling=.false.

      double precision, parameter :: pi=acos(-1.d0)
      double precision, parameter :: pim1=1.d0/pi
      double precision, parameter :: pi2=pi*0.5d0
      double precision, parameter :: year=365.d0*86400.d0
      double precision, parameter :: pc=3.086d18      ! parsec to cm
      double precision, parameter :: g=6.67d-08
      double precision, parameter :: msun=1.989d33


      double precision, parameter :: rgas=8.317d7
      double precision, parameter :: gamma=5.d0/3.d0     ! for the barotropic EOS, internal energy
      double precision, parameter :: gammam1=gamma-1.d0
      double precision, parameter :: gamma_bar=4.d0/3.d0    ! for the barotropic EOS, pressure and sound speed
      double precision, parameter :: gammam1_bar=gamma_bar-1.d0
      double precision, parameter :: cs0=2.d4

      !RESTART parameters
      integer, parameter :: lstart=1
      double precision, parameter :: initial_timestep=1.d-05
      double precision, parameter :: timeinit=0.d0
      character*21, parameter     :: inputfile='InitialZpinch_0.dat'
      !End RESTART parameters

      double precision, parameter :: ro3d_0=1.d-14
      double precision, parameter :: Cour=0.4d0
      double precision, parameter :: C_accel=0.25d0
      double precision, parameter :: tol_temp=0.05d0
      double precision, parameter :: tol_taxis=0.05
      double precision, parameter :: alfaAV=1.d0
      double precision, parameter :: betaAV=2.d0
      double precision, parameter :: alfau=0.05
      double precision, parameter :: alfaB=0.5d0
      double precision, parameter :: hindex=5.d0
      double precision, parameter :: epsdiv=0.2d0
      double precision, parameter :: mass_scaling=1.d0
      double precision, parameter :: rboundary_scaling=1.d0     ! fine tuning of bounday condition
      double precision, parameter :: minvalue=1.d-60



      !double precision, parameter :: mu0=4.d0*pi    ! Magnetic permeability. In IS: 4*pi*1e-7. In cgs:  4*pi (???)
      double precision, parameter :: mu0=1.d0       ! Magnetic permeability in code units. B(cgs)= B(code) x sqrt(4 pi)
      double precision, parameter :: fclean=1.d0    ! To clean divergence. Better <=1
      double precision, parameter :: sigma_c=1.d0   ! Adimensional constant to control the magnetic  decay term during cleaning

      integer signature,value1,signatureC,signatureNC,index,irhomax,irhomin
      integer nvmin,nvmax,nvr,nsatura,nje,nm,npp,nnm,i1,i2,npp0,nmt
      integer invmax,nvSum,npointer1,ntotal,nmmm,ncount,nodo
      integer nnma,l,kl,i,invmin,j,jg,jVec,k,n,nav


      integer, dimension(n1) :: nv,nvvv
      integer, dimension(n2) :: ngn,ngh
      integer, dimension(nt) :: nnod,npa,np,npn,npart,npointer
      integer, dimension(0:nt) :: nh1,nh2
      integer, dimension(ntree) :: naux,nauxa
      integer, dimension(n1,vecMax) :: vecMatrix

      double precision m1,nu2,nindex1,nindex2,pkv,pmt,ro3dij,ro3dmax,ro3dmin
      double precision momentumxx,momentumyy,momentumzz,momentumxxp,momentumyyp,romax
      double precision kernxxi,kernxyi,kernyxi,kernyyi,rr,rrr,ux,uy,v1,val,rhomax
      double precision kernxxj,kernxyj,kernyxj,kernyyj,kernxx,kernyy,kx
      double precision iadxx,iadyy,iadzz,accel,angmomphi,angmomtot,angmomr,angmomz
      double precision lambda,masstot,lcentered,ax,ay,Bcx,Bcz,Bcphi,Bmo12,Bmod
      double precision dtnma,dtnme,tt,dt,xt,dx,dy,h1,d02,d05,pi2v1,w1d,dter
      double precision Bmod2,txx,txy,tyy,hhi,d11,d21,rij_2,rdenom,rinvdenom
      double precision ddmax,ddist,rad,h0,h11,hexp,hfac
      double precision divi,divj,roti,rotj,fbalsi,fbalsj,fij,qiji,qijj,qijiu,qijju,tqij
      double precision proi,proj,proeneri,tp,ti3,fix0,djix,djiy,vjix,vjiy,vjiz
      double precision gradxha,gradyha,gradxhb,gradyhb,Bphiav,Bradav,Bzav,checkdeltar
      double precision dteri,dvx,dvy,dmy,dmy1,dmy2,smdiv,tau_a,beta,hhj,h2,uijx,uijy
      double precision v2,pi2v2,w2d,dterj,dj,vijrij,r3dij,rij,wij,vijsignal,vijsignalu
      double precision dmy01,dmy02,dmy11,dmy12,cm,cm0,d1,d2,dc,di,dtn,du,gradxhab,gradyhab
      double precision aviscox,aviscoy,dpx,dpy,dpxavu,dpyavu,Bij2,Fabha,Fabhb
      double precision vijsignalBi,vijsignalBj,dmyr,dmyz,dmyphy,vijsignaldis
      double precision viscBdis,viscBi,viscBj,dmyB,dmyBclean,Bdist,Btot2,hoop_r,hoop_phi
      double precision hoop,ratio,epsilon,a11,a12,a13,a21,a22,a23,a31,a32,a33
      double precision smdi,egrav,eint,ekin,elost,emag,errdivB,errdivBav,ycma
      double precision time_accel,time_temp,time_Cour,time_axis,errdivBmax,errdivBmin,etot0
      double precision totalmassC,space,rdmax,zdmaxneg,zdmaxpos,rtr
      double precision rtzu,rtzd,aux1,aux2,aux3,dpdro,desplr,desplz,rtr5,rtz5,xcm,ycm,dmass
      double precision Sij_i11,Sij_i12,Sij_i13,Sij_i21,Sij_i22,Sij_i23,Sij_i31,Sij_i32,Sij_i33
      double precision Sij_j11,Sij_j12,Sij_j13,Sij_j21,Sij_j22,Sij_j23,Sij_j31,Sij_j32,Sij_j33
      double precision vij1,vij2,vij3,vmod,vr,x1,xcma,y1,ycm_check,angmomtot_0

      double precision, dimension(4)    :: ccf0
      double precision, dimension(n1)   :: sum3,evisco,eviscu,maxvsignal
      double precision, dimension(n1)   :: divB,cha0,cha
      double precision, dimension(n1)   :: checkInorm,checkdeltax,checkdeltay
      double precision, dimension(n1)   :: divvphir,divvphiz,divvrr,divvrz,divvzr,divvzz
      double precision, dimension(n1)   :: gradrBphi,gradrBr,gradrBz
      double precision, dimension(n1)   :: proener,sum1B,sum2B,sum3B,B2mu0,phi,eviscoB
      double precision, dimension(n1)   :: ap1,ap2,ap1_2,ap2_2,ri_2,ap1p2,ap2p2
      double precision, dimension(n1)   :: lz,lz0
      double precision, dimension(n1)   :: avevel,bfieldr,gradPr

      double precision, dimension(n2)   :: x,y,dudt,dudv,dpdt,h,mass,massa
      double precision, dimension(n2)   :: ro2d,ro3d,c,p,pro,u,u0
      double precision, dimension(n2)   :: sum1,sum2,sum5,va_sroot
      double precision, dimension(n2)   :: cxx,cxy,cyy
      double precision, dimension(n2)   :: divv,curlz,vol,valfven2,phidiv,phidiv0

      double precision, dimension(nt)   :: dxx,xcet,ycet,xce,yce

      double precision, dimension(n2,2) :: a,a0,f
      double precision, dimension(n1,3) :: Bdivclean,Bindis
      double precision, dimension(n2,3) :: v0,v,B,B0

      double precision, dimension(ntree,1)   :: xia,yia
      double precision, dimension(ntree,2)   :: xi,xd,yi,yd

      double precision, dimension(n2,3,3)  :: Sij

      logical, dimension(nt) :: particle
      logical, dimension(n2) :: axis

      character*14 nomfs1,nomfs2
      character*6  prefix,suffix


!     Interpolating kernels are the 2d-sinc kernels with index n (default n=5: i.e. hindex)
      data (ccf0(i),i=1,4)/1.3090244919d-01,1.935848488354d-02,-6.164290621127d-03,5.22450269233d-02/

      pkv=ccf0(1)*hindex+ccf0(2)/hindex+ccf0(3)/hindex**2+ccf0(4)
      dtnma=initial_timestep
      dtnme=dtnma
      dtn=(dtnma+dtnme)*0.5d0
      tt=timeinit

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      open(1,file='InitialZpinch_0.dat')
      open(2,file='tree_check_2.dat')
      open(3,file='timectrl_2.dat')

      v  = 0.d0
      lz  = 0.d0
      lz0 = lz
      phidiv  = 0.d0
      phidiv0 = phidiv
      sum1 = 0.d0
      sum2 = 0.d0
      sum3 = 0.d0


      if(profiling) call profile(0,'')
      !read data
      do i=1,n1
         read(1,*) a(i,1),a(i,2),dmy,mass(i),h(i),dmy,dmy,ro3d(i),&
              & ro2d(i),dmy,dmy,dmy,u(i),p(i),dmy,dmy, &
              &   dmy,dmy,dmy,dmy,dmy,dmy,dmy,nv(i)
      enddo
      if(profiling) call profile(1,'   read_data')

      ! Sets initial values of the magnetic field.
      B  = 0.d0
      B0 = B
      
      !Initial conditions for the Z-pinch test
      Bmod=3.d0
      do i=1,n1
         rr=sqrt(a(i,1)**2)
         rrr=exp(-(rr-0.5d0)**2/0.01d0)
         v(i,1)=0.d0
         v(i,2)=-1.d0
         v(i,3)=0.d0
         B(i,1)=0.d0
         B(i,2)=0.d0
         B(i,3)=Bmod*rrr
         p(i)=1.d0
         u(i)=p(i)/ro3d(i)/(gamma-1.d0)
         u0(i)=u(i)
         lz0(i)=v(i,3)*a(i,1)
         lz(i)=lz0(i)
         v0(i,1)=v(i,1)
         v0(i,2)=v(i,2)
         v0(i,3)=v(i,3)
         B0(i,1)=B(i,1)
         B0(i,2)=B(i,2)
         B0(i,3)=B(i,3)
      enddo

      !Interparticle distance.
      !XXX Only valid for grid. To be review for a glass ICs
      space=abs(a(1,2)-a(2,2))/2.d0

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !Maximum distance in r (radial)  and z (vertical)
      rdmax=maxval(abs(a(:,1)))
      zdmaxpos=maxval(a(:,2))
      zdmaxneg=abs(minval(a(:,2)))
      
      rtr=rdmax+space*rboundary_scaling
      rtzu=abs(zdmaxpos)+space
      rtzd=abs(zdmaxneg)+space
      
      !Initial EOS call (id. gas)
      do  i=1,n1
         p(i)=u(i)*(gamma-1.d0)*ro3d(i)
         c(i)=dsqrt(p(i)/ro3d(i))
         u0(i)=u(i)
      enddo
      
      !Ghosts particles
      desplr=0.d0
      desplz=0.d0
      rtr5=rtr/5.d0         !Copy just a fifth of the domain to save memory
      rtz5=rtzu/5.d0
      
      do i=1,n2
         ngh(i)=i
         axis(i)=.false.
      enddo
      
      do i=1,nc
         call ghosts(desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5,i,n,&
              &  a,v,B,n1,n2,axis,ngh)
      enddo
      desplr=3.d0 
      desplz=3.d0 
      a(:,1)=a(:,1)+desplr
      a(:,2)=a(:,2)+desplz
      a0=a
      v0=v
      do i=1,n
         h(i)=h(ngh(i))
         u(i)=u(ngh(i))
         mass(i)=mass(ngh(i))
         ro3d(i)=ro3d(ngh(i))
      enddo
      u0=u
      !*** CENTER OF MASS ***
      ycm=0.d0
      dmass=0.d0
      do i=1,n1
         ycm=ycm+mass(i)*a(i,2)
         dmass=dmass+mass(i)
      enddo
      !ycm = ycm/dmass
      xcm=desplr
      ycm=desplz
      
      !*** SIZE OF THE SYSTEM ***
      ddmax = 0.0d0
      do i = 1,n
         ddist=sqrt((xcm-a(i,1))**2+(ycm-a(i,2))**2)
         ddmax=max(ddmax,ddist)
      end do
      rad = ddmax
      
      
      !     **************************************************************
      !     Iterations start, nnl is the total number of iterations
      !     kl=1,2 is centering the schemei
      do l=lstart,nnl
         do kl=1,2
            print *,l,kl
            dt = (dtnma / 2.d0) * dble(kl)
            if(kl.eq.1) write(3,'(1x,i5,7(1x,e14.7))') l,time_accel,time_temp,time_Cour,time_axis,dtnma,tt
            
            ! TREE CALCULATION
            if(profiling) call profile(0,'')
            do i = 1,n
               massa(i)=mass(i)
               x(i)=a(i,1)
               y(i)=a(i,2)
               npointer(i)=i
            end do
            do i=1,4
               nauxa(i)=0
            end do
            nnma=0
            ntotal=0
            npn(1)=1
            npn(2)=n
            
            xt=pi*sqrt(xcm**2+ycm**2)/sqrt(7.d0)
            nje=0
            nm=1
            xia(1,1)=0.d0
            yia(1,1)=0.d0
6           nje=nje+1
            npp=0
            dx=xt/2.d0**nje
            dy=dx
            do j=1,nm
               np(2*j-1)=npn(2*j-1)
               np(2*j)=npn(2*j)
               xi(j,1)=xia(j,1)
               yi(j,1)=yia(j,1)
               naux(j)=nauxa(j)
            enddo
            nnm=nm
            if(nnm .ge. nnma) nnma=nnm
            if(nnm .gt. ntree-10) nsatura=1
            nm=0
            do j=1,nnm
               nmmm=0
               do i1=1,2
                  xi(j,i1)=xi(j,1)+(i1-1)*dx
                  xd(j,i1)=xi(j,i1)+dx
                  do i2=1,2
                     yi(j,i2)=yi(j,1)+(i2-1)*dy
                     yd(j,i2)=yi(j,i2)+dy
                     npp0=0
                     xcma=0.d0
                     ycma=0.d0
                     pmt=0.d0
                     do i=np(2*j-1),np(2*j)
                        if((xi(j,i1) .le. x(i)) .and. (x(i).lt.xd(j,i1))) then
                           if((yi(j,i2) .le. y(i)) .and. (y(i).lt.yd(j,i2))) then
                              npp=npp+1
                              npp0=npp0+1
                              xcma=xcma+x(i)*massa(i)
                              ycma=ycma+y(i)*massa(i)
                              pmt=pmt+massa(i)
                              x1=x(npp)
                              y1=y(npp)
                              m1=massa(npp)
                              npointer1=npointer(npp)
                              x(npp)=x(i)
                              y(npp)=y(i)
                              massa(npp)=massa(i)
                              npointer(npp)=npointer(i)
                              x(i)=x1
                              y(i)=y1
                              massa(i)=m1
                              npointer(i)=npointer1
                           endif
                        endif
                     enddo
                     if(npp0-1.lt.0) goto 25
                     if(npp0-1.eq.0) goto 26
27                   ntotal=ntotal+1
                     if(ntotal .gt. nt-50) nsatura=1
                     nm=nm+1
                     nmmm=nmmm+1
                     xcet(ntotal)=(xi(j,i1)+xd(j,i1))/2.d0
                     ycet(ntotal)=(yi(j,i2)+yd(j,i2))/2.d0
                     npa(ntotal)=naux(j)
                     nauxa(nm)=ntotal
                     dxx(ntotal)=dx
                     particle(ntotal)=.false.
                     npn(2*nm-1)=npp-npp0+1
                     npn(2*nm)=npp
                     xia(nm,1)=xi(j,i1)
                     yia(nm,1)=yi(j,i2)
                     goto 28
26                   ntotal=ntotal+1
                     if(ntotal .gt. nt-50) nsatura=1
                     npa(ntotal)=naux(j)
                     nmmm=nmmm+1
                     xcet(ntotal)=x(npp)
                     ycet(ntotal)=y(npp)
                     dxx(ntotal)=0.d0
                     particle(ntotal)=.true.
                     nnod(ntotal)=npointer(npp)
                     npart(npointer(npp))=ntotal
                     
25                   continue
28                   continue
                  end do
               enddo
               nh2(naux(j))=ntotal
               nh1(naux(j))=ntotal-nmmm+1
            end do
            if(nm .ne. 0) goto 6
            
            !Check that all particles have been found
            ncount=0
            do i=1,ntotal
               if(particle(i)) ncount=ncount+1
            enddo
            write(2,'(2x,6(2x,i6))') l,kl,n,ncount,ntotal,nsatura
            
            if(profiling) call profile(1,'   tree_calc')
            
            ! *** initial setting of h with ni0r.
            if(kl.eq.1) then
               do i=1,n1
                  if(nv(i) .le. 20) then
                     nv(i)=20
                     hfac=15.d0
                     hexp=1.d0/4.d0
                  else if(nv(i) .ge. 150) then
                     nv(i)=150
                     hfac=31.d0
                     hexp=1.d0/5.d0
                  else
                     hfac=31.d0
                     hexp=1.d0/5.d0
                  endif
                  h0=h(i)
                  h(i)=h(i)*0.5d0*(1.d0+hfac*dble(ni0r)/dble(nv(i)))**hexp
                  h(i)=0.5d0*(h0+h(i))
               enddo
            endif

            if(l.le.initial_adjust_iterations) then
               do i=1,n1
                  p(i)=1.d0
                  u(i)=p(i)/(gamma-1.d0)/ro3d(i)
                  u0(i)=u(i)
               enddo
            endif

            if(profiling) call profile(0,'')
            if(ncount.eq.n .or. l.le.100) then
               !!BUILDING THE LISTED-LINK of NEIGHBORS
               !!Real and axis-reflective ghosts particles
               nv(1:n1)=0

               !$omp parallel private(i,nvSum,nodo,d1,d2,dc,j,d05,v1,v2)
               !$omp do schedule(guided)
               do i = 1,n1
412               nvSum=0
                  nodo=0
411               nodo=nodo+1
                  d1=abs(a(i,1)-xcet(nodo))
                  d2=abs(a(i,2)-ycet(nodo))
                  dc=4.d0*h(i)+dxx(nodo)/2.d0
                  !test intersection of hierarchical cube with the action
                  !cube of particle i

                  if((d1.lt.dc).and.(d2.lt.dc))goto 451
                  
421               if(nodo.eq.nh2(npa(nodo))) then
                     nodo=npa(nodo)
                     if(nodo.eq.0) goto 4011
                     goto 421
                  else
                     goto 411
                  endif

451               if(particle(nodo)) then
                     j=nnod(nodo)
                     if(i.eq.j) goto 421
                     d05=sqrt(d1**2+d2**2)
                     v1=d05/h(i)
                     v2=d05/h(j)
                     if((v1.ne.0.d0.and.v1.lt.2.d0).or.(v2.lt.2.d0)) then
                     !if(v1.ne.0.d0.and.v1.lt.2.d0) then
                        nvSum=nvSum+1
                        if(nvSum.le.(vecMax-20)) then
                           vecMatrix(i,nvSum)=j
                        else
!                          print *,"Too many neighbors",nvSum,i
!                          stop
                        endif
                     endif
                     goto 421
                  endif
                  nodo=nh1(nodo)-1
                  goto 411
4011              continue
                  nv(i)=min(nvSum,vecMax)
                  !if number of neighbors lesser than 10 change hi
                  if(nv(i) .le. 10) then
                     nv(i)=0
                     h(i)=1.2d0*h(i)
                     goto 412
                  endif
                  !if number of neighbors is larger than 120 change hi
                  if(nv(i) .ge. 120) then
                     nv(i)=0
                     h(i)=0.9d0*h(i)
                     goto 412
                  endif
               end do
               !$omp end do
               !$omp end parallel
            endif      ! close if(ncount.eq.n)

            if(profiling) call profile(1,'  find_neigh')

            if(profiling) call profile(0,'')
            nvmax=maxval(nv)
            nvmin=minval(nv)
            invmax=maxloc(nv,DIM=1)
            invmin=minloc(nv,DIM=1)
            ro3dmax=maxval(ro3d)
            ro3dmin=minval(ro3d)
            irhomax=maxloc(ro3d,DIM=1)
            irhomin=minloc(ro3d,DIM=1)
            
            !Density calculation
            ro2d=0.d0
            !$omp parallel private(i,h1,jVec,j,d1,d2,d02,d05,v1,pi2v1,w1d,dter,signature)
            !$omp do
            do i = 1,n1
               h1=1.d0/h(i)**2
               do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  d1=a(i,1)-a(j,1)
                  d2=a(i,2)-a(j,2)
                  d02=d1**2+d2**2
                  d05=sqrt(d02)
                  v1=d05/h(i)
                  pi2v1=pi2*v1
                  if(v1.lt.2.d0) then
                     call func(pi2v1,w1d)
                  else
                     w1d=0.d0
                  endif
                  w1d=w1d**hindex
                  dter=h1*w1d
                  signature=1
                  if(axis(j)) signature=-1
                  ro2d(i)=ro2d(i)+mass(j)*dter*signature
               enddo
               !Self-contribution and VEs
               ro2d(i)=pkv*ro2d(i)+pkv*mass(i)*h1
               ro3d(i)=ro2d(i)/(2.d0*pi*(a(i,1)-desplr))
               vol(i)=mass(i)/ro2d(i)
            end do
            !$omp end do
            !$omp end parallel

            if(profiling) call profile(1,'density_calc')
            !EOS call (id. gas)
            do i=1,n1
               p(i)=u(i)*(gamma-1.)*ro3d(i)
               c(i)=dsqrt(p(i)/ro3d(i))
            enddo

            !Alfven velocity
            do i=1,n1
               valfven2(i)=(B(i,1)**2+B(i,2)**2+B(i,3)**2)/mu0/ro3d(i)
               va_sroot(i)=dsqrt(c(i)**2+valfven2(i))
            enddo

            !calculate the Stress tensor Sij  (including Magnetic field)
            if(profiling) call profile(0,'')
            !$omp parallel private(k,Bmod2,i,j)
            !$omp do
            do k=1,n1
               Bmod2=(B(k,1)**2+B(k,2)**2+B(k,3)**2)
               B2mu0(k)=Bmod2/(2.d0*mu0)
               pro(k)=(p(k)+B2mu0(k))/ro2d(k)
               proener(k)=p(k)/ro2d(k)
               do i=1,3
                  do j=1,3
                     Sij(k,i,j)=B(k,i)*B(k,j)/mu0
                     Sij(k,i,j)=Sij(k,i,j)/ro2d(k)
                  enddo
               enddo
            enddo
            !$omp end do

            !re-assign magnitudes (real+ ghosts).

            !$omp do
            do i=1,n
               h(i)=h(ngh(i))
               ro2d(i)=ro2d(ngh(i))
               ro3d(i)=ro3d(ngh(i))
               pro(i)=pro(ngh(i))
               c(i)=c(ngh(i))
               u(i)=u(ngh(i))
               p(i)=p(ngh(i))
               vol(i)=vol(ngh(i))
               valfven2(i)=valfven2(ngh(i))
               va_sroot(i)=va_sroot(ngh(i))
               do j=1,3
                  do k=1,3
                     if(i.eq.2 .and. j.eq.3 .and. axis(j)) then
                        Sij(i,j,k)=-Sij(ngh(i),j,k)
                     else
                        Sij(i,j,k)=Sij(ngh(i),j,k)
                     endif
                  enddo
               enddo
            enddo
            !$omp end do
            !$omp end parallel

            if(profiling) call profile(1,' Stress_tens')

            !IAD matrix calculation
            if(profiling) call profile(0,'')
            checkInorm=0.d0
            checkdeltax=0.d0
            checkdeltay=0.d0
            
            !$omp parallel private(i,hhi,h1,txx,txy,tyy,jVec,j,d1,d2,d02,d05,&
            !$omp    v1,pi2v1,w1d,dter,val)
            !$omp do
            do i = 1,n1
               hhi=h(i)*h(i)
               h1=1.d0/hhi
               txx=0.d0
               tyy=0.d0
               txy=0.d0
               do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  d1=abs(a(i,1)-a(j,1))
                  d2=abs(a(i,2)-a(j,2))
                  d02=d1**2+d2**2
                  d05=sqrt(d02)
                  v1=d05/h(i)
                  pi2v1=pi2*v1
                  if(v1.lt.2.d0) then
                     call func(pi2v1,w1d)
                  else
                     w1d=0.d0
                  endif
                  w1d=pkv*w1d**hindex
                  
                  !calculation of  matrix Tau
                  dter=h1*w1d
                  val=vol(j)*dter
                  txx=txx+val*(a(j,1)-a(i,1))**2
                  tyy=tyy+val*(a(j,2)-a(i,2))**2
                  txy=txy+val*(a(j,1)-a(i,1))*(a(j,2)-a(i,2))
                  checkInorm(i)=checkInorm(i)+val
                  checkdeltax(i)=checkdeltax(i)+ (a(j,1)-a(i,1))*val
                  checkdeltay(i)= checkdeltay(i)+(a(j,2)-a(i,2))*val
               enddo

               !Matrix inversion
               cxx(i)=txx-txy**2/tyy
               cxx(i)=1.d0/cxx(i)
               cxy(i)=-cxx(i)*txy/tyy
               cyy(i)=tyy-txy**2/txx
               cyy(i)=1.d0/cyy(i)
               
               !Self-contribution to checkInorm
               checkInorm(i)=checkInorm(i)+mass(i)/ro2d(i)/h(i)**2*pkv
            end do
            !$omp end do
            !$omp end parallel

            !Matrix elements (including ghosts)
            do i=1,n
               cxx(i)=cxx(ngh(i))
               cxy(i)=cxy(ngh(i))
               cyy(i)=cyy(ngh(i))
            enddo
            

            !     Summations to find the viscous stress tensor and divergence and curl
            !     of the velocity
            
            if(profiling) call profile(0,'')
            divv=0.d0
            curlz=0.d0
            divB=0.d0
            divvrr=0.d0
            divvrz=0.d0
            divvzr=0.d0
            divvzz=0.d0
            divvphir=0.d0
            divvphiz=0.d0
            gradrBr=0.d0
            gradrBphi=0.d0
            gradrBz=0.d0
            f=0.d0
            
            !$omp parallel private(i,hhi,h1,h11,jVec,j,djix,djiy,d1,d2,d02,d05,&
            !$omp   v1,pi2v1,w1d,dteri,vjix,vjiy,vjiz,kernxxi,kernxyi,kernyxi,kernyyi,&
            !$omp   gradxha,gradyha,signature,dvx,dvy,dmy,dmy1,dmy2,smdiv)
            !$omp do
            do i = 1,n1
               hhi=h(i)*h(i)
               h1=1.d0/hhi
               h11=h1/hhi
               do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  djix=a(j,1)-a(i,1)
                  djiy=a(j,2)-a(i,2)
                  d1=abs(djix)
                  d2=abs(djiy)
                  d02=d1**2+d2**2
                  d05=sqrt(d02)

                  !Divergence and curl calculation
                  v1=d05/h(i)
                  pi2v1=pi2*v1
                  w1d=0.d0
                  if(v1.lt.2.d0) then
                     call func(pi2v1,w1d)
                  else
                     w1d=0.d0
                  endif
                  w1d=pkv*w1d**hindex
                  dteri=h1*w1d
                  vjix=v(j,1)-v(i,1)
                  vjiy=v(j,2)-v(i,2)
                  vjiz=v(j,3)-v(i,3)
                  kernxxi=cxx(i)*djix
                  kernxyi=cxy(i)*djiy
                  kernyyi=cyy(i)*djiy
                  kernyxi=cxy(i)*djix
                  gradxha=(kernxxi+kernxyi)*dteri
                  gradyha=(kernyxi+kernyyi)*dteri
                  signature=1
                  if(axis(j)) signature=-1
                  dvx=vjix*gradxha
                  dvy=vjiy*gradyha
                  dmy=vol(j)*ro2d(j)*signature
                  divv(i)=divv(i)+(dvx+dvy)*dmy
                  dmy1=gradxha*dmy
                  dmy2=gradyha*dmy
                  curlz(i)=curlz(i)+(vjiy*gradxha-vjix*gradyha)*dmy
                  divvrr(i)=divvrr(i)+vjix*dmy1
                  divvrz(i)=divvrz(i)+vjix*dmy2
                  divvzr(i)=divvzr(i)+vjiy*dmy1
                  divvzz(i)=divvzz(i)+vjiy*dmy2
                  divvphir(i)=divvphir(i)+vjiz*dmy1
                  divvphiz(i)=divvphiz(i)+vjiz*dmy2
                  divB(i)=divB(i)+(B(j,1)-B(i,1))*dmy1+(B(j,2)-B(i,2))*dmy2
                  gradrBr(i)=gradrBr(i)+(B(j,1)-B(i,1))*dmy1
                  gradrBz(i)=gradrBz(i)+(B(j,2)-B(i,2))*dmy1
                  gradrBphi(i)=gradrBphi(i)+(B(j,3)-B(i,3))*dmy1
               enddo
               !Divide by ro2(i)
               divv(i)=divv(i)/ro2d(i)
               divvrr(i)=divvrr(i)/ro2d(i)
               divvrz(i)=divvrz(i)/ro2d(i)
               divvzr(i)=divvzr(i)/ro2d(i)
               divvzz(i)=divvzz(i)/ro2d(i)
               divvphir(i)=divvphir(i)/ro2d(i)
               divvphiz(i)=divvphiz(i)/ro2d(i)
               divB(i)=divB(i)/ro2d(i)
               curlz(i)=curlz(i)/ro2d(i)
               gradrBr(i)=gradrBr(i)/ro2d(i)
               gradrBz(i)=gradrBz(i)/ro2d(i)
               gradrBphi(i)=gradrBphi(i)/ro2d(i)
               !Add Hoop-Stress terms to divB,divv  (smoothed close to the axis)
               smdiv=a(i,1)-desplr+epsdiv*h(i)
               divv(i)=divv(i)+v(i,1)/smdiv
               divB(i)=divB(i)+B(i,1)/smdiv
            end do
            !$omp end do
            !$omp end parallel
            
            !set values for ghosts
            do i=1,n
               divv(i)=divv(ngh(i))
               curlz(i)=curlz(ngh(i))
            enddo

            if(profiling) call profile(1,' Visc_tensor')
            !To clean or not to clean the divergence. That is the question ..
            
            if(cleandiv) then
               !Calculation of the derivative of parameter phidiv necessary to
               !clean the divergence of B with the recipe of Wissing et al (2020).
               if(profiling) call profile(0,'')
               
               !$omp parallel private(i,tau_a)
               !$omp do
               do i=1,n1
                  cha(i)=fclean*va_sroot(i)
                  if(l.eq.lstart) cha0(i)=cha(i)
                  tau_a=h(i)/cha0(i)/sigma_c
                  phidiv(i)=(phidiv0(i)/cha0(i))*cha(i)-(cha(i)*cha0(i)*divB(i)+&
                       &    phidiv0(i)/2.d0*divv(i)*cha(i)/cha0(i)+&
                       &    phidiv0(i)/tau_a*cha(i)/cha0(i))*dt
               enddo
               !$omp end do
               !$omp end parallel
               
               if(profiling) call profile(1,'   Clean_div')
            endif

            !set values for ghosts
            do i=1,n
               phidiv(i)=phidiv(ngh(i))
            enddo

!     *********************************************************
!     calculation of the different terms in the momentum equation
!     **********************************************************

            if(profiling) call profile(0,'')
            sum1=0.d0; sum2=0.d0; sum3=0.d0;
            sum1B=0.d0; sum2B=0.d0; sum3B=0.d0
            sum5=0.d0
            Bindis=0.d0; Bdivclean=0.d0
            maxvsignal=0.d0;
            evisco=0.d0; eviscoB=0.d0; eviscu=0.d0
            

            !$omp parallel private(i,beta,hhi,h1,jVec,j,dmy1,dmy2,d1,d2,d02,hhj,&
            !$omp     h2,d05,uijx,uijy,v1,v2,pi2v1,pi2v2,w1d,w2d,djix,djiy,dteri,dterj,&
            !$omp     signatureC,signatureNC,kernxxi,kernxyi,kernyxi,kernyyi,kernxxj,&
            !$omp     kernxyj,kernyxj,kernyyj,gradxha,gradxhb,gradyha,gradyhb,di,dj,&
            !$omp     vij1,vij2,vij3,vijrij,rij,fij,wij,vijsignal,vijsignalu,ro3dij,&
            !$omp     divi,divj,roti,rotj,fbalsi,fbalsj,qiji,qijiu,qijj,qijju,proi,proj,&
            !$omp     Sij_i11,Sij_i12,Sij_i13,Sij_i21,Sij_i22,Sij_i23,Sij_i31,Sij_i32,Sij_i33,&
            !$omp     Sij_j11,Sij_j12,Sij_j13,Sij_j21,Sij_j22,Sij_j23,Sij_j31,Sij_j32,Sij_j33,&
            !$omp     proeneri,momentumxxp,momentumyyp,momentumxx,momentumyy,momentumzz,&
            !$omp     dmy01,dmy02,dmy11,dmy12,aviscox,aviscoy,dpxavu,dpyavu,iadxx,iadyy,iadzz,&
            !$omp     Bij2,Fabha,Fabhb,vijsignalBi,vijsignalBj,vijsignaldis,viscBi,viscBj,&
            !$omp     viscBdis,dmy,dmyB,dmyBclean,tqij,smdiv,Bdist,gradxhab,gradyhab)
            !$omp do
            do i = 1,n1
               if(B2mu0(i).ne.0.d0) then
                  beta=p(i)/B2mu0(i)
                  phi(i)=0.d0
                  if(beta.lt.1.d0) phi(i)=2.d0
                  if(1.d0.le.beta .and. beta.le.2.d0) phi(i)=2.d0*(2.d0-beta)
               else
                  phi(i)=0.d0
               endif
               hhi=h(i)*h(i)
               h1=1.d0/hhi
               do jVec=1,nv(i)
                  j=vecMatrix(i,jVec)
                  dmy1=a(i,1)-a(j,1)
                  dmy2=a(i,2)-a(j,2)
                  d1=abs(dmy1)
                  d2=abs(dmy2)
                  d02=d1**2+d2**2
                  hhj=h(j)*h(j)
                  h2=1.d0/hhj
                  d05=sqrt(d02)
                  uijx=dmy1/d05
                  uijy=dmy2/d05
                  v1=d05/h(i)
                  v2=d05/h(j)
                  pi2v1=pi2*v1
                  pi2v2=pi2*v2
                  if(v1.lt.2.d0) then
                     call func(pi2v1,w1d)
                  else
                     w1d=0.d0
                  endif
                  if(v2.lt.2.d0) then
                     call func(pi2v2,w2d)
                  else
                     w2d=0.d0
                  endif

                  !(j)-(i) is ihat is used in IAD gradient formula.
                  djix=-dmy1
                  djiy=-dmy2
                  w1d=pkv*w1d**hindex
                  w2d=pkv*w2d**hindex
                  dteri=h1*w1d
                  dterj=h2*w2d
                  signatureC=1
                  signatureNC=1
                  if(axis(j)) then
                     if(crossedMom) then
                        signatureC=-1
                     else
                        signatureNC=-1
                     endif
                  endif
                  kernxxi=cxx(i)*djix*signatureNC
                  kernxyi=cxy(i)*djiy*signatureNC
                  kernxxj=cxx(j)*djix*signatureC
                  kernxyj=cxy(j)*djiy*signatureC
                  kernyyi=cyy(i)*djiy*signatureNC
                  kernyxi=cxy(i)*djix*signatureNC
                  kernyyj=cyy(j)*djiy*signatureC
                  kernyxj=cxy(j)*djix*signatureC
                  gradxha=(kernxxi+kernxyi)*dteri
                  gradxhb=(kernxxj+kernxyj)*dterj
                  gradyha=(kernyxi+kernyyi)*dteri
                  gradyhb=(kernyxj+kernyyj)*dterj
                  gradxhab=0.5*(gradxha+gradxhb)
                  gradyhab=0.5*(gradyha+gradyhb)

                  ! Artificial viscosity calculation. Includes Balsara limiters.

                  di=a(i,1)-desplr
                  dj=a(j,1)-desplr
                  vij1=v(i,1)-v(j,1)
                  vij2=v(i,2)-v(j,2)
                  vij3=v(i,3)-v(j,3)
                  vijrij=vij1*(a(i,1)-a(j,1))+vij2*(a(i,2)-a(j,2))
                  rij=sqrt((a(j,1)-a(i,1))**2+(a(j,2)-a(i,2))**2)
                  fij=1.d0
                  if(vijrij.lt.0.d0) then
                     wij=vijrij/rij
                     !vijsignal=alfaAV*(va_sroot(i)+va_sroot(j))/2.d0-betaAV*wij
                     vijsignal=alfaAV*(va_sroot(i)+va_sroot(j))/2.d0-betaAV*wij
                     ro3dij=0.5d0*(ro3d(i)+ro3d(j))
                     !vijsignalu=sqrt(abs(p(i)-p(j))/ro3dij)
                     vijsignalu=abs(wij)
                     if(vijsignal.lt.0.d0) vijsignal=0.d0
                     maxvsignal(i)=max(maxvsignal(i),vijsignal)
                     divi=abs(divv(i))
                     divj=abs(divv(j))
                     roti=abs(curlz(i))
                     rotj=abs(curlz(j))
                     if(balsara_lim .and.l.ge.initial_adjust_iterations) then
                        fbalsi=divi/(divi+roti+1.d-5*c(i)/h(i))
                        fbalsj=divj/(divj+rotj+1.d-5*c(j)/h(j))
                        fij=max(0.05d0,0.5*(fbalsi+fbalsj))
                        fbalsi=fij
                        fbalsj=fij
                     else
                        fbalsi=1.d0
                        fbalsj=1.d0
                     endif
                     qiji=-vijsignal*wij/2.d0/pi
                     qijiu=alfau*vijsignalu*(u(i)-u(j))
                     qijj=qiji
                     qijju=qijiu
                     qiji=qiji*fbalsi
                     qijj=qijj*fbalsj
                  else
                     qiji=0.d0
                     qijj=0.d0
                     qijiu=0.d0
                     qijju=0.d0
                  endif        !closing the AV algorithm

                  if(.not.crossedMom) then
                     proi=pro(i)*di/ro2d(i)
                     proj=pro(j)*dj/ro2d(j)
                     proeneri=proener(i)*di/ro2d(i)
                     Sij_i11=Sij(i,1,1)/ro2d(i)
                     Sij_i12=Sij(i,1,2)/ro2d(i)
                     Sij_i13=Sij(i,1,3)/ro2d(i)
                     Sij_i21=Sij(i,2,1)/ro2d(i)
                     Sij_i22=Sij(i,2,2)/ro2d(i)
                     Sij_i23=Sij(i,2,3)/ro2d(i)
                     Sij_i31=Sij(i,3,1)/ro2d(i)
                     Sij_i32=Sij(i,3,2)/ro2d(i)
                     Sij_i33=Sij(i,3,3)/ro2d(i)
                     Sij_j11=Sij(j,1,1)/ro2d(j)
                     Sij_j12=Sij(j,1,2)/ro2d(j)
                     Sij_j13=Sij(j,1,3)/ro2d(j)
                     Sij_j21=Sij(j,2,1)/ro2d(j)
                     Sij_j22=Sij(j,2,2)/ro2d(j)
                     Sij_j23=Sij(j,2,3)/ro2d(j)
                     Sij_j31=Sij(j,3,1)/ro2d(j)
                     Sij_j32=Sij(j,3,2)/ro2d(j)
                     Sij_j33=Sij(j,3,3)/ro2d(j)
                  else
                     proi=pro(i)*di/ro2d(j)
                     proj=pro(j)*dj/ro2d(i)
                     proeneri=proener(i)*di/ro2d(j)
                     Sij_i11=Sij(i,1,1)/ro2d(j)
                     Sij_i12=Sij(i,1,2)/ro2d(j)
                     Sij_i13=Sij(i,1,3)/ro2d(j)
                     Sij_i21=Sij(i,2,1)/ro2d(j)
                     Sij_i22=Sij(i,2,2)/ro2d(j)
                     Sij_i23=Sij(i,2,3)/ro2d(j)
                     Sij_i31=Sij(i,3,1)/ro2d(j)
                     Sij_i32=Sij(i,3,2)/ro2d(j)
                     Sij_i33=Sij(i,3,3)/ro2d(j)
                     Sij_j11=Sij(j,1,1)/ro2d(i)
                     Sij_j12=Sij(j,1,2)/ro2d(i)
                     Sij_j13=Sij(j,1,3)/ro2d(i)
                     Sij_j21=Sij(j,2,1)/ro2d(i)
                     Sij_j22=Sij(j,2,2)/ro2d(i)
                     Sij_j23=Sij(j,2,3)/ro2d(i)
                     Sij_j31=Sij(j,3,1)/ro2d(i)
                     Sij_j32=Sij(j,3,2)/ro2d(i)
                     Sij_j33=Sij(j,3,3)/ro2d(i)
                  endif
                  !Separate diagonal and Off-diagonal parts of Sij
                  momentumxxp=proi*gradxha+proj*gradxhb
                  momentumyyp=proi*gradyha+proj*gradyhb
                  momentumxx=Sij_i11*gradxhab*di+Sij_j11*gradxhab*dj+Sij_i12*gradyhab*di+Sij_j12*gradyhab*dj
                  momentumyy=Sij_i21*gradxhab*di+Sij_j21*gradxhab*dj+Sij_i22*gradyhab*di+Sij_j22*gradyhab*dj
                  momentumzz=Sij_i31*gradxhab*di+Sij_j31*gradxhab*dj+Sij_i32*gradyhab*di+Sij_j32*gradyhab*dj
                  dmy01=vol(i)*mass(j)/mass(i)*qiji
                  dmy02=vol(j)*qijj
                  aviscox=0.5d0*(dmy01*gradxha+dmy02*gradxhb)
                  aviscoy=0.5d0*(dmy01*gradyha+dmy02*gradyhb)

                  dmy11=vol(i)*mass(j)/mass(i)*qijiu
                  dmy12=vol(j)*qijju
                  dpxavu=0.5d0*(dmy11*gradxha+dmy12*gradxhb)
                  dpyavu=0.5d0*(dmy11*gradyha+dmy12*gradyhb)

                  sum1(i)= sum1(i) + mass(j)*momentumxxp+aviscox
                  sum2(i)= sum2(i) + mass(j)*momentumyyp+aviscoy

                  sum1B(i)=sum1B(i)+mass(j)*momentumxx
                  sum2B(i)=sum2B(i)+mass(j)*momentumyy
                  sum3B(i)=sum3B(i)+mass(j)*momentumzz
                  !If necessary substracts the extra terms to convert IAD0 to IAD and mitigate the tensile instability
                  if(.not.iad0) then
                     iadxx=mass(j)*phi(i)*(Sij_i11*gradxhab+Sij_i12*gradyhab)*di
                     iadyy=mass(j)*phi(i)*(Sij_i21*gradxhab+Sij_i22*gradyhab)*di
                     iadzz=mass(j)*phi(i)*(Sij_i31*gradxhab+Sij_i32*gradyhab)*di
                     sum1B(i)=sum1B(i)-iadxx
                     sum2B(i)=sum2B(i)-iadyy
                     sum3B(i)=sum3B(i)-iadzz
                  endif
                  Bij2=(B(i,1)-B(j,1))**2+(B(i,2)-B(j,2))**2+(B(i,3)-B(j,3))**2
                  Fabha=(-djix*gradxha-djiy*gradyha)/rij
                  Fabhb=(-djix*gradxhb-djiy*gradyhb)/rij
                  
                  !Dissipation is always acting on B. Obtain (dB/dt)_diss  (cartesian part called Bindis calculated as in  Weiss&Shen2020):
                  vijsignalBi=sqrt(valfven2(i))
                  vijsignalBj=sqrt(valfven2(j))
                  !*******************************************************
                  dmyr=-vij3*dmy2
                  dmyz=vij3*dmy1
                  dmyphy=vij1*dmy2-vij2*dmy1
                  vijsignaldis=2.d0*sqrt(dmyr**2+dmyz**2+dmyphy**2)/rij  !!A variant of the magnetic signal velocity in dissipation
                  !*******************************************************
!!                 vijsignaldis=(vijsignalBi+vijsignalBj)  !! signal velocity = alfven velocity
                  viscBdis=0.5d0*alfaB*vijsignaldis
                  viscBi=0.5*alfaB*vijsignaldis
!!                 viscBi=0.5d0*alfaB*vijsignalBi
!!                 viscBj=0.5d0*alfaB*vijsignalBj
                  dmy=vol(j)*viscBdis*(Fabha+Fabhb)/2.d0
                  dmyB=dmy
                  dmyBclean=vol(j)*(phidiv(i)+phidiv(j))
                  Bindis(i,1)=Bindis(i,1)+(B(i,1)-B(j,1))*dmyB
                  Bindis(i,2)=Bindis(i,2)+(B(i,2)-B(j,2))*dmyB
                  Bindis(i,3)=Bindis(i,3)+(B(i,3)-B(j,3))*dmyB
                  Bdivclean(i,1)=Bdivclean(i,1)+dmyBclean*gradxhab
                  Bdivclean(i,2)=Bdivclean(i,2)+dmyBclean*gradyhab
                  Bdivclean(i,3)=0.d0

                  dmyB=-dmy/ro3d(i)/mu0*Bij2  !It looks like there is a 1/2 constant missing here but it is placed in the energy Eq.

                  eviscoB(i)=eviscoB(i)+dmyB

                  tqij=mass(j)*(vij1*gradxha+vij2*gradyha)
                  sum5(i)=sum5(i)+tqij*proeneri
                  evisco(i)=evisco(i)+ aviscox*vij1+aviscoy*vij2
                  eviscu(i)=eviscu(i)+dpxavu*uijx+dpyavu*uijy
               enddo   ! Close momentum neighbors loop

               !add non-axial part to Bindis
               smdiv=epsdiv*h(i)
               Bdist=min(di,h(i))
               Bindis(i,1)=Bindis(i,1)+viscBi*Bdist*(gradrBr(i)-B(i,1)/di)/(di+smdiv)
               Bindis(i,2)=Bindis(i,2)+viscBi*Bdist*(gradrBz(i)/(di+smdiv))
               Bindis(i,3)=Bindis(i,3)+viscBi*Bdist*(gradrBphi(i)-B(i,3)/di)/(di+smdiv)
            enddo !Close main momentum loop
            !$omp end do
            !$omp end parallel

            if(profiling) call profile(1,'    mom_ener')
            
            !         HOOP-STRESS TERMS
            if(profiling) call profile(0,'')
            !$omp parallel private(i,di,Btot2,hoop_r,hoop_phi,hoop)
            !$omp do
            do i=1,n1
               di=a(i,1)-desplr
               !with magnetic field in the Hoop-stress
               Btot2=B(i,1)**2+B(i,2)**2+B(i,3)**2
               hoop_r=2.d0*pi*((p(i)+Btot2/2.d0/mu0)-B(i,3)**2/mu0)/ro2d(i)
               sum1(i)=hoop_r-2.d0*pi*sum1(i)+2.d0*pi*sum1B(i)
               sum2(i)=-2.d0*pi*sum2(i)+2.d0*pi*sum2B(i)
               hoop_phi=2.d0*pi*(B(i,1)*B(i,3)/mu0)/ro2d(i)
               sum3(i)=2.d0*pi*sum3B(i)+hoop_phi
               evisco(i)=max(0.d0,evisco(i))
                hoop=-2.d0*2.d0*pi*p(i)/ro2d(i)*v(i,1)
                sum5(i)=hoop+2.d0*2.d0*pi*sum5(i)+ 2.d0*pi*evisco(i)+ 2.d0*eviscu(i)
             enddo
             !$omp end do
             !$omp end parallel
             if(profiling) call profile(1,' Hoop_stress')
             

             !Average of gradient_r  around the axis
             do i=1,n
                if(axis(i)) then
                   sum1(i)=-sum1(ngh(i))
                else
                   sum1(i)=sum1(ngh(i))
                endif
             enddo


             !           INDUCTION EQUATION
             if(profiling) call profile(0,'')
             !$omp parallel private(i,di,ratio,smdiv,Bcx,Bcz,Bcphi,lcentered,smdi,&
             !$omp    a11,a12,a13,a21,a22,a23,a31,a32,a33)
             !$omp do
             do i=1,n1
                di=abs(a(i,1)-desplr)
                ratio=di/h(i)
                smdiv=epsdiv*h(i)
                a11=-(divvzz(i)+v(i,1)/(di+smdiv))
                a12=divvrz(i)
                a13=-v(i,3)/(di+smdiv)
                a21=divvzr(i)
                a22=-(divvrr(i)+v(i,1)/(di+smdiv))
                a23=0.d0
                a31=divvphir(i)
                a32=divvphiz(i)
                a33=-(divvrr(i)+divvzz(i))
                Bcx=B(i,1)
                Bcz=B(i,2)
                Bcphi=B(i,3)
                B(i,1)=B0(i,1)+(a11*Bcx+a12*Bcz+a13*Bcphi+Bindis(i,1))*dt-Bdivclean(i,1)*dt
                B(i,2)=B0(i,2)+(a21*Bcx+a22*Bcz+a23*Bcphi+Bindis(i,2))*dt-Bdivclean(i,2)*dt
                B(i,3)=B0(i,3)+(a31*Bcx+a32*Bcz+a33*Bcphi+Bindis(i,3))*dt-Bdivclean(i,3)*dt
                lcentered=lz(i)
                smdi=di+epsdiv*h(i)
                v(i,1)=v0(i,1)+sum1(i)*dt+lcentered**2/smdi**3*dt
                v(i,2)=v0(i,2)+sum2(i)*dt
                lz(i)=lz0(i)+smdi*sum3(i)*dt
                v(i,3)=lz(i)/smdi
                u(i)=max(1.d-6, u0(i)+0.5d0*sum5(i)*dt+0.5*eviscoB(i)*dt)
             enddo
             !$omp end do
             !$omp end parallel
             if(profiling) call profile(1,'   Induction')


             !make the average of radial velocity near the axis
             if(averVelr) then
                if(profiling) call profile(0,'')
                !$omp parallel private(i,di,ratio,jVec,j,d1,d2,d02,d05,v1)
                !$omp do
                do i=1,n1
                   di=a(i,1)-desplr
                   ratio=abs(di)/h(i)
                   nvvv(i)=1
                   avevel(i)=0.d0
                   bfieldr(i)=0.d0
                   if(ratio.le.2.d0) then
                      do jVec=1,nv(i)
                         j=vecMatrix(i,jVec)
                         d1=abs(a(i,1)-a(j,1))
                         d2=abs(a(i,2)-a(j,2))
                         d02=d1**2+d2**2
                         d05=sqrt(d02)
                         v1=d05/h(i)
                         if(v1.ne.0.d0.and.v1.lt.2.d0) then
                            nvvv(i)=nvvv(i)+1
                            avevel(i)=avevel(i)+v(j,1)
                            bfieldr(i)=bfieldr(i)+B(j,1)
                         endif
                      enddo
                      avevel(i)=avevel(i)+v(i,1)  !Self-contribution
                      bfieldr(i)=bfieldr(i)+B(i,1)
                   endif
                enddo
                !$omp end do

                !$omp do
                do i=1,n1
                   di=a(i,1)-desplr
                   if(nvvv(i).gt.1) then
                      v(i,1)=avevel(i)/nvvv(i)
                      B(i,1)=bfieldr(i)/nvvv(i)
                      if(di/h(i).le.0.5) B(i,1)=0.d0
                   endif
                enddo
                !$omp end do
                !$omp end parallel
                if(profiling) call profile(1,' Aver_Bfield')

             endif
             do i=1,n1
                a(i,1) = a0(i,1) + ((v(i,1)+v0(i,1))/2.d0)*dt
                a(i,2) = a0(i,2) + ((v(i,2)+v0(i,2))/2.d0)*dt
                if(a(i,1)-desplr .lt.0.d0) then
                   print*, 'crossing the axis !!'
                   print*, l,kl,i
                   print*, a(i,1)-desplr,a(i,2)-desplz
                   print*, a0(i,1)-desplr,a0(i,2)-desplz,ro3d(i)
                   print*, u(i),p(i),c(i)
                   print*, u0(i),h(i)
                   print*, v(i,1),v(i,2)
                   print*, v0(i,1),v0(i,2)
                   stop
                endif
             enddo


             !Apply PBC
             do i=1,n1
                if(a(i,2)-desplz .gt. rtzu) a(i,2)=a(i,2)-(rtzu+rtzd)
                if(a(i,2)-desplz .lt. -rtzd) a(i,2)=a(i,2)+(rtzu+rtzd)
             enddo
             
             do i=1,n2
                ngh(i)=i
                axis(i)=.false.
             enddo

             do jg=1,nc
                call ghosts(desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5,jg,n,&
                     &    a,v,B,n1,n2,axis,ngh)
                do i=1,n
                   h(i)=h(ngh(i))
                   mass(i)=mass(ngh(i))
                   u(i)=u(ngh(i))
                   ro2d(i)=ro2d(ngh(i))
                enddo
             enddo
          end do

          !     **************************************************************
          !      END OF MODEL CALCULATION
          !    *************************************************************
          do i=1,n1
             ro3d(i)=ro2d(i)/(2.d0*pi*(a(i,1)-desplr))
          enddo

          ! Center of mass ***
          
          ycm = 0.0d0
          dmass=0.d0
          do i = 1,n1
             dmass=dmass+mass(i)
             ycm = ycm + a(i,2)*mass(i)
          end do
          ycm_check=ycm/dmass
          ycm=desplz
          xcm=desplr
          
          !Find max density
          rhomax=maxval(ro3d, DIM=1)
          irhomax=maxloc(ro3d, DIM=1)
          
          !Energy conservation
          
          ekin=0.d0
          eint=0.d0
          emag=0.d0
          egrav=0.d0
          romax=-1.d0
          errdivBav=0.d0
          errdivBmax=0.d0
          errdivBmin=100000.d0
          angmomr=0.d0
          angmomz=0.d0
          angmomphi=0.d0
          Bradav=0.d0
          Bzav=0.d0
          Bphiav=0.d0
          masstot=0.d0
          nav=0
          
          do i=1,n1
             angmomr=angmomr-mass(i)*lz0(i)/(a(i,1)-xcm)*(a(i,2)-ycm)
             angmomz=angmomz+mass(i)*lz0(i)
             angmomphi=angmomphi+mass(i)*((a(i,2)-ycm)*v(i,1)-(a(i,1)-xcm)*v(i,2))
             vmod=(v(i,1)**2+v(i,2)**2+v(i,3)**2)
             ekin=ekin+0.5d0*mass(i)*vmod
             eint=eint+mass(i)*u(i)
             Bmo12=dsqrt(B(i,1)**2+B(i,2)**2+B(i,3)**2)
             emag=emag+0.5d0*Bmo12**2/mu0*mass(i)/ro3d(i)
             if(Bmo12 .ne. 0.d0) then
                errdivB=h(i)*abs(divB(i))/Bmo12
                Bradav=Bradav+abs(B(i,1))*mass(i)
                Bzav=Bzav+abs(B(i,2))*mass(i)
                Bphiav=Bphiav+abs(B(i,3))*mass(i)
             else
                errdivB=0.d0
             endif
             if(errdivB.ge.errdivBmax) errdivBmax=errdivB
             if(errdivB.le.errdivBmin) errdivBmin=errdivB
             errdivBav=errdivBav+errdivB*mass(i)
             nav=nav+1
             masstot=masstot+mass(i)
             if(ro3d(i).gt.romax) romax=ro3d(i)
          end do
          if(l.le.initial_adjust_iterations) then
             etot0=ekin+eint+emag
             angmomtot_0=dsqrt(angmomr**2+angmomz**2+angmomphi**2)
             cm0=abs(ycm_check-desplz)
          endif
          angmomtot=dsqrt(angmomr**2+angmomz**2+angmomphi**2)
          cm=abs(ycm_check-desplz)
          open(4,file='ZPinchnb60n5Ener_02',access='append')
          elost=(ekin+eint+emag)-etot0
          write(4,137)l,tt,ekin,eint,egrav,emag,elost/etot0,&
               &   errdivBmin,errdivBmax,errdivBav,(cm-cm0)/4.6d17,&
               &  abs(angmomtot-angmomtot_0)/angmomtot_0,angmomtot_0,&
               &  angmomtot,romax,Bradav,Bzav,Bphiav
          close(4)
137       format(2x,i8,17(2x,e16.9))
          

          !Timestep control
          ! In acceleration
          time_accel=1.d50
          do i = 1,n1
             accel=sqrt(sum1(i)**2+sum2(i)**2)
             if(accel.eq.0.d0) cycle
             dmy=C_accel*sqrt(h(i)/accel)
             time_accel=min(time_accel,dmy)
          end do

          ! In temperature
          time_temp=1.d50
          if(timeTemp) then
             do i = 1,n1
                du=abs(u(i)-u0(i))
                if(du.eq.0.d0) cycle
                dmy=(u(i)*dtnma*tol_temp)/du
                time_temp=min(time_temp,dmy)
             end do
          endif

          ! Courant condition
          time_Cour=1.d50
          do i = 1,n1
             if(maxvsignal(i).gt.0.d0) then
                dmy=min(Cour*h(i)/va_sroot(i),Cour*h(i)/maxvsignal(i))
             else
                dmy=Cour*h(i)/va_sroot(i)
             endif
             time_Cour=min(time_Cour,dmy)
          end do
          
          ! In the axis-crossing time in convergent particles.
          time_axis=1.d50
          if(timeaxis) then
             do i = 1,n1
                if(v(i,1).ge.0.d0) cycle
                dmy=tol_taxis*(a(i,1)-desplr)/abs(v(i,1))
                time_axis=min(time_axis,dmy)
             end do
          endif

          ! Final time-step
          dtnme=dtnma
          dtnma = min(time_accel,time_temp,time_Cour,time_axis)
          if(dtnma .ge. 1.05d0*dtnme) dtnma=1.05d0*dtnme
          if(dtnma .le. 0.2d0*dtnme) dtnma=0.2d0*dtnme
          if(dtnma.le.2.d-08) dtnma=2.d-08
          dtn=0.5d0*(dtnme+dtnma)
          tt=tt+dtnma
          ! Init for next time-step
          a0 = a
          v0 = v
          u0 = u
          B0 = B
          phidiv0 = phidiv
          cha0 = cha
          lz0 = lz
          
          ! *************************************************************
          !  Output
          ! *************************************************************
          
          if(mod(l,iodelay).eq.0 .or. l.eq.lstart) then
             prefix='Zpi_2.'
             !value1=idint(tt*1.e-8)
             value1=l
             write(suffix,'(i6.6)') value1
             nomfs1=prefix//suffix
             
             !call nomfils(nomfs1,nomfs2,value1,prefix)
             open(10,file=nomfs1)
             do i =1,n1
                ax=a(i,1)-xcm
                ay=a(i,2)-ycm
                rrr=sqrt(ax**2+ay**2)
                ux=ax/rrr
                uy=ay/rrr
                vr=v(i,1)*ux+v(i,2)*uy
                checkdeltar=sqrt(checkdeltax(i)**2+checkdeltay(i)**2)
                !B(i,1)=max(B(i,1),minvalue)
                !B(i,2)=max(B(i,2),minvalue)
                !B(i,3)=max(B(i,3),minvalue)
                if(abs(B(i,1)) .le.minvalue) B(i,1)=minvalue
                if(abs(B(i,2)) .le.minvalue) B(i,2)=minvalue
                if(abs(B(i,3)) .le.minvalue) B(i,3)=minvalue
                divB(i)=max(divB(i),minvalue)
                write(10,'(39(1x,e23.16),1x,i4)') ax,ay,rrr,mass(i),h(i),&
                     &       c(i),ro3d(i),ro2d(i),v(i,1),v(i,2),&
                     &       v(i,3),u(i),p(i),sum1(i),sum2(i),&
                     &       checkInorm(i),divv(i),divvrr(i),curlz(i),sum5(i),&
                     &       evisco(i),eviscu(i),checkdeltar,B(i,1),B(i,2),&
                     &       B(i,3),phi(i),sum1B(i),sum2B(i),Bindis(i,1),&
                     &       Bindis(i,2),Bindis(i,3),divB(i),phidiv(i),divvzz(i),&
                     &       gradrBphi(i),vr,lz(i),&
                     &       sum3(i),nv(i)
             enddo
             close(10)
          endif

       enddo !main loop
       stop
     END PROGRAM AxisSPHYNX_MHD


     subroutine func(pi2x,w)
       IMPLICIT NONE
       double precision, intent(in)::pi2x
       double precision, intent(out)::w
       if(pi2x.eq.0.d0) then
          w=1.d0
       else
          w=sin(pi2x)/pi2x
       endif
       RETURN
     END subroutine func


     subroutine dfunc(pi2,pi2x,dw)
       IMPLICIT NONE
       double precision, intent(in):: pi2,pi2x
       double precision, intent(out)::dw
       if(pi2x.eq.0.d0)then
          dw=0.d0
       else
          dw=pi2*(1.d0/tan(pi2x)-1.d0/(pi2x))
       endif
       RETURN
     END subroutine dfunc


!                                SUBROUTINE NOMFILS
!    --------------------------------------------------------------
!
     subroutine nomfils(nomfs,nomfs2,i,prefix)
       implicit none
       integer,INTENT(in)::i
       integer decmillar,millar,cent,dec,uni,rest
       character,INTENT(out)::nomfs*14,nomfs2*14
       character,INTENT(in)::prefix*6
       character sufix*5
       decmillar=i/10000+48
       rest=mod(i,10000)
       millar=rest/1000+48
       rest=mod(i,1000)
       cent=rest/100+48
       rest=mod(i,100)
       dec=rest/10+48
       uni=mod(i,10)+48
       sufix=char(decmillar)//char(millar)//char(cent)//char(dec)//&
            &    char(uni)
       nomfs=prefix//sufix
       nomfs2='sed.'//sufix
       return
     end subroutine nomfils

     subroutine ghosts(desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5,jg,ng,&
          &    a1,v1,B1,n1,nt,axis,ngh)
       implicit none
       
       integer i
       double precision dx,dy
       
       integer, intent(in)  :: jg,n1,nt
       integer, intent(out) :: ng
       integer, dimension(nt), intent(out) :: ngh
       
       double precision, intent(in) :: desplr,desplz,rtr,rtzu,rtzd,rtr5,rtz5
       
       logical, dimension(nt), intent(out) :: axis
       
       double precision, dimension(nt,2), intent(inout) :: a1
       double precision, dimension(nt,3), intent(inout) :: v1,B1

       !         axis-reflective ghosts
       if(jg.eq.1) then
          ng=n1
          do i=1,n1
             dx=a1(i,1)-desplr
             if(dx.le.rtr5) then
                ng=ng+1
                a1(ng,1)=-a1(i,1)+2.d0*desplr
                a1(ng,2)=a1(i,2)
                v1(ng,1)=-v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=-v1(i,3)
                B1(ng,1)=-B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)=-B1(i,3)
                axis(ng)=.true.
                ngh(ng)=i
             endif
          end do
       endif
       !           box 2  (periodic ghosts up)
       if(jg.eq.2) then
          do i=1,n1
             dy=a1(i,2)-desplz
             if(dy.ge.rtz5*4.d0) then
                ng=ng+1
                a1(ng,1)=a1(i,1)
                a1(ng,2)=a1(i,2)-(rtzu+rtzd)
                v1(ng,1)=v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=v1(i,3)
                B1(ng,1)=B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)=B1(i,3)
                axis(ng)=.false.
                ngh(ng)=i
             endif
          end do
       endif
       !           box 3 (periodic ghosts down)
       if(jg.eq.3) then
          do i=1,n1
             dy=a1(i,2)-desplz
             if(dy.le.-rtz5*4.d0) then
                ng=ng+1
                a1(ng,1)=a1(i,1)
                a1(ng,2)=(rtzu+rtzd)+a1(i,2)
                v1(ng,1)=v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=v1(i,3)
                B1(ng,1)=B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)=B1(i,3)
                axis(ng)=.false.
                ngh(ng)=i
             endif
          end do
       endif
       !           Box 4 (Reflective ghosts right.
       !                  Have to be reflective to preserve the mass of the particles)
       if(jg.eq.4) then
          do i=1,n1
             dx=a1(i,1)-desplr
             if(dx.ge.rtr5*4.d0) then
                ng=ng+1
                a1(ng,1)=2.d0*rtr-a1(i,1)+2.d0*desplr
                a1(ng,2)=a1(i,2)
                v1(ng,1)=-v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=v1(i,3)
                B1(ng,1)=B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)=B1(i,3)
                axis(ng)=.false.
                ngh(ng)=i
             endif
          end do
       endif
       !           Box 5 (periodic ghosts upper-left)
       if(jg.eq.5) then
          do i=1,n1
             dx=a1(i,1)-desplr
             dy=a1(i,2)-desplz
             if(dx.le.rtr5 .and. dy.ge.rtz5*4.d0) then
                ng=ng+1
                a1(ng,1)=-a1(i,1)+2.d0*desplr
                a1(ng,2)=a1(i,2)-(rtzu+rtzd)
                v1(ng,1)=-v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=-v1(i,3)
                B1(ng,1)=-B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)=-B1(i,3)
                axis(ng)=.true.
                ngh(ng)=i
             endif
          end do
       endif
       !           Box 6 (reflective ghosts upper-right. First periodic up then make reflective left.
       !                  (have to be reflective to preserve the mass of the particle))
       if(jg.eq.6) then
          do i=1,n1
             dx=a1(i,1)-desplr
             dy=a1(i,2)-desplz
             if(dx.ge.rtr5*4.d0 .and. dy.ge.rtz5*4.d0)then
                ng=ng+1
                a1(ng,1)=2.d0*rtr-a1(i,1)+2.d0*desplr
                a1(ng,2)=a1(i,2)-(rtzu+rtzd)
                v1(ng,1)=-v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=v1(i,3)
                B1(ng,1)=B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)=B1(i,3)
                axis(ng)=.false.
                ngh(ng)=i
             endif
          end do
       endif
       !        Box 7 (Periodic ghosts lower-left)
       if(jg.eq.7) then
          do i=1,n1
             dx=a1(i,1)-desplr
             dy=a1(i,2)-desplz
             if(dx.le.rtr5 .and. dy.le.-rtz5*4.d0) then
                ng=ng+1
                a1(ng,1)=-a1(i,1)+2.d0*desplr
                a1(ng,2)=(rtzu+rtzd)+a1(i,2)
                v1(ng,1)=-v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=-v1(i,3)
                B1(ng,1)=-B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)=-B1(i,3)
                axis(ng)=.true.
                ngh(ng)=i
             endif
          end do
       endif
       !           Box 8 (Reflective ghosts lower-right. First periodic down. Then reflective left.)
       if(jg.eq.8) then
          do i=1,n1
             dx=a1(i,1)-desplr
             dy=a1(i,2)-desplz
             if(dx.ge.rtr5*4.d0 .and. dy.le.-rtz5*4.d0)then
                ng=ng+1
                a1(ng,1)=2.d0*rtr-a1(i,1)+2.d0*desplr
                a1(ng,2)=(rtzu+rtzd)+a1(i,2)
                v1(ng,1)=-v1(i,1)
                v1(ng,2)=v1(i,2)
                v1(ng,3)=v1(i,3)
                B1(ng,1)=B1(i,1)
                B1(ng,2)=B1(i,2)
                B1(ng,3)= B1(i,3)
                axis(ng)=.false.
                ngh(ng)=i
             endif
          end do
       endif
       return
     end subroutine ghosts

     subroutine gravity(n_grav_rows,tp,gi1,gi2,gi3)

       implicit none
       integer i
       integer, intent(in) :: n_grav_rows
       double precision, dimension(n_grav_rows), intent(out) :: tp,gi1,gi2,gi3
       
       open(33,file='I1I2.txt')
       do i=1,n_grav_rows
          read(33,*) tp(i),gi1(i),gi2(i),gi3(i)
       enddo
       close(33)
       return
     end subroutine gravity

     subroutine profile(mode,label)
       integer, intent(in):: mode
       character :: label*12
       double precision start
       double precision omp_get_wtime
       double precision end
       
       save start
       
       if(mode.eq.0) then
          start=omp_get_wtime()
       else if(mode.eq.1) then
       end if=omp_get_wtime()
       print '(a12,a2,es12.5,a2)',label(1:),': ',end-start,' s'
    endif
  end subroutine profile
