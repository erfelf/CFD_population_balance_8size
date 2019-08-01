	Program FV_PBE_lax_fred
	Use List_of_Variables	
	Implicit none
	!********************************************************!
	!	****	Population Balance Modeling	****	 !
	!********************************************************!
	write (*,*) 'Please choose the test:'
	write (*,*) '(1) Sets of particle in the center'
	write (*,*) '(2) Full domain particle'
	write (*,*) '(3) particle jet'
c	read (*,*) Flag
	Flag=2
	IF (Flag .EQ. 3) then	
		write (*,*) 'Please tell us the input velocity "Uin"'
c		read (*,*) Upin
		Upin=0.5
        	Num_mesh_in=160
        	max_Den_p_in=4000
	End IF
	write (*,*) 'Please choose the behaviour of the background fluid'	
	write (*,*) '(1) Stationary fluid'
	write (*,*) '(2) Constant velocity in y direction'
	write (*,*) '(3) Constant velocity in x direction'
	write (*,*) '(4) Couette Flow NW'	
	write (*,*) '(5) Channel Flow NW'
c	read (*,*) Flag4
	Flag4=1
	IF (Flag4 .GT. 1) then	
		write (*,*) 'Please tell us the input velocity "Uin"'
c		read (*,*) Uin
		Uin=-0.1
	End IF
	write (*,*) 'Please choose particle forces'
	write (*,*) '(1) Just Gravity'
	write (*,*) '(2) Gravity and Buoyancy'
	write (*,*) '(3) Gravity and Drag '
	write (*,*) '(4) Gravity, Buoyancy and Drag'
	write (*,*) '(5) Gravity, Buoyancy, Drag and lift NW'
c   	read (*,*) Flag1
	Flag1=2

	write (*,*) 'Please choose if you want to activate aggregation'
	write (*,*) '(1)',Num_family,'sets of particles, no aggregation'
	write (*,*) '(2)',Num_family,'sets of particles and aggregation'
	write (*,*) '(3)',Num_family,'sets of particles+ aggregation and
     1  dissaggergation'	
	write (*,*) 'Note:for more or less particle open List_of_Varible' ! update the List_of_Varible if more particle sets or less is needed
c  	read (*,*) Flag2
	Flag2=1
	write (*,*) 'Please choose the option for background fluid'	
	write (*,*) '(1) Deactivating fluid flow in the background'
	write (*,*) '(2) Activating fluid flow but no Interaction force'
	write (*,*) '(3) Activating fluid flow with Interaction forces'	
c	read (*,*) Flag3
	Flag3=1

	write (*,*) 'Please choose the hight of the domian'	
c	read (*,*) H
	H=0.0002
	write (*,*) 'Please choose the length of the domian'	
c	read (*,*) L
	L=0.0002

	write (*,*) 'Please specify the time of the simulation'
	write (*,*) 'Note: choose small number for no drag (0.1-0.2 s)'
c	read (*,*) t_final
	t_final=1.0
        deltat_start=0.01
C       !Creating files to write the data
	Open(unit=1, file='Particles.dat')
C       Should be updated if the number of classes changed
	Write (1,*) ' VARIABLES = "x","y","alpha","rho1","rho2","rho3"'   !Update

	Open(unit=3, file='Fluid.dat')
	Write (3,*) ' VARIABLES = "x","y","alpha","U","V","P" ' 

	Call read_mesh
  	Call Set_Initial_Condition
	Call Set_test_Condition
	if (flag2 .EQ. 3) then
	    CFL=0.001
	End If 

	current_time = 0.0
        
	Save_time=0.0 ! time when the results will be written to the files above

C	!-----------Time Loop 
	Do while (t_final .GT. current_time)
	  deltat=deltat_start ! will be updated according to particle deltat 
	  If (current_time+deltat .GT. t_final) then   ! last step
	    deltat = t_final - current_time
	  End If

	  ! Initialization for new time step 
	  Source(:,:)=0.0
	  dUdt(:,:)=0.0
	  dUdx(:,:)=0.0
	  dUdy(:,:)=0.0
	  dAggdt(:,:)=0.0
	  pressure(:)=0.0
	
C       !----Defining the artificial Pressure
	! p=alpha/(10^8*(max_pack-alpha)^9)  
	  Do i=1,nxp
	    Do j=1,nyp
	      index1 = i + nxp*(j-1)
  	      pressure(index1)=alpha(index1)/
     1          ((1.0D8*(max_Pack-alpha(index1)))**9)
	      di_pressure_di_alpha(index1)=alpha(index1)/
     1          (2.0D8*(alpha(index1)-max_Pack)**10)-1/
     2	        (1.0D8*(alpha(index1)-max_Pack)**9)
	    End Do
	  End Do
C-------------------------Solving Particles-------------------------------------
C-------------------------------------------------------------------------------
C         ! minmod
c	Do k=1,Num_family
c	   Do i=1,nxp
c	        j=1
c		index1= i + nxp*(j-1)
c		Ud(1)=U(index1,3*k-2)
c		Ud(2)=U(index1,3*k-1)
c		Ud(3)=-1*U(index1,3*k)
c		fdiff=(U(index1+nxp,3*k-2:3*k)-U(index1,3*k-2:3*k))/Sns(j+1)
c		bdiff=(U(index1,3*k-2:3*k)-Ud)/(0.5*Sns(j+1))
c		Call minmod
c		dUdy(index1,3*k-2:3*k)=d_ud_loc
c		Do j=2,nyp-1
c		   index1= i + nxp*(j-1)
cc		   fdiff=(U(index1+nxp,3*k-2:3*k)-U(index1,3*k-2:3*k))/Sns(j+1)
c		   bdiff=(U(index1,3*k-2:3*k)-U(index1-nxp,3*k-2:3*k))/Sns(j)
c		   Call minmod
c		   dUdy(index1,3*k-2:3*k)=d_ud_loc
c		End Do
c		j=nyp
c		index1= i + nxp*(j-1)
c		Uu(1)=U(index1,3*k-2)
c		Uu(2)=U(index1,3*k-1)
c		Uu(3)=-1*U(index1,3*k)
c		fdiff=(Uu-U(index1,3*k-2:3*k))/(0.5*Sns(j))
c		bdiff=(U(index1,3*k-2:3*k)-U(index1-nxp,3*k-2:3*k))/Sns(j)
c		Call minmod
c		dUdy(index1,3*k-2:3*k)=d_ud_loc
c	    End Do	
c	    Do j=1,nyp
c		i=1
c		index1= i + nxp*(j-1)
c		Ul(1)=U(index1,3*k-2)
c		Ul(2)=-1*U(index1,3*k-1)
c		Ul(3)=U(index1,3*k)
c		fdiff=(U(index1+1,3*k-2:3*k)-U(index1,3*k-2:3*k))/Sew(i+1)
c		bdiff=(U(index1,3*k-2:3*k)-Ul)/(0.5*Sew(i+1))
c		Call minmod
c		dUdx(index1,:)=d_ud_loc
c		Do i=2,nxp-1
c		    index1= i + nxp*(j-1)
c		    fdiff=(U(index1+1,3*k-2:3*k)-U(index1,3*k-2:3*k))/Sew(i+1)
c		    bdiff=(U(index1,3*k-2:3*k)-U(index1-1,3*k-2:3*k))/Sew(i)
c		    Call minmod
c		    dUdx(index1,3*k-2:3*k)=d_ud_loc
c		End Do
c		i=nxp
c		index1= i + nxp*(j-1)
c		Ur(1)=U(index1,3*k-2)
c		Ur(2)=-1*U(index1,3*k-1)
c		Ur(3)=U(index1,3*k)
c		fdiff=(Ur-U(index1,3*k-2:3*k))/(0.5*Sew(i))
c		bdiff=(U(index1,3*k-2:3*k)-U(index1-1,3*k-2:3*k))/Sew(i)
c		Call minmod
c		dUdx(index1,3*k-2:3*k)=d_ud_loc
c	    End Do
c	End Do



C         ! LAX_Friedrichs method 
C         ! East Boundary condition for Particles 
	  i=1		
	  Do j=1,nyp
	    index1 = i + nxp*(j-1)
	    Do k=1,Num_family
		if (Flag .EQ.3 .AND. j>nyp/2 .AND. J<nyp/2+num_mesh_in) then
			sigma2=real(num_mesh_in)/3.0
		       j_start=real(j)-real(nyp/2.0)-real(Num_mesh_in/2.0)
		      Ul(1)=max_Den_p_in*exp(-real(j_start)**2.0/(2.0*sigma2)
     1			)/sqrt(2.0*pi*sigma2)
		       Ul(2)=Upin*Ul(1)
		       Ul(3)=0		
		else		
		       Ul(1)=U(index1,3*k-2)
		       Ul(2)=-1*U(index1,3*k-1)
		       Ul(3)=U(index1,3*k)
		End If
		Ur=U(index1,3*k-2:3*k)!-dUdx(index1,3*k-2:3*k)*(Sew(i+1)/2.0)
		pl=pressure(index1)
		pr=pressure(index1)
		d_pl=di_pressure_di_alpha(index1)
                d_pr=di_pressure_di_alpha(index1)
		Call Local_Lax_Friedrichs_x
		dUdt(index1,3*k-2:3*k)=dUdt(index1,3*k-2:3*k)+Fx*Sns(j+1)
 		deltat = min(deltat, CFL*Sns(j+1)/fastest_wave)
	    End Do
	  End Do
C	! Sweeping Particle Domain left to right
	  Do i=2,nxp 
	    Do j=1,nyp
		index1 = i + nxp*(j-1)
	        Do k=1,Num_family
		   Ul=U(index1-1,3*k-2:3*k)!+dUdx(index1,3*k-2:3*k)*(Sew(i)/2.0)
		   Ur=U(index1,3*k-2:3*k)!-dUdx(index1+1,3*k-2:3*k)*(Sew(i+1)/2.0)
		   pl=pressure(index1-1)
		   pr=pressure(index1)
		   d_pl=di_pressure_di_alpha(index1-1)
                   d_pr=di_pressure_di_alpha(index1)
		   Call Local_Lax_Friedrichs_x
		dUdt(index1,3*k-2:3*k)=dUdt(index1,3*k-2:3*k)+Fx*Sns(j+1)
	    dUdt(index1-1,3*k-2:3*k)=dUdt(index1-1,3*k-2:3*k)-Fx*Sns(j+1) 
 		   deltat = min(deltat, CFL*Sns(j+1)/fastest_wave)
		End Do
	    End Do
	  End Do
C	! West boundary Condition for Particles 
	  i=nxp  
	  Do j=1,nyp
	     index1 = i + nxp*(j-1)
	     Do k=1,Num_family
	     if (Flag .EQ. 3 ) then     
		Ur(1)=0
		Ur(2)=0
		Ur(3)=0
	     Else 		
		Ur(1)=U(index1,3*k-2)
		Ur(2)=-1*U(index1,3*k-1)
		Ur(3)=U(index1,3*k)
	     End If
		Ul=U(index1,3*k-2:3*k)!+dUdx(index1,3*k-2:3*k)*(Sew(i)/2.0)
		pl=pressure(index1)
		pr=pressure(index1)
		d_pl=di_pressure_di_alpha(index1)
                d_pr=di_pressure_di_alpha(index1)	
		Call Local_Lax_Friedrichs_x
		dUdt(index1,3*k-2:3*k)=dUdt(index1,3*k-2:3*k)-Fx*Sns(j+1)
 		deltat = min(deltat, CFL*Sns(j+1)/fastest_wave)
	     End Do
	  End Do


C	 ! South Boundary condition for Particles
	  j=1 
	  Do i=1,nxp
	     index1 = i + nxp*(j-1)
	     Do k=1,Num_family
	     If (Flag .EQ. 3 ) then     
		Ud(1)=0
		Ud(2)=0
		Ud(3)=0
	     Else 	
		Ud(1)=U(index1,3*k-2)
		Ud(2)=U(index1,3*k-1)
		Ud(3)=-1*U(index1,3*k)
	     End If
		Uu=U(index1,3*k-2:3*k)!-dUdy(index1,3*k-2:3*k)*(Sns(i+1)/2.0)
		pd=pressure(index1)
		pu=pressure(index1)
		d_pd=di_pressure_di_alpha(index1)
                d_pu=di_pressure_di_alpha(index1)
		Call Local_Lax_Friedrichs_y
		dUdt(index1,3*k-2:3*k)=dUdt(index1,3*k-2:3*k)+Fy*Sew(i+1)
 		deltat = min(deltat, CFL*Sew(i+1)/fastest_wave)
	     End Do
	  End Do
C	!Sweeping Domain South to North
	  Do j=2,nyp 
	     Do i=1,nxp
		index1 = i + nxp*(j-1)
	        Do k=1,Num_family
		   Ud=U(index1-nxp,3*k-2:3*k)!+dUdy(index1,3*k-2:3*k)*(Sns(i)/2.0)
		 Uu=U(index1,3*k-2:3*k)!-dUdy(index1+nxp,3*k-2:3*k)*(Sns(i+1)/2.0)
		   pd=pressure(index1-nxp)
		   pu=pressure(index1)
  		   d_pd=di_pressure_di_alpha(index1-nxp)
                   d_pu=di_pressure_di_alpha(index1)
		   Call Local_Lax_Friedrichs_y
		   dUdt(index1,3*k-2:3*k)=dUdt(index1,3*k-2:3*k)
     1	           +Fy*Sew(i+1)
 	           dUdt(index1-nxp,3*k-2:3*k)=dUdt(index1-nxp,3*k-2:3*k)
     1             -Fy*Sew(i+1) 
		   deltat = min(deltat, CFL*Sew(i+1)/fastest_wave)
		End Do
	      End Do	
	   End Do
C  	! North boundary condition for particles
	   j=nyp
	   Do i=1,nxp
	      index1 = i + nxp*(j-1)
	      Do k=1,Num_family
	     If (Flag .EQ. 3 ) then     
		Uu(1)=0
		Uu(2)=0
		Uu(3)=0
	     Else 
		 Uu(1)=U(index1,3*k-2)
		 Uu(2)=U(index1,3*k-1)
		 Uu(3)=-1*U(index1,3*k)
	     End If
		 Ud=U(index1,3*k-2:3*k)!+dUdy(index1,3*k-2:3*k)*(Sns(i)/2.0)
		 pd=pressure(index1)
		 pu=pressure(index1)
		 d_pd=di_pressure_di_alpha(index1)
                 d_pu=di_pressure_di_alpha(index1)	
		 Call Local_Lax_Friedrichs_y
		dUdt(index1,3*k-2:3*k)=dUdt(index1,3*k-2:3*k)-Fy*Sew(i+1)
 		 deltat = min(deltat, CFL*Sew(i+1)/fastest_wave)
	      End Do
	   End Do

	   Call Get_Source

	   if (Flag2 .GT. 1) then
	       Call Aggregate
           End if

	if (Flag4 .EQ. 2 .OR. Flag4 .EQ. 3) then
	   Do i=1,nxp
	       Do j=1,nyp
		 index1 = i + nxp*(j-1)
	         Ugg(i,j)=(U_fluid(i+1,j+1)+U_fluid(i+2,j+1))/2.0
		 Vgg(i,j)=(V_fluid(i+1,j+1)+V_fluid(i+1,j+2))/2.0

		 Do k=1,Num_family
	   	   U(index1,3*k-2)=U(index1,3*k-2)+deltat*
     1           dUdt(index1,3*k-2)/(Sew(i+1)*Sns(j+1))+deltat*
     2           dAggdt(index1,3*k-2)+deltat*Source(index1,3*k-2)

	   U(index1,3*k-1)=(U(index1,3*k-1)+deltat*dUdt(index1,3*k-1)/
     1     (Sew(i+1)*Sns(j+1))+deltat*dAggdt(index1,3*k-1)+deltat*
     2     Source(index1,3*k-1)-Drag_t(k)*deltat*U(index1,3*k-2)
     3     *Ugg(i,j))/(1.0+Drag_t(k)*deltat)

	   U(index1,3*k)=(U(index1,3*k)+deltat*dUdt(index1,3*k)/
     1     (Sew(i+1)*Sns(j+1))+deltat*dAggdt(index1,3*k)+
     2     deltat*Source(index1,3*k)-Drag_t(k)*deltat*U(index1,3*k-2)
     3     *Vgg(i,j))/(1.0+Drag_t(k)*deltat)
                  End Do
	      End Do
	   End Do
	Else
	  Do j=1,nyp 
	     Do i=1,nxp
		index1 = i + nxp*(j-1)	
	        Do k=1,Num_family
	   	   U(index1,3*k-2)=U(index1,3*k-2)+deltat*
     1           dUdt(index1,3*k-2)/(Sew(i+1)*Sns(j+1))+deltat*
     2           dAggdt(index1,3*k-2)+deltat*Source(index1,3*k-2)
	 	If (Flag3 .LT. 3 ) then 
	   	   U(index1,3*k-1)=(U(index1,3*k-1)+deltat*
     1             dUdt(index1,3*k-1)/(Sew(i+1)*Sns(j+1))+deltat*
     2             dAggdt(index1,3*k-1)+deltat*Source(index1,3*k-1))/
     3 		   (1.0+Drag_t(k)*deltat)
	   	   U(index1,3*k)=(U(index1,3*k)+deltat*dUdt(index1,3*k)/
     1             (Sew(i+1)*Sns(j+1))+deltat*dAggdt(index1,3*k)+
     2             deltat*Source(index1,3*k))/(1.0+Drag_t(k)*deltat)
		else
		   U_old(index1,k)=U(index1,3*k-1);
		   V_old(index1,k)=U(index1,3*k);
		End If
		End Do
	      End Do
	   End Do
	End If
C---------------------------------------------------------------------------
C	! Calculating the volume fraction alpha
	  alpha(:)=0.0
	  Do i=1,nxp
	    Do j=1,nyp
	      index1 = i + nxp*(j-1)
	      Do k=1,Num_family
	       alpha(index1)=alpha(index1)+U(index1,3*k-2)/(RBC_dens)
	      End Do
	    End Do
	  End Do
C        !Solving Fluid 	  
	   If (Flag3 .GT. 1) then 
		Call get_fluid_Velocity
	   End If
C 	! Updating the time
	   current_time  = current_time+deltat
C	! Writing the results and updating the save time
           If (current_time .GE. Save_time) then
	        write (*,*) 'current time=', current_time, 'dt=',deltat
 		Call Save_data
	        Save_time = Save_time+t_final/15	
           End if
	End Do
	Close (unit=1)
	Close (unit=3)
	End Program FV_PBE_lax_fred
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	  Opening the mesh file   	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

	Subroutine read_mesh
	Use List_of_Variables
	 Open(5,File='_Mesh.Dat') 
	 Read(5,*)nx,ny,Dxg,Dyg
C        Fluid parameters
	 Allocate( Dxep (nx)  )
	 Allocate( Dxpw (nx)  )
	 Allocate( Dynp(ny)   )
	 Allocate( Dyps(ny)   )
	 Allocate( Sns(ny)    )
	 Allocate( Sew(nx)    )
	 Allocate( Xu(nx)     )
	 Allocate( Yv(ny)     )
	 Allocate( Dxepu(nx)  )
	 Allocate( Dxpwu(nx)  )
	 Allocate( Sewu (nx)  )
	 Allocate( Dynpv(ny)  )
	 Allocate( Dypsv(ny)  )
	 Allocate( Snsv(ny)   )
	 Allocate( X(nx)      )
	 Allocate( Y(ny)      )
	 Allocate( Du(nx,ny)  )    
	 Allocate( Dv(nx,ny)  )    
	 Allocate( Ppp(nx,ny) )     
	 Allocate( Ap(nx,ny)  )    
	 Allocate( An(nx,ny)  )    
	 Allocate( As(nx,ny)  )    
	 Allocate( Ae(nx,ny)  )    
	 Allocate( Aw(nx,ny)  )    
	 Allocate( Su(nx,ny)  )    
	 Allocate( Sp(nx,ny)  )
	 Allocate( Ug(nx,ny)  ) 
	 Allocate( Vg(nx,ny)  )
	 Allocate( P(nx,ny)   ) 
	 Allocate( PP(nx,ny)  )
	 Allocate( Densit(nx,ny)    )
	 Allocate( U_fluid(nx,ny)   ) 
	 Allocate( V_fluid(nx,ny)   )     
	 Allocate( Densit_old(nx,ny)  ) 
	 Allocate( U_fluid_old(nx,ny) ) 
	 Allocate( V_fluid_old(nx,ny) ) 
	 Allocate( Up(nx,ny,Num_family)  ) 
	 Allocate( Vp(nx,ny,Num_family)  )
	 Allocate( Den_p(nx,ny,Num_family)  )
C        Particle parameters
	 nxp=nx-2
	 nyp=ny-2
	 Allocate( pressure(nxp*nyp)  )  
	 Allocate( alpha(nxp*nyp)     ) 
	 Allocate( Ugg(nxp,nyp)  ) 
	 Allocate( Vgg(nxp,nyp)  )	 
	 Allocate( U(nxp*nyp,Num_family*3)    )  
	 Allocate( dUdt(nxp*nyp,Num_family*3) )
	 Allocate( dUdx(nxp*nyp,Num_family*3) )
	 Allocate( dUdy(nxp*nyp,Num_family*3) )
	 Allocate( U_old(nxp*nyp,Num_family)  )    
	 Allocate( V_old(nxp*nyp,Num_family)  )
	 Allocate( Source(nxp*nyp,Num_family*3)   )  
	 Allocate( dAggdt(nxp*nyp,Num_family*3)   )  
	 Allocate( di_pressure_di_alpha(nxp*nyp)  ) 

C       read mesh and scale it using H and L
	Do I=1,nx;
		Read(5,*)X(I)
		X(I)=X(I)*L
	End Do
	Do J=1,ny
		Read(5,*)Y(J)
		Y(J)=Y(J)*H
	End Do
	Close(Unit=5, Status='Keep')

C       Defining used dx and dy 
	! ------------- General
	Dxpw(1)=0.0
	Dxep(nx)=0.0
	
	Do I=1,nx-1
	    Dxep(I)=X(I+1)-X(I)
	    Dxpw(I+1)=Dxep(I)
	End Do

	Dyps(1)=0.0
	Dynp(ny)=0.0
	Do J=1,ny-1
		Dynp(J)=Y(J+1)-Y(J)
		Dyps(J+1)=Dynp(J)
	End Do

	Sew(1)=0.0        !Not Found
	Sew(nx)=0.0       !Not Found
	Do I=2,ny-1
		Sew(I)=0.5*(Dxep(I)+Dxpw(I))
	End Do
		
	Sns(1)=0.0        !Not Found
	Sns(ny)=0.0       !Not Found
	Do J=2,ny-1
		Sns(J)=0.5*(Dynp(J)+Dyps(J))
	End Do

	!------------------ U
	Xu(1)=0.0         !Not Found
	Xu(2)=0.0         !Not Found
	Do I=3,nx
		Xu(I)=0.5*(X(I)+X(I-1))
	End Do
	
	Dxpwu(1)=0.0  !Not Found
	Dxpwu(2)=0.0  !Not Found
	Dxepu(1)=0.0  !Not Found
	Dxepu(nx)=0.0 !Not Found
	Do I=2,nx-1
		Dxepu(I)=Xu(I+1)-Xu(I)
		Dxpwu(I+1)=Dxepu(I)
	End Do

	Sewu(1)=0.0       !Not Found
	Sewu(2)=0.0       !Is Boundary Condition  U=0
	Do I=3,nx-1
		Sewu(I)=0.5*(Dxepu(I)+Dxpwu(I))
	End Do

	!------------------ V
	Yv(1)=0.0               !Not Found
	Yv(2)=0.0               !Not Found
	Do J=3,ny
		Yv(J)=0.5*(Y(J)+Y(J-1))
	End Do

	Dypsv(1)=0.0    !Not Found
	Dypsv(2)=0.0    !Not Found
	Dynpv(ny)=0.0   !Not Found
	Do J=2,ny-1
		Dynpv(J)=Yv(J+1)-Yv(J)
		Dypsv(J+1)=Dynpv(J)
	End Do

	Snsv(1)=0.0             !Not Found
	Snsv(ny)=0.0    !Is Boundary Condition V=0
	Do J=3,ny-1
		Snsv(J)=0.5*(Dynpv(J)+Dypsv(J))
	End Do

	End Subroutine read_mesh

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	    Initialization      	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	
	Subroutine Set_Initial_Condition  
	Use List_of_Variables
	alpha=0.0;
	Do k=1,Num_family
		vol_fam(k)=Vol_RBC*k
		mass_fam(k)=vol_fam(k)*RBC_dens
		radius_fam(k)=(vol_fam(k)*3.0/(4.0*pi))**(1.0/3.0)
	End Do
	Do i = 1,nxp
	    Do j = 1,nyp
		index1 = i + nxp*(j-1)
		U(index1,:)=0.0
	    End Do
	End Do
C	Making B_Agg and D_Agg
	B_Agg=0.0
	D_Agg=0.0
	if (Flag2 .EQ. 2) then 
	   k=1
	   Do i=1,Num_family-1
	      Do j=1,Num_family-1
	         If (j .GE. i) then
		     B_Agg(i,j)=unknown(k) !parameter unknown is coming from experiment look inside the list_of_variable.f file to find out these numbers
		     k=k+1
	         End If	
      End Do
	   End Do
	   Do i=1,Num_family-1
	      Do j=1,Num_family-1
	         If (i .GT. j) then
	            ! Update if the volume relation changed
		    B_Agg(i,j)=1.0/(2.0**(Real(i-j)))*B_Agg(j,i) 
	         End If
	      End Do
	    End Do
	    Do i=2,Num_family
	       Do j=1,i-1
	          D_Agg(i,j)=B_Agg(i-1,j)
	       End Do
	    End Do
	End If

C	Fluid initialization
	Do i = 1,nx
	    Do j = 1,ny

		If (Flag4 .EQ. 2) then
		V_fluid_old(i,j)=Uin
		U_fluid_old(i,j)=0.0
		U_fluid(i,j)=0.0
		V_fluid(i,j)=Uin
		elseif (Flag4 .EQ. 3) then
		U_fluid_old(i,j)=Uin
		V_fluid_old(i,j)=0.0
		U_fluid(i,j)=Uin
		V_fluid(i,j)=0.0
		else
		U_fluid_old(i,j)=0.0
		V_fluid_old(i,j)=0.0
		U_fluid(i,j)=0.0
		V_fluid(i,j)=0.0
		End If
		Pp(i,j)=0.0
		P(i,j)=0.0
		Ap(i,j)=0.0
		An(i,j)=0.0
		As(i,j)=0.0
		Ae(i,j)=0.0
		Aw(i,j)=0.0
		Sp(i,j)=0.0
		Su(i,j)=0.0
		Dv(i,j)=0.0
		Du(i,j)=0.0
	    End Do
	End Do
	Up(:,:,:)=0.0
	Vp(:,:,:)=0.0
	Den_p(:,:,:)=0.0
	
	Sormax =0.00001


	End Subroutine Set_Initial_Condition 

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	    Test Condition       	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	
	Subroutine Set_test_Condition 
	Use List_of_Variables
	If (Flag .EQ. 1) then    ! Test 1 Particles in Center
	   CFL=0.7
	   epsilon1=0.000000000001
 	   Do i = nxp/2-20,nxp/2+20
	      Do j =nyp/2-20, nyp/2+20
	     	index1 = i + nxp*(j-1)
	        Do k=1,Num_family
		   U(index1,3*k-2) = 0.02*RBC_dens 
		   U(index1,3*k-1:3*k) = 0.0
		End Do
	      End Do
	   End Do
	Else If (Flag .EQ. 2) then   ! Test 2 Full Domain of Particles
	   epsilon1=0.000000000001
	   CFL=0.7
 	   Do i = 1,nxp
	      Do j =1, nyp
	     	index1 = i + nxp*(j-1)
	        Do k=1,Num_family
		   U(index1,3*k-2) = 0.02*RBC_dens 
		   U(index1,3*k-1:3*k) = 0.0
		End Do
	      End Do
	   End Do
	Else If (Flag .GT. 2) then !test 3 No particle
	   epsilon1=0.000000000001
	   CFL=0.5
 	   Do i = 1,nxp
	      Do j =1, nyp
	     	index1 = i + nxp*(j-1)
		   U(index1,:) = 0.0 
	      End Do
	   End Do
	End if
	Do i=1,nxp
	  Do j=1,nyp
	    index1 = i + nxp*(j-1)
	    Do k=1,Num_family
	    alpha(index1)=alpha(index1)+U(index1,3*k-2)/(RBC_dens)
	    End Do
	  End Do
	End Do
	Densit (:,:)=ro_fluid
        Do I=1,nxp
           Do J=1,nyp
	      index1 = I + nxp*(J-1)
	      Densit_old(I+1,J+1)=(1.0-alpha(index1))*ro_fluid	
	   End Do
 	End Do
        Do I=1,nx
	       Densit_old(I,1)=Densit_old(I,2)
	       Densit_old(I,ny)=Densit_old(I,ny-1)
 	End Do
        Do J=1,ny
	       Densit_old(1,J)=Densit_old(2,J)
	       Densit_old(nx,J)=Densit_old(nx-1,J)
 	End Do
	Densit_old(nx,ny)=Densit_old(nx-1,ny-1)
        Densit_old(nx,1)=Densit_old(nx-1,2)
	Densit_old(1,1)=Densit_old(2,2)
        Densit_old(1,ny)=Densit_old(2,ny-1)

	
	End Subroutine Set_test_Condition

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	    Lax Friedrichs Method      	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	
	Subroutine Local_Lax_Friedrichs_x
	Use List_of_Variables
	If (Ul(1) .LT. 0.01) then
	d_pl=0.0
	End If
	If (Ur(1) .LT. 0.01) then
	d_pr=0.0
	End If
	fastest_wave = MAX(abs(Ul(2))/(Ul(1)+epsilon1)
     2	+ sqrt(d_pl/(Ul(1)+epsilon1)), abs(Ur(2))/(Ur(1)+epsilon1)
     3	+ sqrt(d_pr/(Ur(1)+epsilon1)), abs(Ul(2))/(Ul(1)+epsilon1))

 	Fl(1) = Ul(2)
 	Fl(2) = Ul(2)/(Ul(1)+epsilon1)*Ul(2)+pl
 	Fl(3) = Ul(2)/(Ul(1)+epsilon1)*Ul(3)        
 	Fr(1) = Ur(2)
 	Fr(2) = Ur(2)/(Ur(1)+epsilon1)*Ur(2)+pr
 	Fr(3) = Ur(2)/(Ur(1)+epsilon1)*Ur(3) 
        Fx = 0.5*(Fl+Fr) - 0.5*fastest_wave*(Ur-Ul)
	End  Subroutine Local_Lax_Friedrichs_x
	
	Subroutine Local_Lax_Friedrichs_y
	Use List_of_Variables
	If (Uu(1) .LT. 0.01) then
	d_pu=0.0
	End If
	If (Ud(1) .LT. 0.01) then
	d_pd=0.0
	End If
	fastest_wave = max(abs(Ud(3))/(Ud(1)+epsilon1)
     2  + sqrt(d_pd/(Ud(1)+epsilon1)), abs(Uu(3))/(Uu(1)+epsilon1)
     3  + sqrt(d_pu/(Uu(1)+epsilon1)), abs(Ud(3))/(Ud(1)+epsilon1))
 	Fd(1) = Ud(3)
 	Fd(2) = Ud(3)/(Ud(1)+epsilon1)*Ud(2)
 	Fd(3) = Ud(3)/(Ud(1)+epsilon1)*Ud(3)+pd        
 	Fu(1) = Uu(3)
 	Fu(2) = Uu(3)/(Uu(1)+epsilon1)*Uu(2)
 	Fu(3) = Uu(3)/(Uu(1)+epsilon1)*Uu(3)+pu
        Fy = 0.5*(Fu+Fd) - 0.5*fastest_wave*(Uu-Ud)
	End  Subroutine Local_Lax_Friedrichs_y

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	       minmod method       	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
c	Subroutine minmod
c	Use List_of_Variables
c	  d_ud_loc(:)=0.0	
c	  if (fdiff(1)*bdiff(1) .Lt. 0.0) then
c	   	d_ud_loc(1) = 0.0
c	  elseif (abs(fdiff(1)) .Lt. abs(bdiff(1))) then
c		d_ud_loc(1) = fdiff(1)
c	  else
c		d_ud_loc(1) = bdiff(1)
c 	  end if
c
c	  if (fdiff(2)*bdiff(2) .Lt. 0.0) then
c   		d_ud_loc(2) = 0.0
c  	  elseif (abs(fdiff(2)) .Lt. abs(bdiff(2))) then
c    		d_ud_loc(2) = fdiff(2)
c  	  else
c   		d_ud_loc(2) = bdiff(2)
c  	  end if
c
c  	  if (fdiff(3)*bdiff(3) .Lt. 0.0) then
c    		d_ud_loc(3) = 0.0
c  	  elseif (abs(fdiff(3)) .Lt. abs(bdiff(3))) then
c    	        d_ud_loc(3) = fdiff(3)
c  	  else
c    	        d_ud_loc(3) = bdiff(3)
c  	  end if
		
c	End  Subroutine minmod	
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	    source term in PBE      	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
C	Here the source term in PBE will be calculated
	Subroutine Get_Source
	Use List_of_Variables
	If (Flag1 .EQ. 1) Then   ! Just Gravity
	Drag_t(:)=0.0   
	Do i=1,nxp
	      Do j=1,nyp
		 index1 = i + nxp*(j-1)
		 Do k=1,Num_family
		    Source(index1,3*k-2)=0.0
		    Source(index1,3*k-1)=0.0
		    Source(index1,3*k)=U(index1,3*k-2)*Gravity
		 End Do
	      End Do
	   End Do
	Else If (Flag1 .EQ. 2) Then ! Gravity and Buoyancy
	Drag_t(:)=0.0   
	Do i=1,nxp
	      Do j=1,nyp
		 index1 = i + nxp*(j-1)
		 Do k=1,Num_family
		    Source(index1,3*k-2)=0.0
		    Source(index1,3*k-1)=0.0
		    Source(index1,3*k)=(U(index1,3*k-2)-alpha(index1)*
     1              ro_fluid)*Gravity
		 End Do
	      End Do
	   End Do
 	Else If (Flag1 .EQ. 3) Then ! with Drag no buoyancy
           Do k=1,Num_family
	      Drag_t(k)=6.0*pi*miu*radius_fam(k)/mass_fam(k)
	   End Do
	   Do i=1,nxp
	      Do j=1,nyp
		 index1 = i + nxp*(j-1)
		 Do k=1,Num_family	    
		    Source(index1,3*k-2)=0.0
		    Source(index1,3*k-1)=0.0
	            Source(index1,3*k)=U(index1,3*k-2)*Gravity !alpha*(U(index1,1)-ro_fluid)  
		 End Do
	      End Do
	   End Do
 	Else If (Flag1 .EQ. 4) Then ! with Drag + buoyancy
           Do k=1,Num_family
	      Drag_t(k)=6.0*pi*miu*radius_fam(k)/mass_fam(k)
	   End Do
	   Do i=1,nxp
	      Do j=1,nyp
		 index1 = i + nxp*(j-1)
		 Do k=1,Num_family	   
		    Source(index1,3*k-2)=0.0
		    Source(index1,3*k-1)=0.0
	            Source(index1,3*k)=(U(index1,3*k-2)-alpha(index1)*
     1		    ro_fluid)*Gravity  
		 End Do
	      End Do
	   End Do
 	Else If (Flag1 .EQ. 5) Then ! with Drag + buoyancy + lift later
           Do k=1,Num_family
	      Drag_t(k)=6.0*pi*miu*radius_fam(k)/mass_fam(k)
	   End Do
	   Do i=1,nxp
	      Do j=1,nyp
		 index1 = i + nxp*(j-1)
		 Do k=1,Num_family	   
		    Source(index1,3*k-2)=0.0
		    Source(index1,3*k-1)=0.0
	            Source(index1,3*k)=(U(index1,3*k-2)-alpha(index1)*
     1  	    ro_fluid)*Gravity  
		 End Do
	      End Do
	   End Do
	Else If (Flag1 .EQ. 6) then
	Drag_t=0.0   
	Do i=1,nxp
	      Do j=1,nyp
		 index1 = i + nxp*(j-1)
		 Do k=1,Num_family
		    Source(index1,3*k-2)=0.0
		    Source(index1,3*k-1)=0.0
		    Source(index1,3*k)=0.0
		 End Do
	      End Do
	   End Do 
	End If
	End Subroutine Get_Source
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	      aggregation         	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
C	here the calculation for aggregation is done
	Subroutine Aggregate
	Use List_of_Variables
	dAggdt=0.0
C	Calculating dAggdt
	Do k=1,Num_family
	   If (k .EQ. 1) then
	      Do i=1,Num_family
	          dAggdt(:,3*k-2)=dAggdt(:,3*k-2)-1.0*
     1            B_Agg(k,i)*U(:,3*k-2)*U(:,3*i-2)
	          dAggdt(:,3*k-1)=dAggdt(:,3*k-1)-1.0*
     1            B_Agg(k,i)*U(:,3*k-2)*U(:,3*i-1)
	          dAggdt(:,3*k)=dAggdt(:,3*k)-1.0*
     1            B_Agg(k,i)*U(:,3*k-2)*U(:,3*i)
	      End Do
	   Else
	      Do i=1,Num_family
		dAggdt(:,3*k-2)=dAggdt(:,3*k-2)+D_Agg(k,i)*
     1	         U(:,3*k-5)*U(:,3*i-2)-B_Agg(k,i)*U(:,3*k-2)*U(:,3*i-2)
		dAggdt(:,3*k-1)=dAggdt(:,3*k-1)+D_Agg(k,i)*
     1	         U(:,3*k-5)*U(:,3*i-1)-B_Agg(k,i)*U(:,3*k-2)*U(:,3*i-1)
		dAggdt(:,3*k)=dAggdt(:,3*k-2)+D_Agg(k,i)*
     1	         U(:,3*k-5)*U(:,3*i)-B_Agg(k,i)*U(:,3*k-2)*U(:,3*i)
	      End Do
           End If
        End Do
	If (Flag2 .EQ. 3) then
		! ADD DisAGG
	End If

	End  Subroutine Aggregate

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	    Fluid Flow Solver      	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
C	Here velocities in x and y direction will be calculated
	Subroutine get_fluid_Velocity  ! (X,Y,deltax,deltay,deltat,U_fluid,V_fluid,Pp,source)
	Use List_of_Variables	

	 Kr=0
	 Source1=10.0
         Do I=1,nxp
            Do J=1,nyp
	        index1 = I + nxp*(J-1)
	        Densit(I+1,J+1)=(1.0-alpha(index1))*ro_fluid
		Do k=1,Num_family
		  Den_p(I+1,J+1,k)=U(index1,3*k-2)
	       	  Up(I+1,J+1,k)=U(index1,3*k-1)
	       	  Vp(I+1,J+1,k)=U(index1,3*k)
		End Do
	    End Do
 	 End Do
        Do I=1,nx
	       Densit(I,1)=Densit(I,2)
	       Densit(I,ny)=Densit(I,ny-1)
	       Up(I,1,:)=Up(I,2,:)
	       Vp(I,1,:)=Vp(I,2,:)
	       Den_p(I,1,:)=Den_p(I,2,:)
	       Up(I,ny,:)=Up(I,ny-1,:)
	       Vp(I,ny,:)=Vp(I,ny-1,:)
	       Den_p(I,ny,:)=Den_p(I,ny-1,:)
 	End Do
        Do J=1,ny
	       Densit(1,J)=Densit(2,J)
	       Densit(nx,J)=Densit(nx-1,J)
	       Up(1,J,:)=Up(2,J,:)
	       Vp(1,J,:)=Vp(2,J,:)
	       Den_p(1,J,:)=Den_p(2,J,:)
	       Up(nx,J,:)=Up(nx-1,J,:)
	       Vp(nx,J,:)=Vp(nx-1,J,:)
	       Den_p(nx,J,:)=Den_p(nx-1,J,:)
 	End Do
	Densit(nx,ny)=Densit(nx-1,ny-1)
        Densit(nx,1)=Densit(nx-1,2)
	Densit(1,1)=Densit(2,2)
        Densit(1,ny)=Densit(2,ny-1)


	 Do While (Source1 .GT. Sormax)
           Kr=Kr+1 
	if (Flag3 .EQ. 3 .AND. Flag1 .GT. 2) then
	   Do i=1,nxp
	       Do j=1,nyp
		 index1 = i + nxp*(j-1)
	         Ugg(i,j)=(U_fluid(i+1,j+1)+U_fluid(i+2,j+1))/2.0
		 Vgg(i,j)=(V_fluid(i+1,j+1)+V_fluid(i+1,j+2))/2.0

		 Do k=1,Num_family
	   U(index1,3*k-1)=(U_old(index1,k)+deltat*dUdt(index1,3*k-1)/
     1     (Sew(i+1)*Sns(j+1))+deltat*dAggdt(index1,3*k-1)+deltat*
     2     Source(index1,3*k-1)-Drag_t(k)*deltat*U(index1,3*k-2)
     3     *Ugg(i,j))/(1.0+Drag_t(k)*deltat)
	   U(index1,3*k)=(V_old(index1,k)+deltat*dUdt(index1,3*k)/
     1     (Sew(i+1)*Sns(j+1))+deltat*dAggdt(index1,3*k)+
     2     deltat*Source(index1,3*k)-Drag_t(k)*deltat*U(index1,3*k-2)
     3     *Vgg(i,j))/(1.0+Drag_t(k)*deltat)
                  End Do
	      End Do
	   End Do
	End If

C   ------------------Calc U ---------------------    		
           Do I=3,nx-1
              Do J=2,ny-1
                    Asu=Sewu(I)
                    Awe=Yv(J+1)-Yv(J) 
                    Aea=Awe
                    Ano=Asu
                    Vol=Asu*Awe
C                  -----Calculate Convection Coefficients--------------  
               Ce=Densit(I,J)*.5*(U_fluid(I,J)+U_fluid(I+1,J))*Aea
               Cw=Densit(I-1,J)*.5*(U_fluid(I,J)+U_fluid(I-1,J))*Awe
               Cs=0.25*(Densit(I,J)+Densit(I,J-1)+Densit(I-1,J)+
     1         Densit(I-1,J-1))*.5*(V_fluid(I,J)+V_fluid(I-1,J))*Asu
               Cn=0.25*(Densit(I,J)+Densit(I,J+1)+Densit(I-1,J)+
     1         Densit(I-1,J+1))*.5*(V_fluid(I,J+1)+V_fluid(I-1,J+1))*Ano
C                 -----Calculate Diffusion Coefficients
                  Dn=Viscos*Ano/(Y(J+1)-Y(J))
                  Ds=Viscos*Asu/(Y(J)-Y(J-1))
                  De=Viscos*(Aea/(Dxepu(I)))
                  Dw=Viscos*(Awe/(Dxpwu(I)))
C                 -----Assemble Main Coefficients
                  An(I,J)= max(-Cn, (Dn-Cn/2.0), 0.0)
                  As(I,J)= max( Cs, (Ds+Cs/2.0), 0.0)
                  Ae(I,J)= max(-Ce, (De-Ce/2.0), 0.0)
                  Aw(I,J)= max( Cw, (Dw+Cw/2.0), 0.0)
C                 -----Calculate Coefficients Of Source Terms
		  Gue=(U_fluid(I,J)-U_fluid(I+1,J))/(Xu(I)-Xu(I+1))
                  Guw=(U_fluid(I,J)-U_fluid(I-1,J))/(Xu(I)-Xu(I-1))
                  Gvs=(V_fluid(I,J)-V_fluid(I-1,J))/(X(I)-X(I-1))
                  Gvn=(V_fluid(I,J+1)-V_fluid(I-1,J+1))/(X(I)-X(I-1))

                  Gen=(Viscos*Gue-Viscos*Guw)*Awe+(Viscos*Gvn-
     1            Viscos*Gvs)*Ano

                  Spp=-Awe*(P(I,J)-P(I-1,J))

                  Smp=Cn-Cs+Ce-Cw
                  Cp=max(0.0,Smp)
                  Sp(I,J)=-Cp
                  Su(I,J)=Cp*U_fluid(I,J)+Spp
                  Su(I,J)=Su(I,J)+Gen				  
	       End Do
	    End Do
C        %%%%%%%%%% BOUNDARY U %%%%%%%%% (learn more)
C 	    TOP WALL
            J=ny-1
	    Do I=2,nx-1
		Asu=Sewu(I)
                U_fluid(I,J)=0.0 ! check later
		Tmult=Viscos*Asu/(Yv(J+1)-Y(J))
                Sp(I,J)=Sp(I,J)-Tmult
                Su(I,J)=Su(I,J)+Tmult*U_fluid(I,J+1)
                U_fluid(I,J+1)=U_fluid(I,J)
  		An(I,J)=0.0
            End Do
C  	    SOUTH WALL
            J=2
	    Do I=2,nx-1
		Asu=Sewu(I)
                U_fluid(I,J)=0.0
                Tmult=Viscos*Asu/(Y(J)-Yv(J))
                Sp(I,J)=Sp(I,J)-Tmult
                Su(I,J)=Su(I,J)+Tmult*U_fluid(I,J-1)
	        As(I,J)=0.0
		U_fluid(I,J-1)=U_fluid(I,J)
            End Do
C           East WALL
            I=nx-1
	    Do J=2,ny-1
                U_fluid(I,J)=0.0
                U_fluid(nx,J)=U_fluid(I,J)
	        Ae(I,J)=0.0
            End Do
C  	    West WALL
            I=2
	    Do J=2,ny-1
               U_fluid(I,J)=0.0
	       Aw(I,J)=0.0
            End Do
	
            Resoru=0.0
            Do I=3,nx-1
               Do J=2,ny-1
                  Vol=Sewu(I)*(Yv(J+1)-Yv(J))
                  Ap(I,J)=An(I,J)+As(I,J)+Ae(I,J)+Aw(I,J)-Sp(I,J)
C                 ----- transient 
                  Su(I,J)=Su(I,J)+0.5*(Densit_old(I-1,J)+
     1            Densit_old(I,J))*Vol/deltat*U_fluid_old(I,J)
             Ap(I,J)=Ap(I,J)+0.5*(Densit(I-1,J)+Densit(I,J))*Vol/deltat
	          if (Flag3 .EQ. 3) then
		     F1=0
		     F2=0
	   	     Do k=1,Num_family
			F1=F1+Drag_t(k)*0.5*(Den_p(I,J,k)*
     1        Up(I,J,k)+Den_p(I+1,J,k)*Up(I+1,J,k))
			F2=F2+Drag_t(k)*0.5*(Den_p(I,J,k)
     1                  +Den_p(I+1,J,k))
                     End Do
			Su(I,J)=Su(I,J)+Sew(I)*Sns(J)*deltat*F1
			Ap(I,J)=Ap(I,J)+Sew(I)*Sns(J)*deltat*F2
		   End If
                  Du(I,J)=Ap(I,J)
                  Resor=An(I,J)*U_fluid(I,J+1)+As(I,J)*U_fluid(I,J-1)+
     1            Ae(I,J)*U_fluid(I+1,J)+Aw(I,J)*U_fluid(I-1,J)-
     2            Ap(I,J)*U_fluid(I,J)+Su(I,J)
                  Resoru=Resoru+abs(Resor)

C               --- Under-Relaxation
                  Ap(I,J)=Ap(I,J)/Urfu
                  Su(I,J)=Su(I,J)+(1.0-Urfu)*Ap(I,J)*U_fluid(I,J)
C                 Du(I,J)=Du(I,J)*Urfu !Simple
                  Du(I,J)=Du(I,J)*(1.0-Urfu)/Urfu-Sp(I,J)	!Simplec
               End Do
             End Do


	     Do  Nu=1,NswpU
		Call Liner_Solver(3,2,U_fluid)
	     End Do

C   ------------------ End Calc U ---------------------  

C   ------------------ Calc V ---------------------  
           Do I=2,nx-1
              Do J=3,ny-1
                  Awe=Snsv(J);
                  Asu=Xu(I+1)-Xu(I);
                  Aea=Awe;
                  Ano=Asu;
                  Vol=Snsv(J)*(Xu(I+1)-Xu(I))
                Ce=0.25*(Densit(I,J)+Densit(I+1,J)+Densit(I,J-1)+
     1    Densit(I+1,J-1))*0.5*(U_fluid(I+1,J)+U_fluid(I+1,J-1))*Aea
                Cw=0.25*(Densit(I,J)+Densit(I-1,J)+Densit(I,J-1)+ 
     1    Densit(I-1,J-1))*0.5*(U_fluid(I,J)+U_fluid(I,J-1))*Awe
                Cn=Densit(I,J)*0.5*(V_fluid(I,J)+V_fluid(I,J+1))*Ano
                Cs=Densit(I,J-1)*0.5*(V_fluid(I,J)+V_fluid(I,J-1))*Asu

                Dn=Viscos*Ano/Dynpv(J)
                Ds=Viscos*Asu/Dypsv(J)
                De=Viscos*Aea/(X(I+1)-X(I))
                Dw=Viscos*Awe/(X(I)-X(I-1))

               Spp=-Ano*(P(I,J)-P(I,J-1))!+Gravity*.5*(Densit(I,J)+Densit(I,J-1))*Vol
                An(I,J)= max(-Cn, (Dn-Cn/2.0), 0.0)
                As(I,J)= max( Cs, (Ds+Cs/2.0), 0.0)
                Ae(I,J)= max(-Ce, (De-Ce/2.0), 0.0)
                Aw(I,J)= max( Cw, (Dw+Cw/2.0), 0.0)
                Smp=Cn-Cs+Ce-Cw

                Guw=(U_fluid(I,J)-U_fluid(I,J-1))/(Y(J)-Y(J-1))
                Gvn=(V_fluid(I,J)-V_fluid(I,J+1))/(Yv(J)-Yv(J+1))
                Gvs=(V_fluid(I,J)-V_fluid(I,J-1))/(Yv(J)-Yv(J-1))
                Gue=(U_fluid(I+1,J)-U_fluid(I+1,J-1))/(Y(J)-Y(J-1))

	        Gen=0.0 !Gravity*.5*(Densit(I,J)+Densit(I,J-1))*Vol;
	        Gen=Gen+(Viscos*Gue-Viscos*Guw)*Awe+(Viscos*Gvn
     1      -Viscos*Gvs)*Ano
                Cp=max(0.0,Smp)
                Sp(I,J)=-Cp
                Su(I,J)=Cp*V_fluid(I,J)+Spp
               Su(I,J)=Su(I,J)+Gen
              End Do
           End Do
C        %%%%%%%%%% BOUNDARY V %%%%%%%%%
	!-----North
	   J=ny-1
           Do I=2,nx-1
              An(I,J)=0.0
              V_fluid(I,J)=0.0
           End Do
    	!-----South 
	   J=2
           Do I=2,nx-1
              As(I,J)=0.0
              V_fluid(I,J)=0.0
           End Do
        !-----West
	    I=2
            Do J=2,ny-1
                 Awe=Snsv(J)
                 V_fluid(I,J)=0.0
                 V_fluid(I-1,J)=V_fluid(I,J)
                 Tmult=Viscos*Awe/(X(I)-Xu(I))
                 Sp(I,J)=Sp(I,J)-Tmult
                 Su(I,J)=Su(I,J)+Tmult*V_fluid(I-1,J)
                 Aw(I,J)=0.0
            End Do
        !-----Eest
	    I=nx-1
            Do J=2,ny-1
                 Awe=Snsv(J)
                 V_fluid(I,J)=0.0
                 V_fluid(I+1,J)=V_fluid(I,J)
                 Tmult=Viscos*Awe/(Xu(I+1)-X(I))
                 Sp(I,J)=Sp(I,J)-Tmult
                 Su(I,J)=Su(I,J)+Tmult*V_fluid(I+1,J)
                 Ae(I,J)=0.0
            End Do
            Resorv=0.0

            Do I=2,nx-1
               Do J=3,ny-1
                  Vol=Snsv(J)*(Xu(I+1)-Xu(I))
                  Ap(I,J)=An(I,J)+As(I,J)+Ae(I,J)+Aw(I,J)-Sp(I,J)
C        ! Transient
               Su(I,J)=Su(I,J)+0.5*(Densit_old(I,J-1)+Densit_old(I,J))*
     1          Vol/deltat*V_fluid_old(I,J)

             Ap(I,J)=Ap(I,J)+0.5*(Densit(I,J-1)+Densit(I,J))*Vol/deltat
	          if (Flag3 .EQ. 3) then
		     F1=0
		     F2=0
 	   	     Do k=1,Num_family
		        F1=F1+Drag_t(k)*0.5*(Den_p(I,J,k)*Vp(I,J,k)+
     1                  Den_p(I,J+1,k)*Vp(I,J+1,k))
		        F2=F2+Drag_t(k)*0.5*(Den_p(I,J,k)
     1                  +Den_p(I,J+1,k))
                     End Do
			Su(I,J)=Su(I,J)+Sew(I)*Sns(J)*deltat*F1
			Ap(I,J)=Ap(I,J)+Sew(I)*Sns(J)*deltat*F2
		  End If

                  Dv(I,J)=Ap(I,J)

                  Resor=An(I,J)*V_fluid(I,J+1)+As(I,J)*V_fluid(I,J-1)+
     1            Ae(I,J)*V_fluid(I+1,J)+Aw(I,J)*V_fluid(I-1,J)-
     2            Ap(I,J)*V_fluid(I,J)+Su(I,J)
                  Resorv=Resorv+abs(Resor)
                  Ap(I,J)=Ap(I,J)/Urfv
                  Su(I,J)=Su(I,J)+(1.0-Urfv)*Ap(I,J)*V_fluid(I,J)
C                 Dv(I,J)=Dv(I,J)*Urfv !Simple
                  Dv(I,J)=Dv(I,J)*(1.0-Urfv)/Urfv-Sp(I,J) !Simplec
               End Do
            End Do

	
	    Do  Nv=1,NswpV 
		Call Liner_Solver(2,3,V_fluid)
	    End Do
C   ------------------ Endof Calc V ---------------------  

C   ------------------ Calc P --------------------- (learn more) 	
            Resorm=0.0
            Do I=2,nx-1
               Do J=2,ny-1
                   Aea=Sns(J)
                   Awe=Aea
                   Ano=Sew(I)
                   Asu=Ano
                   Vol=Sns(J)*Sew(I)
                !----- North
                  if (J .EQ. ny-1) then
                        An(I,J)=0.0
                        G2starn=0.0
                  else
         An(I,J)=-0.5*(Densit(I,J)+Densit(I,J+1))*(Ano)**2.0/Dv(I,J+1)
            G2starn=0.5*(Densit(I,J)+Densit(I,J+1))*V_fluid(I,J+1)*Ano
		  end if

                !----- South
                   if (J .EQ. 2) then
                        As(I,J)=0.0
                        G2stars=0.0
                    else
           As(I,J)=-0.5*(Densit(I,J)+Densit(I,J-1))*(Asu)**2.0/Dv(I,J)
              G2stars=0.5*(Densit(I,J)+Densit(I,J-1))*V_fluid(I,J)*Asu
		    end if
                !----- East
                    if (I .EQ. nx-1) then
                       Ae(I,J)=0.0
                        G1stare=0.0
                    else
          Ae(I,J)=-0.5*(Densit(I,J)+Densit(I+1,J))*(Awe)**2.0/Du(I+1,J)
             G1stare=0.5*(Densit(I,J)+Densit(I+1,J))*U_fluid(I+1,J)*Aea
		    end if
                !----- West
                    if (I .EQ. 2) then
                       Aw(I,J)=0.0
                        G1starw=0.0
                    else
           Aw(I,J)=-0.5*(Densit(I,J)+Densit(I-1,J))*(Aea)**2.0/Du(I,J)
              G1starw=0.5*(Densit(I,J)+Densit(I-1,J))*U_fluid(I,J)*Aea
		    end if
                
                  Smp=G1stare-G1starw+G2starn-G2stars
                  Su(I,J)=Smp+(Densit(I,J)-Densit_old(I,J))/deltat*Vol
		  Sp(I,J)=0.0
                  Resorm=Resorm+abs(Smp)
               End Do
            End Do
            Do I=2,nx-1
               Do J=2,ny-1
                  Ap(I,J)=Ae(I,J)+Aw(I,J)+An(I,J)+As(I,J) !-Sp(I,J) was added
               End Do
            End Do
            Do  Np=1,NswpP
            	Call Liner_Solver(2,2,Pp)
            End Do
		
C   ------------------ Endof Calc P ---------------------  
C   ------------------ Correction --------------------- 
C    -------------------- U
            Do  I=3,nx-1
                Do  J=2,ny-1
		Awe=Sns(J)
        U_fluid(I,J)=U_fluid(I,J)-Awe/Du(I,J)*(Pp(I,J)-Pp(I-1,J))
                End Do
            End Do
C    -------------------- V
            Do  I=2,nx-1
                Do  J=3,ny-1
		Asu=Sew(I)
        V_fluid(I,J)=V_fluid(I,J)-Asu/Dv(I,J)*(Pp(I,J)-Pp(I,J-1))
                End Do
            End Do
C    --------------------- Pressures (what is pressure refrence?) 
            Ppref=Pp(2,2)
            Do  I=2,nx-1
                Do J=2,ny-1
                   P(I,J)=P(I,J)+Urfp*(Pp(I,J)-Ppref)
                   Pp(I,J)=0.0
                End Do
            End Do
C   ------------------ Endof Correction --------------------- 	

	    Source1=max(Resoru,Resorv,Resorm)  

	If (current_time .GT. 0.01) then
	max_iter=100
	End If

	    If ( Kr .GT. max_iter) then  ! 500 is max itteration per time
c	    write (*,*)  Source1, Kr
		Exit
	    End If
	 End Do
         U_fluid_old=U_fluid
	 V_fluid_old=V_fluid
	 Densit_old=Densit	 

	
	End Subroutine get_fluid_Velocity
	
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
	!********************************************************!
	!--------------------------------------------------------!
	!	****	Saving the data for tecplot    	****	 !
	!--------------------------------------------------------!
	!********************************************************!
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
C       
	Subroutine Save_data
	Use List_of_Variables
	Write (1,*) 'Zone I=', nxp, '  J=', nyp, 't="', current_time, '"'   
	Do j=1,nyp
		Do i=1,nxp
		   index1 = i + nxp*(j-1)
		   Write (1,*) X(i+1), Y(j+1), alpha(index1), U(index1,1) !Update when Num_family is changed
     1            , U(index1,2), U(index1,3)
		End Do  
	End Do
c	Write (2,*) SUM(U(:,1))/(nx*ny),SUM(U(:,2))/(nx*ny),
c     2  SUM(U(:,4))/(nx*ny),SUM(U(:,5))/(nx*ny), 
c     3  SUM(U(:,7))/(nx*ny),SUM(U(:,8))/(nx*ny)	

	 Write (3,*) 'Zone I=', nx,' J=',ny,'t="',current_time,'"' 
        PP=P
	Pp(nx,ny)=P(nx-1,ny-1)
        Pp(nx,1)=P(nx-1,2)
	Pp(1,1)=P(2,2)
        Pp(1,ny)=P(2,ny-1)
	Do j=1,ny
	    Do i=1,nx
		 if (i .EQ. 1) then
		     Ug(i,j)=U_fluid(1,j)
		     Pp(1,J)=P(2,J)
		 elseif (i .EQ. nx) then
		     Ug(i,j)=U_fluid(nx,j)
		     Pp(nx,J)=P(nx-1,J)
		 else 
	             Ug(i,j)=(U_fluid(i,j)+U_fluid(i+1,j))/2.0
		 end if
		 if (j .EQ. 1) then
		     Vg(i,j)=V_fluid(i,1)
                     Pp(I,1)=P(I,2);
		 elseif (j .EQ. ny) then
		     Ug(i,j)=U_fluid(i,ny)
		     Pp(I,ny)=P(I,ny-1)
		 else
		     Vg(i,j)=(V_fluid(i,j)+V_fluid(i,j+1))/2.0
		 end if
		 
	       Write (3,*) X(i), Y(j), Densit(i,j)/ro_fluid, Ug(i,j),
     1         Vg(i,j), Pp(i,j)
	    End Do  
	 End Do

	End  Subroutine Save_data



