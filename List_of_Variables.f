	Module List_of_Variables

	Integer i,j,k,index1,index2, t, Kr, N, Maxts, Np,Nu,Nv,t1
	Integer Flag, Flag1, Flag2, Flag3, Flag4, Flag5
 	Integer nx, ny, nxp, nyp
	Integer :: Maxit=500, NswpMax=20
	Integer :: num_mesh_in
	Integer :: j_start 
	Integer, parameter :: Num_family=8 ! Can be updated
	Integer, parameter :: NswpV=3, NswpU=3, NswpP=10

	Double Precision :: RBC_dens=1125.0, Gravity=-9.81
	Double Precision :: max_Pack=0.9
	Double Precision :: ro_fluid=1025.0, miu=3.0D-3, pi=3.1415
	Double Precision :: Viscos=3.0D-3
	Double Precision :: Vol_RBC=90.0D-18 !90 femto liter-18 -11 work for moving 1cm
	Double Precision :: Urfu=.8, Urfv=.8, Urfp=.9	!under relaxation
	Double Precision :: Great=10D30
	Double Precision :: t_final, Upin, Uin
	Double Precision L, H, sigma2
	Double Precision max_Den_p_in, deltat_start
	Double Precision current_time, max_deltat_x, max_deltat_y, u1, v1
	Double Precision Save_time, max_iter
	Double Precision pl, pr, pu, pd , F1, F2,CFL,epsilon1
	Double Precision d_pl, d_pr, d_pu, d_pd,Tmult
	Double Precision deltat, fastest_wave
	Double Precision, dimension(3) :: Fx, Fy, Ul, Ur, Ud, Uu
	Double Precision, dimension(3) :: Fl, Fr, Fu, Fd 
	Double Precision, dimension(3) :: fdiff, bdiff, d_ud_loc 
	Double Precision, dimension(Num_family) :: mass_fam, vol_fam 
	Double Precision, dimension(Num_family) ::  Drag_t, radius_fam  
	Double Precision, dimension(Num_family,Num_family):: B_Agg,D_Agg  

c	Double Precision, dimension(3) :: unknown=(/2.4D-4, 1.26D-6
c     1  ,1.6D-6/)
	Double Precision, dimension(8) :: unknown=(/2.4D-4, 1.26D-6,
     1   1.6D-6,5.6D-6,7.79D-6,5.79D-5,8.96D-4,2.2D-3/)

	Double Precision :: Asu, Awe, Aea, Ano,	Vol,Sormax, Source1
	Double Precision :: Ce, Cw, Cn, Cs, De, Dw, Dn, Ds
	Double Precision ::  Gue, Guw, Gvn, Gvs	
	Double Precision :: Gen, Spp, Smp, Cp
	Double Precision:: Resoru, Resorv, Resorm, Resor
	Double Precision:: G1stare, G1starw, G2starn, G2stars, ppref	

	Double Precision, Allocatable::  Dxep (:)
	Double Precision, Allocatable::  Dxpw (:)
	Double Precision, Allocatable::   Dynp(:)
	Double Precision, Allocatable::   Dyps(:)
	Double Precision, Allocatable::   Sns(:)
	Double Precision, Allocatable::   Sew(:)
	Double Precision, Allocatable::   Xu(:)
	Double Precision, Allocatable::   Yv(:)
	Double Precision, Allocatable::   Dxepu(:)
	Double Precision, Allocatable::   Dxpwu(:)
	Double Precision, Allocatable::   Sewu (:)
	Double Precision, Allocatable::   Dynpv(:)
	Double Precision, Allocatable::   Dypsv(:)
	Double Precision, Allocatable::   Snsv(:)
	Double Precision, Allocatable::   X(:)
	Double Precision, Allocatable::   Y(:)
	Double Precision, Allocatable::   Du(:,:)               
	Double Precision, Allocatable::   Dv(:,:)        
	Double Precision, Allocatable::   Ppp(:,:)        
	Double Precision, Allocatable::   Ap(:,:)        
	Double Precision, Allocatable::   An(:,:)        
	Double Precision, Allocatable::   As(:,:)        
	Double Precision, Allocatable::   Ae(:,:)        
	Double Precision, Allocatable::   Aw(:,:)        
	Double Precision, Allocatable::   Su(:,:)        
	Double Precision, Allocatable::   Sp(:,:)
	Double Precision, Allocatable::   Ug(:,:)   
	Double Precision, Allocatable::   Vg(:,:)
	Double Precision, Allocatable::   Up(:,:,:)
	Double Precision, Allocatable::   Den_p(:,:,:)   
	Double Precision, Allocatable::   Vp(:,:,:)
	Double Precision, Allocatable::   P (:,:)    
	Double Precision, Allocatable::   Pp(:,:)
	Double Precision, Allocatable::   Densit(:,:)
	Double Precision, Allocatable::   U_fluid(:,:)
	Double Precision, Allocatable::   V_fluid(:,:)
	Double Precision, Allocatable::   Densit_old(:,:)
	Double Precision, Allocatable::   U_fluid_old(:,:)
	Double Precision, Allocatable::   V_fluid_old(:,:)

	Double Precision, Allocatable::   Ugg(:,:)   
	Double Precision, Allocatable::   Vgg(:,:)
	Double Precision, Allocatable::   pressure(:)
	Double Precision, Allocatable::   alpha(:)
	Double Precision, Allocatable::   U (:,:) 
	Double Precision, Allocatable::   dUdt(:,:)
	Double Precision, Allocatable::   dUdx(:,:)
	Double Precision, Allocatable::   dUdy(:,:)
	Double Precision, Allocatable::   U_old(:,:)
	Double Precision, Allocatable::   V_old(:,:)
	Double Precision, Allocatable::   Source(:,:)
	Double Precision, Allocatable::   dAggdt(:,:)
	Double Precision, Allocatable::   di_pressure_di_alpha(:)
	
	End Module List_of_Variables
