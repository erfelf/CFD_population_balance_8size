	Subroutine Liner_Solver(Istart,Jstart,Phi)
	Use List_Of_Variables

	  Integer Jstm1,Jstart,Istart,Jj
	  Double Precision  AA(ny),BB(ny),CC(ny),DD(ny),Term
	  Double Precision,Intent(inOut), Dimension(nx,ny):: Phi
	  
	Jstm1=Jstart-1
	AA(Jstm1)=0.0

!-----Commence W-E Sweep
	Do i=Istart,nx-1
           CC(Jstm1)=Phi(i,Jstm1)

!-----Commence S-N Traverse
	   Do J=Jstart,ny-1
!-----Assemble Toma Coefficients
	      AA(J)=An(i,J)
	      BB(J)=As(i,J)
	      CC(J)=Ae(i,J)*Phi(i+1,J)+Aw(i,J)*Phi(i-1,J)
     1        +Su(i,J)
	      DD(J)=Ap(i,J)
!-----Calculate Coefficients Of Recorrence Formula
              Term=1.0/(DD(J)-BB(J)*AA(J-1))
	      AA(J)=AA(J)*Term
              CC(J)=(CC(J)+BB(J)*CC(J-1))*Term
            End Do

!-----Obtain New Phi's
            Do Jj=Jstart,ny-1
	       J=ny+Jstm1-Jj
               Phi(i,J)=AA(J)*Phi(i,J+1)+CC(J)
            End Do
	End Do

	End Subroutine Liner_Solver

