////////////////////////////////////////////////////////////////////
//  SVD_Test_ND.C
//  Dec. 2018 - J. Hill
////////////////////////////////////////////////////////////////////
//  Performs a stochastic singular value decomposition on arbitrarily-rotated N-dimensional square matrices
//  Returns eigenvalues and eigenvectors as well as visualizations of the original, rotated, and reduced matrices

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <chrono>
#include <vector>
#include <string>

//generates a visualization of the matrix A formatted either for viewing in terminal or LaTeX
void PrintMatrix(gsl_matrix* A, std::string format)
{
  if (format == "standard")
    {
      for (int i = 0; i < A->size1; i++)
	{
	  std::cout<<std::endl;
	  std::cout.width(20);
	  std::cout.fill(' ');
	  std::cout<<std::right<<"[";
	  
	  for (int j = 0; j < A->size2; j++)
	    {
	      std::cout.width(20);
	      std::cout.fill(' ');
	      std::cout<<std::right<<gsl_matrix_get(A,i,j);
	      
	    }
	  
	  std::cout.width(10);
	  std::cout.fill(' ');
	  std::cout<<std::right<<"]"<<std::endl;
	}
    }

  else if (format == "latex")
    {
      std::cout.precision(3);

      for (int i = 0; i < A->size1; i++)
	{		 
	  for (int j = 0; j < A->size2; j++)
	    {
	      std::cout<<std::right<<gsl_matrix_get(A,i,j);

	      if (j < A->size2-1)
		{
		  std::cout<<" & ";	      
		}
	    }
	  
	  std::cout<<std::right<<"\\\\"<<std::endl;
	}

      std::cout.precision(5);
    }

  return;
}

//rotates A about the plane formed by axes i and j by angle theta
void Rotate_ijPlane(gsl_matrix* A, int i, int j, double theta)
{
  if (A->size1 != A->size2)
    {
      std::cout<<"Basis rotation is only possible for square matrices"<<std::endl;
      return;
    }

  double val = 0;
  int stat;
  int N = A->size1;
  gsl_matrix* A_prime = gsl_matrix_calloc(N,N);
  double s = sin(theta);
  double c = cos(theta);

  stat = gsl_matrix_memcpy(A_prime,A);

  if (i >= N || j >= N)
    {
      std::cout<<"Specified rotation indices out of range"<<std::endl;
      return;
    }

  for (int ii = 0; ii < N; ii++)
    {
      if (ii == i)
	{
	  gsl_matrix_set(A_prime,i,i,gsl_pow_2(c)*gsl_matrix_get(A,i,i) - c*s*(gsl_matrix_get(A,i,j) + gsl_matrix_get(A,j,i)) + gsl_pow_2(s)*gsl_matrix_get(A,j,j));
	  gsl_matrix_set(A_prime,i,j,c*s*gsl_matrix_get(A,i,i) - gsl_pow_2(s)*gsl_matrix_get(A,j,i) + s*c*gsl_matrix_get(A,i,j) - s*c*gsl_matrix_get(A,j,j));
	}

      else if (ii == j)
	{
	  gsl_matrix_set(A_prime,j,i,s*c*gsl_matrix_get(A,i,i) + gsl_pow_2(c)*gsl_matrix_get(A,j,i) - gsl_pow_2(s)*gsl_matrix_get(A,i,j) - s*c*gsl_matrix_get(A,j,j));
	  gsl_matrix_set(A_prime,j,j,gsl_pow_2(s)*gsl_matrix_get(A,i,i) + s*c*(gsl_matrix_get(A,j,i) + gsl_matrix_get(A,i,j)) + gsl_pow_2(c)*gsl_matrix_get(A,j,j));
	}

      else
	{
	  gsl_matrix_set(A_prime,ii,i,c*gsl_matrix_get(A,ii,i) - s*gsl_matrix_get(A,ii,j));
	  gsl_matrix_set(A_prime,i,ii,c*gsl_matrix_get(A,i,ii) - s*gsl_matrix_get(A,j,ii));

	  gsl_matrix_set(A_prime,ii,j,s*gsl_matrix_get(A,ii,i) + c*gsl_matrix_get(A,ii,j));
	  gsl_matrix_set(A_prime,j,ii,s*gsl_matrix_get(A,i,ii) + c*gsl_matrix_get(A,j,ii));
	}
    }


    gsl_matrix_memcpy(A,A_prime);
    gsl_matrix_free(A_prime);
  return;
}

//orthonormalizes the column vectors of Y via the Graham-Schmidt process
gsl_matrix* Orthonormalize(gsl_matrix* Y)
{
  int M = Y->size1;
  int N = Y->size2;
  gsl_matrix* Q = gsl_matrix_calloc(M,N);
  gsl_vector* p = gsl_vector_calloc(M);
  gsl_vector* q = gsl_vector_calloc(M);

  double magnitude;
  double mag_p1 = 0;
  double mag_p2 = 0;
  double p_dot_q = 0;

  for (int i = 0; i < N; i++)
    {
      magnitude = 0;
      p_dot_q = 0;

      for (int j = 0; j < M; j++)
	{
          gsl_vector_set(p,j,gsl_matrix_get(Y,j,i));
          gsl_vector_set(q,j,0);
        }

      for (int k = 0; k < i; k++)
        {
          p_dot_q = 0;

	  for (int m = 0; m < M; m++)
            {
              p_dot_q += gsl_vector_get(p,m)*gsl_matrix_get(Q,m,k);
            }

          for (int m = 0; m < M; m++)
            {
              gsl_vector_set(p,m,gsl_vector_get(p,m) - gsl_matrix_get(Q,m,k)*p_dot_q);                                                                                  
            }
        }//end of k loop                                                                                                                                                                                

      for (int j = 0; j < M; j++)
        {
          magnitude += gsl_pow_2(gsl_vector_get(p,j));
	}

      magnitude = pow(magnitude,0.5);

      if (magnitude == 0)
        {
          magnitude = 1;
        }

      for (int j = 0; j < M; j++)
        {
	  gsl_matrix_set(Q,j,i,gsl_vector_get(p,j)/magnitude);
        }

    }//end of i loop                                                                                                                                                                                    

  return Q;

}//end of Orthonormalize                                                                                                                                                                                
//finds the reduced rank of tensor H and stores it as m
gsl_matrix*  ReduceHilbertSpace(gsl_matrix* H, int& m)
{
  bool complete = false;
  int N = H->size1;

  m = 1;
  int num_samples = 1;
  double val,tr,prev,tol;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  double sigma = 1;
  clock_t seed = clock();
  gsl_rng_set(r,seed);

  // prev is a metric of the previous rank reduction's completeness //initialized as a large number to guarantee first iteration is completed
  prev = 1e4;
  tol = 1e-1;

 checkpoint:
  gsl_matrix* W = gsl_matrix_alloc(N,m);
  val = 0;

  //generate the random gaussian matrix W
  for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < m; j++)
	{
	  for (int k = 0; k < N; k++)
	    {
	      for (int l = 0; l < N; l++)
		{
		  for (int n = 0; n < num_samples; n++)
		    {
		      val += gsl_matrix_get(H,i,k)*gsl_matrix_get(H,l,k)*
			     gsl_matrix_get(H,l,k)*gsl_ran_gaussian(r,sigma);
		    }
		}
	    }//end of k loop

	  gsl_matrix_set(W,i,j,val);
	  val = 0;
	}
    }

  gsl_matrix* Q_prime = Orthonormalize(W);
  tr = 0;
  val = 0;

  //calculate the trace of Q_prime^2 //should approach unity as m->N
  for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < m; j++)
	{
	  val += gsl_pow_2(gsl_matrix_get(Q_prime,i,j));
	}

      tr += val;
      val = 0;
    }

  tr /= m;
  
  //check to see if the diagonal components of Q_prime^2 are normalized to within fraction tol
  //if not, incrment the number of components m and restart
  if (fabs(1-tr) > tol && fabs(1-tr) < prev && m < N-1)
    {
      prev = fabs(1-tr);
      gsl_matrix_free(W);
      gsl_matrix_free(Q_prime);
      m++;
      
      goto checkpoint;
    }
    
  gsl_matrix_free(W);
  gsl_rng_free(r);

  return Q_prime;
}

//builds the rank-reduced approximation of H from the singular vectors U and the sign-corrected singular values (i.e. eigenvalues) s
gsl_matrix* ReconstructH(gsl_matrix* U, gsl_vector* s)
{
  int M = s->size;
  int N = U->size1;
  gsl_matrix* H = gsl_matrix_calloc(N,N);
  double val = 0;

  for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < N; j++)
	{
	  for (int k = 0; k < M; k++)
	    {
	      val += gsl_matrix_get(U,i,k)*gsl_matrix_get(U,j,k)*gsl_vector_get(s,k);
	    }

	  gsl_matrix_set(H,i,j,val);
	  val = 0;
	}
    }

  return H;
}

int main()
{
  int N = 7;//dimension of the test space
  int n = 4;//number of nonzero eigenvalues
  std::vector<double>E;//energy eiegenvalues of the test space
  E.resize(N,0);


  //initialize the eigenvalues //
  for (int i = 0; i < n; i++)
    {
      E[i] = 1;
    }

  //unique eigenvalues can be specified via E[i]
  E[0] = -50;
  E[3] = 50;
  ///////////////////////////////

  double val = 0;

  //current rotation scheme is to rotate the (i-1, i) plane by theta_ij*ratio^i for i in the range [1,N-1]
  double theta_ij = 45;
  double ratio = 0.8;
  theta_ij *= M_PI/180;

  int status,m;

  gsl_matrix* M_prime = gsl_matrix_alloc(N,N);
  gsl_matrix* M = gsl_matrix_alloc(N,N);

  //set the matrix elements according to the values in E
  for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < N; j++)
	{
	  if (i != j)
	    {
	      gsl_matrix_set(M,i,j,0);
	    }

	  else
	    {
	      gsl_matrix_set(M,i,j,E[i]);
	    }
	}
    }

  std::cout<<"------M------"<<std::endl;
  PrintMatrix(M,"standard");

  for (int i = 1; i < N; i++)
    {
      Rotate_ijPlane(M,i-1,i,theta_ij);
      theta_ij *= ratio;
    }

  /////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
  std::cout<<"Rotating M"<<std::endl;
  std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////


  std::cout<<"------M------"<<std::endl;
  PrintMatrix(M,"latex");

  gsl_matrix* Q = ReduceHilbertSpace(M,m);
  std::cout<<"m = "<<m<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
  std::cout<<"generating random orthonormal basis"<<std::endl;
  std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////

  gsl_matrix* QHQ = gsl_matrix_alloc(m,m);
  gsl_matrix* H;
  gsl_matrix* U = gsl_matrix_alloc(N,m);
  gsl_matrix* V = gsl_matrix_alloc(m,m);
  gsl_vector* s = gsl_vector_alloc(m);
  gsl_vector* work = gsl_vector_alloc(m);

  std::cout<<"------Q------"<<std::endl;
  PrintMatrix(Q,"standard");

  val = 0;

  //compute rank-reduced approximation of M
  for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < m; j++)
	{
	  for (int k = 0; k < N; k++)
	    {
	      for (int l = 0; l < N; l++)
		{
		  val += gsl_matrix_get(Q,k,i)*gsl_matrix_get(M,k,l)*gsl_matrix_get(Q,l,j);
		}
	    }

	  gsl_matrix_set(QHQ,i,j,val);
	  val = 0;
	}//end of j loop
    }

  std::cout<<"------Q†HQ-----"<<std::endl;
  PrintMatrix(QHQ,"latex");

  /////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
  std::cout<<"Calculating SVD"<<std::endl;
  std::cout<<"/////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
  /////////////////////////////////////////////////////////////////////////////////////////

  //SV-decompose the rank-reduced matrix QHQ
  //note that gsl stores left singular vectors U in the input QHQ
  status = gsl_linalg_SV_decomp(QHQ,V,s,work);

  std::cout<<"------U-----"<<std::endl;
  PrintMatrix(QHQ,"standard");

  std::vector<double> signs(m);
  val = 0;

  std::cout<<"eigenvalues: ";

  //determine the sign corrections of the singular values via the relationship between left and right singular vectors U and V
  for (int i = 0; i < m; i++)
    {
      signs[i] = gsl_matrix_get(QHQ,i,i)/gsl_matrix_get(V,i,i);
      std::cout<<signs[i]*gsl_vector_get(s,i);
      gsl_vector_set(s,i,signs[i]*gsl_vector_get(s,i));

      if (i < m-1)
	{
	  std::cout<<" || ";
	}
    }

  //generate full-rank approximate eigenvectors of H and store them in U
  for (int i = 0; i < N; i++)
    {
      for (int j = 0; j < m; j++)
	{
	  for (int k = 0; k < m; k++)
	    {
	      val += gsl_matrix_get(Q,i,k)*gsl_matrix_get(QHQ,k,j);
	    }

	  gsl_matrix_set(U,i,j,val);
	  val = 0;
	}
    }

  //generate full-rank approximation of H from U and s
  H = ReconstructH(U,s);

  std::cout<<std::endl;

  std::cout<<"------H-----"<<std::endl;
  PrintMatrix(H,"standard");


  gsl_matrix_free(H);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_vector_free(s);
  gsl_vector_free(work);
  gsl_matrix_free(QHQ);
  gsl_matrix_free(M);
  gsl_matrix_free(Q);

  return 0;
}
