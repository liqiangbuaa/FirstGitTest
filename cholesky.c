
//调用方式为     DSPF_sp_cholesky_cmplx(N, Rxx, L);
//        		 DSPF_sp_cholesky_solver_cmplx(N, L, y, b, x);
//	其中 Rxx 为输入12*12自相关矩阵，x并不是输出的自相关矩阵的逆矩阵，而是
//  输出自相关矩阵的第一列
//
//float	Rxx[2*N*N], L[2*N*N], y[2*N], x[2*N];

int DSPF_sp_cholesky_cmplx(const int Nrows, float *A, float *L)
{
  short i,j,k,Ncols;
  float xreal,ximag,yreal,yimag;
  float sum_real,sum_imag;
  float inv_mag_sq;
  /* generate lower diagonal matrix L */
  Ncols=2*Nrows;
  for (j=0;j<Nrows;j++) {
	sum_real=0.0;

	for (k=0;k<=2*j-1;k+=2) {
	  //xreal=L[j*Ncols+k  ];
	  //ximag=L[j*Ncols+k+1];
	  //sum_real+=xreal*xreal+ximag*ximag;
		sum_real+=L[j*Ncols+k]*L[j*Ncols+k]+L[j*Ncols+k+1]*L[j*Ncols+k+1];
	}//end

	xreal=A[j*Ncols+2*j  ]-sum_real;
	ximag=A[j*Ncols+2*j+1];

	  float x_angle,z_angle;
	  /* magnitude */
	  //sum_real=sqrt(xreal*xreal+ximag*ximag);
	  //sum_imag=sqrt(sum_real);
	  sum_real=sqrtsp_i(xreal*xreal+ximag*ximag);
	  sum_imag=sqrtsp_i(sum_real);
	  /* angle */
	  x_angle=atan2(ximag,xreal);
	  z_angle=x_angle/2;
	  /* results */
	  L[j*Ncols+2*j  ]=cos(z_angle)*sum_imag;
	  L[j*Ncols+2*j+1]=sin(z_angle)*sum_imag;

    for (i=j+1;i<Nrows;i++) {
      sum_real=0.0;
      sum_imag=0.0;
      for (k=0;k<=2*j-1;k+=2) {
    	//xreal= L[i*Ncols+k  ];
        //ximag= L[i*Ncols+k+1];
        //yreal= L[j*Ncols+k  ];
        //yimag=-L[j*Ncols+k+1];
    	//sum_real+=xreal*yreal-ximag*yimag;
    	//sum_imag+=ximag*yreal+yimag*xreal;
    	  sum_real+=L[i*Ncols+k]*L[j*Ncols+k]+L[i*Ncols+k+1]*L[j*Ncols+k+1];
    	  sum_imag+=L[i*Ncols+k+1]*L[j*Ncols+k]-L[j*Ncols+k+1]*L[i*Ncols+k];
      }
      xreal=A[i*Ncols+2*j  ]-sum_real;
      ximag=A[i*Ncols+2*j+1]-sum_imag;
      yreal=L[j*Ncols+2*j  ];
      yimag=L[j*Ncols+2*j+1];
      inv_mag_sq=1/(yreal*yreal+yimag*yimag);
      L[i*Ncols+2*j  ]=(xreal*yreal+ximag*yimag)*inv_mag_sq;
      L[i*Ncols+2*j+1]=(ximag*yreal-xreal*yimag)*inv_mag_sq;
 	}
  }//end
  return 0;
}

// LL*inv(R)b = b; force  y = L*inv(R)b;then the above can be changed to Ly = b;
//from the above function "DSPF_sp_cholesky_cmplx" we can get L, then y can be otbained by the equation Ly=b 
//since L has been getted, L* can be obtained obviously. then according to the equation y = L*inv(R)b, we can easily
//get the inv(R)b, which is the numerator of Wopt, the denominator of Wopt is the first number of inv(R)b.
//OK,this function "DSPF_sp_cholesky_solver_cmplx" gives the result x which represents the result inv(R)b.


int DSPF_sp_cholesky_solver_cmplx(const int Nrows,float *L,float *y,float *b,float *x)
{
	short i,k,Ncols;
	float sum_real,sum_imag,xreal,ximag,yreal,yimag;
	float inv_mag_sq;

  /* solve L*y=b for y using forward substitution */
	Ncols=2*Nrows;
	inv_mag_sq=1/(L[0]*L[0]+L[1]*L[1]);
	y[0]=(b[0]*L[0]+b[1]*L[1])*inv_mag_sq;
	y[1]=(b[1]*L[0]-b[0]*L[1])*inv_mag_sq;

	for (i=1;i<Nrows;i++) {
		sum_real=0;
		sum_imag=0;
		for (k=0;k<=i-1;k++) {
			//xreal=L[i*Ncols+2*k  ];
			//ximag=L[i*Ncols+2*k+1];
			//yreal=y[2*k  ];
			//yimag=y[2*k+1];
			//sum_real+=xreal*yreal-ximag*yimag;
			//sum_imag+=xreal*yimag+ximag*yreal;
			sum_real+=L[i*Ncols+2*k]*y[2*k]-L[i*Ncols+2*k+1]*y[2*k+1];
			sum_imag+=L[i*Ncols+2*k]*y[2*k+1]+L[i*Ncols+2*k+1]*y[2*k];
		}
		xreal=b[2*i  ]-sum_real;
		ximag=b[2*i+1]-sum_imag;
		yreal=L[i*Ncols+2*i  ];
		yimag=L[i*Ncols+2*i+1];
		inv_mag_sq=1/(yreal*yreal+yimag*yimag);
		y[2*i  ]=(xreal*yreal+ximag*yimag)*inv_mag_sq;
		y[2*i+1]=(ximag*yreal-xreal*yimag)*inv_mag_sq;
	}

  /* solve U*x=y for x using backward substitution */
	inv_mag_sq=1/(L[286]*L[286]+L[287]*L[287]);
	x[22]=(y[22]*L[286]-y[23]*L[287])*inv_mag_sq;
	x[23]=(y[23]*L[286]+y[22]*L[287])*inv_mag_sq;

	for (i=Nrows-2;i>=0;i--) {
  		sum_real=0;
  		sum_imag=0;
  		for (k=Nrows-1;k>=i+1;k--) {
  			//xreal= L[k*Ncols+2*i  ];
  			//ximag=-L[k*Ncols+2*i+1];
  			//yreal=x[2*k  ];
  			//yimag=x[2*k+1];
  			//sum_real+=xreal*yreal-ximag*yimag;
  			//sum_imag+=xreal*yimag+ximag*yreal;
  			sum_real+=L[k*Ncols+2*i]*x[2*k]+L[k*Ncols+2*i+1]*x[2*k+1];
  			sum_imag+=L[k*Ncols+2*i]*x[2*k+1]-L[k*Ncols+2*i+1]*x[2*k];
  		}
  		xreal=y[2*i  ]-sum_real;
  		ximag=y[2*i+1]-sum_imag;
  		yreal= L[i*Ncols+2*i  ];
  		yimag=-L[i*Ncols+2*i+1];
  		inv_mag_sq=1/(yreal*yreal+yimag*yimag);
  		x[2*i  ]=(xreal*yreal+ximag*yimag)*inv_mag_sq;
  		x[2*i+1]=(ximag*yreal-xreal*yimag)*inv_mag_sq;
	}
	return 0;
}

//#pragma CODE_SECTION(sqrtsp_i, ".text:optci");

static inline float sqrtsp_i (float a)
{
  const float  Half  = 0.5f;
  const float  OneP5 = 1.5f;
  float  x, y;
  int i;

  x = _rsqrsp(a);

  #pragma UNROLL(1)
  for(i=0; i< 2; i++){
    x = x * (OneP5 - (a * x * x * Half));
  }
  y = a * x;

  if (a <= 0.0f) {
    y = 0.0f;
  }
  if (a > FLT_MAX) {
    y = FLT_MAX;
  }

  return (y);
}




//调用方式为     DSPF_sp_cholesky_cmplx(N, Rxx, L);
//        		 DSPF_sp_cholesky_solver_cmplx(N, L, y, b, x);
//	其中 Rxx 为输入12*12自相关矩阵，x并不是输出的自相关矩阵的逆矩阵，而是
//  输出自相关矩阵的第一列
//
//float	Rxx[2*N*N], L[2*N*N], y[2*N], x[2*N];

int DSPF_sp_cholesky_cmplx(const int Nrows, float *A, float *L)
{
  short i,j,k,Ncols;
  float xreal,ximag,yreal,yimag;
  float sum_real,sum_imag;
  float inv_mag_sq;
  /* generate lower diagonal matrix L */
  Ncols=2*Nrows;
  for (j=0;j<Nrows;j++) 
  {
	sum_real=0.0;
	for (k=0;k<=2*j-1;k+=2) 
	{
		sum_real+=L[j*Ncols+k]*L[j*Ncols+k]+L[j*Ncols+k+1]*L[j*Ncols+k+1];
	}

	xreal=A[j*Ncols+2*j  ]-sum_real;
	ximag=A[j*Ncols+2*j+1];

	  float x_angle,z_angle;
	  /* magnitude */
	  sum_real=sqrtsp_i(xreal*xreal+ximag*ximag);
	  sum_imag=sqrtsp_i(sum_real);
	  /* angle */
	  x_angle=atan2(ximag,xreal);
	  z_angle=x_angle/2;
	  /* results */
	  L[j*Ncols+2*j  ]=cos(z_angle)*sum_imag;
	  L[j*Ncols+2*j+1]=sin(z_angle)*sum_imag;

    for (i=j+1;i<Nrows;i++) 
	{
      sum_real=0.0;
      sum_imag=0.0;
      for (k=0;k<=2*j-1;k+=2) 
	  {
    	  sum_real+=L[i*Ncols+k]*L[j*Ncols+k]+L[i*Ncols+k+1]*L[j*Ncols+k+1];
    	  sum_imag+=L[i*Ncols+k+1]*L[j*Ncols+k]-L[j*Ncols+k+1]*L[i*Ncols+k];
      }
      xreal=A[i*Ncols+2*j  ]-sum_real;
      ximag=A[i*Ncols+2*j+1]-sum_imag;
      yreal=L[j*Ncols+2*j  ];
      yimag=L[j*Ncols+2*j+1];
      inv_mag_sq=1/(yreal*yreal+yimag*yimag);
      L[i*Ncols+2*j  ]=(xreal*yreal+ximag*yimag)*inv_mag_sq;
      L[i*Ncols+2*j+1]=(ximag*yreal-xreal*yimag)*inv_mag_sq;
 	}
  }//end
  return 0;
}

// LL*inv(R)b = b; force  y = L*inv(R)b;then the above can be changed to Ly = b;
//from the above function "DSPF_sp_cholesky_cmplx" we can get L, then y can be otbained by the equation Ly=b 
//since L has been getted, L* can be obtained obviously. then according to the equation y = L*inv(R)b, we can easily
//get the inv(R)b, which is the numerator of Wopt, the denominator of Wopt is the first number of inv(R)b.
//OK,this function "DSPF_sp_cholesky_solver_cmplx" gives the result x which represents the result inv(R)b.


int DSPF_sp_cholesky_solver_cmplx(const int Nrows,float *L,float *y,float *b,float *x)
{
	short i,k,Ncols;
	float sum_real,sum_imag,xreal,ximag,yreal,yimag;
	float inv_mag_sq;

  /* solve L*y=b for y using forward substitution */
	Ncols=2*Nrows;
	inv_mag_sq=1/(L[0]*L[0]+L[1]*L[1]);
	y[0]=(b[0]*L[0]+b[1]*L[1])*inv_mag_sq;
	y[1]=(b[1]*L[0]-b[0]*L[1])*inv_mag_sq;

	for (i=1;i<Nrows;i++) {
		sum_real=0;
		sum_imag=0;
		for (k=0;k<=i-1;k++) 
		{
			sum_real+=L[i*Ncols+2*k]*y[2*k]-L[i*Ncols+2*k+1]*y[2*k+1];
			sum_imag+=L[i*Ncols+2*k]*y[2*k+1]+L[i*Ncols+2*k+1]*y[2*k];
		}
		xreal=b[2*i  ]-sum_real;
		ximag=b[2*i+1]-sum_imag;
		yreal=L[i*Ncols+2*i  ];
		yimag=L[i*Ncols+2*i+1];
		inv_mag_sq=1/(yreal*yreal+yimag*yimag);
		y[2*i  ]=(xreal*yreal+ximag*yimag)*inv_mag_sq;
		y[2*i+1]=(ximag*yreal-xreal*yimag)*inv_mag_sq;
	}

  /* solve U*x=y for x using backward substitution */
	inv_mag_sq=1/(L[286]*L[286]+L[287]*L[287]);
	x[22]=(y[22]*L[286]-y[23]*L[287])*inv_mag_sq;
	x[23]=(y[23]*L[286]+y[22]*L[287])*inv_mag_sq;

	for (i=Nrows-2;i>=0;i--) {
  		sum_real=0;
  		sum_imag=0;
  		for (k=Nrows-1;k>=i+1;k--) 
		{
  			sum_real+=L[k*Ncols+2*i]*x[2*k]+L[k*Ncols+2*i+1]*x[2*k+1];
  			sum_imag+=L[k*Ncols+2*i]*x[2*k+1]-L[k*Ncols+2*i+1]*x[2*k];
  		}
  		xreal=y[2*i  ]-sum_real;
  		ximag=y[2*i+1]-sum_imag;
  		yreal= L[i*Ncols+2*i  ];
  		yimag=-L[i*Ncols+2*i+1];
  		inv_mag_sq=1/(yreal*yreal+yimag*yimag);
  		x[2*i  ]=(xreal*yreal+ximag*yimag)*inv_mag_sq;
  		x[2*i+1]=(ximag*yreal-xreal*yimag)*inv_mag_sq;
	}
	return 0;
}

//#pragma CODE_SECTION(sqrtsp_i, ".text:optci");

static inline float sqrtsp_i (float a)
{
  const float  Half  = 0.5f;
  const float  OneP5 = 1.5f;
  float  x, y;
  int i;

  x = _rsqrsp(a);

  #pragma UNROLL(1)
  for(i=0; i< 2; i++){
    x = x * (OneP5 - (a * x * x * Half));
  }
  y = a * x;

  if (a <= 0.0f) {
    y = 0.0f;
  }
  if (a > FLT_MAX) {
    y = FLT_MAX;
  }

  return (y);
}
