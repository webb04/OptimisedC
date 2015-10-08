/**
 * The function to optimise as part of the coursework.
 *
 * l0, l1, l2 and l3 record the amount of time spent in each loop
 * and should not be optimised out. :)
 */
void compute() {

	double t0, t1;

	// Loop 0.
	t0 = wtime();
  #pragma omp parallel for num_threads(4)
	for (int i = 0; i < N; i+=4) {
		__m128 m0 = _mm_set1_ps(0.0f);
		_mm_store_ps(ax+i, m0);
		_mm_store_ps(ay+i, m0);
		_mm_store_ps(az+i, m0);
		// check to see if unrolling is just quicker for loop 0
	}

	t1 = wtime();
	l0 += (t1 - t0);http:

	// Loop 1.
	t0 = wtime();
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < N; i++) {
		float valuesx[4];
		__m128 xi = _mm_set1_ps(x[i]); // vector xi is x[i]
		__m128 yi = _mm_set1_ps(y[i]); // vector yi is y[i]
		__m128 zi = _mm_set1_ps(z[i]);


		for (int j = 0; j < N; j+=4) {

			//storing is exspensive
			__m128 vec_rx = _mm_load_ps(x+j);
			vec_rx = _mm_sub_ps(vec_rx,xi); // first 4 j-i

			__m128 vec_ry = _mm_load_ps(y+j);
			vec_ry = _mm_sub_ps(vec_ry,yi); // first 4 j-i

			__m128 vec_rz = _mm_load_ps(z+j);
			vec_rz = _mm_sub_ps(vec_rz,zi); // first 4 j-i

		  __m128 vec_eps = _mm_set1_ps(eps);
			__m128 vec_r2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(vec_rx,vec_rx),_mm_mul_ps(vec_ry,vec_ry)), _mm_add_ps(_mm_mul_ps(vec_rz,vec_rz), vec_eps));

			__m128 vec_r2inv = _mm_rsqrt_ps(vec_r2);

		  __m128 vec_r6inv =_mm_mul_ps(vec_r2inv,_mm_mul_ps(vec_r2inv, vec_r2inv));

			__m128 vec_mj = _mm_load_ps(m+j);
			__m128 vec_s = _mm_mul_ps(vec_mj, vec_r6inv);

			__m128 vec_10 = _mm_mul_ps(vec_s, vec_rx);
			_mm_store_ps(valuesx, vec_10);
			ax[i] += valuesx[0];
			ax[i] += valuesx[1];
			ax[i] += valuesx[2];
			ax[i] += valuesx[3];

			__m128 vec_16 = _mm_mul_ps(vec_s, vec_ry);
			_mm_store_ps(valuesx, vec_16);
			ay[i] += valuesx[0];
			ay[i] += valuesx[1];
		        ay[i] += valuesx[2];
			ay[i] += valuesx[3];

			__m128 vec_22 = _mm_mul_ps(vec_s, vec_rz);
			_mm_store_ps(valuesx, vec_22);
			az[i] += valuesx[0];
			az[i] += valuesx[1];
		  az[i] += valuesx[2];
			az[i] += valuesx[3];

		}

	}

	t1 = wtime();
	l1 += (t1 - t0);

	//Loop 2.
	t0 = wtime();
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < N; i++) {

		vx[i] += dmp * (dt * ax[i]);
		x[i] += dt * vx[i];
		if (x[i] >= 1.0f || x[i] <= -1.0f) vx[i] *= -1.0f;

		vy[i] += dmp * (dt * ay[i]);
		y[i] += dt * vy[i];
		if (y[i] >= 1.0f || y[i] <= -1.0f) vy[i] *= -1.0f;

		vz[i] += dmp * (dt * az[i]);
		z[i] += dt * vz[i];
		if (z[i] >= 1.0f || z[i] <= -1.0f) vz[i] *= -1.0f;
	}
	
	t1 = wtime();
	l2 += (t1 - t0);

	// Loop 3.
	t0 = wtime();

	t1 = wtime();
	l3 += (t1 - t0);

}
