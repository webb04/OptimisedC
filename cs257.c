#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>

#ifdef __APPLE__
	#include <GLUT/glut.h>
	#include <xmmintrin.h>
#else
	#include <GL/glut.h>
	#include <immintrin.h>
#endif

#define FLOPS (20)

// Simulation parameters.
int N;
int steps;
float eps = 0.00125f;
float dmp = 0.995f;
float dt = 0.001f;
int vis = 0;

// Solution arrays.
float *x, *y, *z;
float *ax, *ay, *az;
float *vx, *vy, *vz;
float *m;
float *c;

// Timers.
double l0 = 0;
double l1 = 0;
double l2 = 0;
double l3 = 0;

// Current timestep.
int t = 0;

// Function prototypes.
double wtime();
void init();
void render();
void compute();
void update();
void check_exit(int t);
void verify();
void cleanup();

// Include coursework function.
#include "coursework.c"

/**
 * Return time in seconds in Epoch.
 */
double wtime() {
	struct rusage r;
	struct timeval t;
	gettimeofday( &t, (struct timezone *)0 );
	return t.tv_sec + t.tv_usec*1.0e-6;
}

/**
 * Set stuff up.
 */
void init() {

	// Allocate arrays, aligned to 16 bytes.
    // openMp threading here? #pragma omp parallel
	x = (float*) _mm_malloc(N * sizeof(float), 16);
	y = (float*) _mm_malloc(N * sizeof(float), 16);
	z = (float*) _mm_malloc(N * sizeof(float), 16);
	ax = (float*) _mm_malloc(N * sizeof(float), 16);
	ay = (float*) _mm_malloc(N * sizeof(float), 16);
	az = (float*) _mm_malloc(N * sizeof(float), 16);
    vx = (float*) _mm_malloc(N * sizeof(float), 16);
	vy = (float*) _mm_malloc(N * sizeof(float), 16);
    vz = (float*) _mm_malloc(N * sizeof(float), 16);
	m = (float*) _mm_malloc(N * sizeof(float), 16);
	
	// Color!
	c = (float*) malloc(3 * N * sizeof(float));

	// Initialise array contents.
	// srand(42) simply to make results repeatable. :)
	srand(42);
	float scale = 1.0f;
	float vscale = 1.0f;
	float mscale = 1.0f;
	for (int i = 0; i < N; i++) {
		x[i] = ((rand() / (float) RAND_MAX) * 2 - 1) * scale;
		y[i] = ((rand() / (float) RAND_MAX) * 2 - 1) * scale;
		z[i] = ((rand() / (float) RAND_MAX) * 2 - 1) * scale;
		ax[i] = 0.0f;
		ay[i] = 0.0f;
		az[i] = 0.0f;
		vx[i] = ((rand() / (float) RAND_MAX) * 2 - 1) * vscale;
		vy[i] = ((rand() / (float) RAND_MAX) * 2 - 1) * vscale;
		vz[i] = ((rand() / (float) RAND_MAX) * 2 - 1) * vscale;
		m[i] = (rand() / (float) RAND_MAX) * mscale;
		c[3*i+0] = (rand() / (float) RAND_MAX);
		c[3*i+1] = (rand() / (float) RAND_MAX);
		c[3*i+2] = (rand() / (float) RAND_MAX);				
	}
		    
    // Initialise visualisation stuff.
    if (vis) {
    	glutInitWindowPosition(100, 100);
    	glutInitWindowSize(500, 500);
    	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    	glutCreateWindow("CS257 Coursework - Visualisation Window");
    	
    	glClearColor(0, 0, 0, 0);
    	glEnable(GL_DEPTH_TEST);
        glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);    	
        glHint(GL_CLIP_VOLUME_CLIPPING_HINT_EXT, GL_FASTEST);	
    	
    	glutDisplayFunc(render); 
    	glutIdleFunc(update);
    }

}

/**
 * Render stuff.
 * (Disks in place of spheres accelerates rendering.)
 */
void render() {   
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);    
    glLoadIdentity();
    for (int i = 0; i < N; i++) {
        glPushMatrix();
        glTranslatef(x[i], y[i], z[i]);        
        glColor3f(c[3*i+0], c[3*i+1], c[3*i+2]);
        glBegin(GL_QUADS);
            glVertex3f(-m[i]*0.01f, -m[i]*0.01f, 0);
            glVertex3f(-m[i]*0.01f, m[i]*0.01f, 0);
            glVertex3f(m[i]*0.01f, m[i]*0.01f, 0);
            glVertex3f(m[i]*0.01f, -m[i]*0.01f, 0);
        glEnd();
        glPopMatrix();
    }

	glutSwapBuffers();

	return;

}

/**
 * Call the coursework compute() function to update the frame.
 * Sleep to maintain 60fps.
 */
#define FPS 60.0f
void update() {
	
	double t_start = wtime();
	compute();
	double t_end = wtime();
	if (t_end - t_start < (1.0f/FPS)) {
	    unsigned int usecs = (unsigned int) (((1.0f/FPS) - (t_end - t_start)) * 1E6);
            usleep(usecs);
	}
	
	if (vis) glutPostRedisplay();
	check_exit(t++);

}

/**
 * Check for exit condition.
 */
void check_exit(int t) {

	if (t >= steps) {
	
	    // Print results and stuff.
	    printf("\n");	
	    printf(" Loop 0 = %f seconds.\n", l0);
	    printf(" Loop 1 = %f seconds.\n", l1);
	    printf(" Loop 2 = %f seconds.\n", l2);
	    printf(" Loop 3 = %f seconds.\n", l3);
	    printf(" Total  = %f seconds.\n", l0 + l1 + l2 + l3);
	    printf("\n");
	    double flops = 20.0f * (double) N * (double) (N-1) * (double) steps;
	    printf(" GFLOP/s = %f\n", flops / 1000000000.0f / (l0 + l1 + l2 + l3));
	    
	    double bytes = 4.0f * (double) N * 10.0f * (double) steps;
	    printf(" GB/s = %f\n", bytes / 1000000000.0f / (l0 + l1 + l2 + l3));
	    printf("\n");	

	    // Verify solution.
	    verify();
	    printf("\n");	

	    // Tidy up.
	    cleanup();
            exit(0);
	    
	}
	
}

/**
 * Verify the solution.
 */
void verify() {
    // loop blocking?
    // loop interchange?
    float phi = 0.0f;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) { 
			float rx = x[j] - x[i];
			float ry = y[j] - y[i];
			float rz = z[j] - z[i];
			float r2 = rx*rx + ry*ry + rz*rz + eps;
			float r2inv = 1.0 / sqrt(r2);
			float r6inv = r2inv * r2inv * r2inv;
			phi += m[j] * r6inv;
		}
	}
	
	printf(" Answer = %f\n", phi);
	
}

/**
 * Free all allocated memory.
 */
void cleanup() {

	_mm_free(x);
	_mm_free(y);
	_mm_free(z);
	_mm_free(ax);
	_mm_free(ay);
	_mm_free(az);
	_mm_free(vx);
	_mm_free(vy);
	_mm_free(vz);
	_mm_free(m);
	free(c);

}

/**
 * Main function.
 */
int main(int argc, char* argv[]) {

	// Parse command line options.
	if (argc != 4) {
		printf("Usage: ./cs257 [# stars] [# timesteps] [visualisation (0 or 1)]\n");
		printf("Eg   : ./cs257 10000 100 1\n");
		exit(0);
	}
	N = atoi(argv[1]);
	steps = atoi(argv[2]);
	vis = atoi(argv[3]);

    // Print a banner.
    printf("\n");
	printf("CS257 Optimisation Coursework\n");
	printf("-----------------------------------------\n");
	printf(" Starting simulation of %d stars for %d timesteps.\n",  N, steps);

    // Initialise and run stuff.
    if (vis) glutInit(&argc, argv);   
	init();
	if (vis) {
	    glutMainLoop();
	} else {
        while (1) update();
	}

}
