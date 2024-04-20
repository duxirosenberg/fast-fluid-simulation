#include "LBM.h"


#ifdef __x86_64__
#include "tsc_x86.h"
#endif






#ifdef LBM_STRUCT 
double perform_Measurement(int N, struct LBM* S, int  time)
#else
double perform_Measurement(int N, int nX, int nY, int nZ, int direction_size, 
                        int time, double tau, double gamma_dot, double c_s, 
                        int boundary_condition,
                        double* density_field,
                        double* velocity_field,
                        double* previous_particle_distributions,
                        double* particle_distributions,
                        const int* directions,
                        const double* weights,
                        int* reverse_indexes)
#endif
{
    clock_t start, end;


    // CALIBRATE; warm up CPU dont know ...

    int n = fmin(50, N/10);
    int num_runs = N-n;
    int i = time;



    for(; i <= n; ++i) {       

        #ifdef LBM_STRUCT 
		perform_timestep(S, i);
		#else
		perform_timestep(nX, nY, nZ, direction_size, i, tau, gamma_dot, c_s, boundary_condition, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);
        #endif		
	}

    start = clock();
    for(; i <= N; i++) {
        

        #ifdef LBM_STRUCT 
		perform_timestep(S, i);
		#else
		perform_timestep(nX, nY, nZ, direction_size, i, tau, gamma_dot, c_s, boundary_condition, density_field, velocity_field, previous_particle_distributions, particle_distributions, directions, weights, reverse_indexes);
        #endif
		
	}
    end = clock();


    return (double)(end-start)/num_runs/CLOCKS_PER_SEC*1000;

}












// double c_clock(double A[], double B[], int n) {
//     int i, num_runs;
//     double cycles;
//     clock_t start, end;

//     num_runs = NUM_RUNS;
// #ifdef CALIBRATE
//     while(num_runs < (1 << 14)) {
//         start = clock();
//         for (i = 0; i < num_runs; ++i) {
//             compute(A, B, n);
//         }
//         end = clock();

//         cycles = (double)(end-start);

//         // Same as in c_clock: CYCLES_REQUIRED should be expressed accordingly to the order of magnitude of CLOCKS_PER_SEC
//         if(cycles >= CYCLES_REQUIRED/(FREQUENCY/CLOCKS_PER_SEC)) break;

//         num_runs *= 2;
//     }
// #endif

//     start = clock();
//     for(i=0; i<num_runs; ++i) { 
//         compute(A, B, n);
//     }
//     end = clock();

//     return (double)(end-start)/num_runs;
// }










// #define NUM_RUNS 100
// #ifdef __x86_64__
// double rdtsc() {
//     int i, num_runs;
//     myInt64 cycles;
//     myInt64 start;
//     num_runs = NUM_RUNS;

//     start = start_tsc();
//     for (i = 0; i < num_runs; ++i) {
//         // doo
//     }

//     cycles = stop_tsc(start)/num_runs;
//     return (double) cycles;
// }
// #endif





// Timing function based on the TimeStep Counter of the CPU.
// #ifdef __x86_64__
// double rdtsc(double A[], double B[], int n) {
//     int i, num_runs;
//     myInt64 cycles;
//     myInt64 start;
//     num_runs = NUM_RUNS;

//     /* 
//      * The CPUID instruction serializes the pipeline.
//      * Using it, we can create execution barriers around the code we want to time.
//      * The calibrate section is used to make the computation large enough so as to 
//      * avoid measurements bias due to the timing overhead.
//      */
// #ifdef CALIBRATE
//     while(num_runs < (1 << 14)) {
//         start = start_tsc();
//         for (i = 0; i < num_runs; ++i) {
//             compute(A, B, n);
//         }
//         cycles = stop_tsc(start);

//         if(cycles >= CYCLES_REQUIRED) break;

//         num_runs *= 2;
//     }
// #endif

//     start = start_tsc();
//     for (i = 0; i < num_runs; ++i) {
//         compute(A, B, n);
//     }

//     cycles = stop_tsc(start)/num_runs;
//     return (double) cycles;
// }
// #endif






// #ifndef WIN32
// double timeofday(double A[], double B[], int n) {
//     int i, num_runs;
//     double cycles;
//     struct timeval start, end;

//     num_runs = NUM_RUNS;
// #ifdef CALIBRATE
//     while(num_runs < (1 << 14)) {
//         gettimeofday(&start, NULL);
//         for (i = 0; i < num_runs; ++i) {
//             compute(A, B, n);
//         }
//         gettimeofday(&end, NULL);

//         cycles = (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)*FREQUENCY;

//         if(cycles >= CYCLES_REQUIRED) break;

//         num_runs *= 2;
//     }
// #endif

//     gettimeofday(&start, NULL);
//     for(i=0; i < num_runs; ++i) {
//         compute(A, B, n);
//     }
//     gettimeofday(&end, NULL);

//     return (double)((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec)/1e6)/ num_runs;
// }

// #else

// #endif