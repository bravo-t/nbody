#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdbool.h>
#include "thread_barrier.h"
#include "thread_control.h"
 
typedef struct {
	double x,y,z;
} vector;

typedef struct {
    int id;
    int number_of_threads;
    double ts;
    thread_barrier_t* acc_barrier;
    thread_barrier_t* vel_barrier;
    thread_barrier_t* pos_barrier;
    thread_barrier_t* col_barrier;
    ThreadControl* handle;
    double* mass;
    vector* acc;
    vector* vel;
    vector* pos;
} SlaveArgs;
 
int bodies = 0;
int timeSteps = 0;
float timescale = 0.01;
double GravConstant = 1;
double* masses;
vector* positions;
vector* velocities;
vector* accelerations;

int calc_i_start(int id, int bodies,int number_of_threads) {
    if (bodies < number_of_threads) {
        if (id == 0) {
            return 0;
        } else {
            return bodies + 1;
        }
    } else {
        return (id*bodies/number_of_threads);
    }
}

int calc_i_end(int id, int bodies,int number_of_threads) {
    int end = ((id+1)*bodies/number_of_threads)-1;
    if (bodies < number_of_threads) {
        if (id == 0) {
            return bodies - 1;
        } else {
            return -1;
        }
    } else {
        if (end < bodies) {
            return end;
        } else {
            // Return a value that makes no sense to prevent functions from over-writing
            return -1;
        }
    }
}

double get_time() {
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t,&tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

void timer(double t, char* msg) {
    double current_t = get_time();
    double elapsed_time = current_t - t;
    printf("Time elapsed in %s: %f\n",msg,elapsed_time);
}

 
vector addVectors(vector a,vector b){
	vector c = {a.x+b.x,a.y+b.y,a.z+b.z};
	return c;
}
 
vector scaleVector(double b,vector a){
	vector c = {b*a.x,b*a.y,b*a.z};
	return c;
}
 
vector subtractVectors(vector a,vector b){
    double x = a.x-b.x;
    double y = a.y-b.y;
    double z = a.z-b.z;
    vector c;
    if (x == 0 && y == 0 && z == 0) {
        c.x = 1e-6;
        c.y = 1e-6;
        c.z = 1e-6;
    } else {
        c.x = x;
        c.y = y;
        c.z = z;
    }
	return c;
}
 
double mod(vector a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}
 
void initiateSystem(char* fileName){
	int i;
	FILE* fp = fopen(fileName,"r");

	masses = (double*)malloc(bodies*sizeof(double));
	positions = (vector*)malloc(bodies*sizeof(vector));
	velocities = (vector*)malloc(bodies*sizeof(vector));
	accelerations = (vector*)malloc(bodies*sizeof(vector));
 
	for(i=0;i<bodies;i++){
		fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf",&masses[i],&(positions[i].x),&(positions[i].y),&(positions[i].z),&(velocities[i].x),&(velocities[i].y),&(velocities[i].z));
	}
 
	fclose(fp);
}

void calcVelocityAfterCollision(vector* v1, vector* v2, double m1, double m2) {
    double m1s2 = m1 - m2;
    double m1p2 = m1 + m2;
    double v1x = v1->x;
    double v2x = v2->x;
    double v1y = v1->y;
    double v2y = v2->y;
    double v1z = v1->z;
    double v2z = v2->z;
    v1->x = (2*m2*v2x+v1x*m1s2)/m1p2;
    v2->x = (2*m1*v1x-v2x*m1s2)/m1p2;
    v1->y = (2*m2*v2y+v1y*m1s2)/m1p2;
    v2->y = (2*m1*v1y-v2y*m1s2)/m1p2;
    v1->z = (2*m2*v2z+v1z*m1s2)/m1p2;
    v2->z = (2*m1*v1z-v2z*m1s2)/m1p2;
}
 
void resolveCollisions(double* masses, vector* positions,vector* velocities, int i_start, int i_end){
	int i,j;
 
	for(i=i_start;i<i_end;i++)
		for(j=i+1;j<bodies;j++){
			if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){
                calcVelocityAfterCollision(&(velocities[i]),&(velocities[j]),masses[i],masses[j]);
			}
		}
}
 
void computeAccelerations(double* masses, vector* positions, vector* accelerations, int i_start, int i_end){
	int i,j;
 
	for(i=i_start;i<=i_end;i++){
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for(j=0;j<bodies;j++){
			if(i!=j){
				accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*masses[j]/pow(mod(subtractVectors(positions[i],positions[j])),3),subtractVectors(positions[j],positions[i])));
			}
		}
        //fprintf(stderr,"A: %d,%f,%f,%f\n",i,accelerations[i].x,accelerations[i].y,accelerations[i].z);
	}
}
 
void computeVelocities(vector* velocities, vector* accelerations, double timescale, int i_start, int i_end){
	int i;
	for(i=i_start;i<=i_end;i++) {
		velocities[i] = addVectors(velocities[i],scaleVector(timescale,accelerations[i]));
        //fprintf(stderr,"V: %d,%f,%f,%f\n",i,velocities[i].x,velocities[i].y,velocities[i].z);
    }
}
 
void computePositions(vector* positions,vector* velocities, vector* accelerations,double timescale, int i_start, int i_end){
	int i;
	for(i=i_start;i<=i_end;i++) {
        vector v = addVectors(velocities[i],scaleVector(0.5*timescale,accelerations[i]));
		positions[i] = addVectors(positions[i],scaleVector(timescale,v));
        //fprintf(stderr,"P: %d,%f,%f,%f\n",i,positions[i].x,positions[i].y,positions[i].z);
    }
}
 
void simulator_slave(
        double* masses,
        vector* positions,
        vector* velocities,
        vector* accelerations,
        double timescale,
        thread_barrier_t* acc_barrier,
        thread_barrier_t* vel_barrier,
        thread_barrier_t* pos_barrier,
        thread_barrier_t* col_barrier,
        int number_of_threads,
        int i_start,
        int i_end) {
    //double start_time = get_time();
    //printf("DEBUG: Calc accelaration\n");
	computeAccelerations(masses,positions,accelerations,i_start,i_end);
    thread_barrier_wait_reinit(acc_barrier,number_of_threads);
    //timer(start_time,"computeAccelerations");
    //start_time = get_time();
    //printf("DEBUG: Calc position\n");
	computePositions(positions, velocities, accelerations,timescale,i_start,i_end);
    thread_barrier_wait_reinit(pos_barrier,number_of_threads);
    //timer(start_time,"computePositions");
    //start_time = get_time();
    //printf("DEBUG: Calc velocity\n");
	computeVelocities(velocities, accelerations,timescale,i_start,i_end);
    thread_barrier_wait_reinit(vel_barrier,number_of_threads);
    //timer(start_time,"computeVelocities");
    //start_time = get_time();
    //printf("DEBUG: resolve collision\n");
	resolveCollisions(masses, positions, velocities,i_start,i_end);
    thread_barrier_wait_reinit(col_barrier,number_of_threads);
    //timer(start_time,"resolveCollisions");
}

void* simulator(void* args) {
    SlaveArgs* a = (SlaveArgs*) args;

    double* masses = a->mass;
    vector* positions = a->pos;
    vector* velocities = a->vel;
    vector* accelerations = a->acc;
    double timescale = a->ts;
    int id = a->id;
    int number_of_threads = a->number_of_threads;
    thread_barrier_t* acc_barrier = a->acc_barrier;
    thread_barrier_t* vel_barrier = a->vel_barrier;
    thread_barrier_t* pos_barrier = a->pos_barrier;
    thread_barrier_t* col_barrier = a->col_barrier;
    ThreadControl* handle  = a->handle;

    int i_start = calc_i_start(id,bodies,number_of_threads);
    int i_end = calc_i_end(id,bodies,number_of_threads);
    while(1) {
        threadController_slave(handle,CONTROL_WAIT_INST);
        //printf("DEBUG: Run simulation on %d to %d\n",i_start,i_end);
        simulator_slave(masses,positions,velocities,accelerations,timescale,acc_barrier,vel_barrier,pos_barrier,col_barrier,number_of_threads,i_start,i_end);
        threadController_slave(handle,CONTROL_EXEC_COMPLETE);
    }
}

void run_simulation(ThreadControl* handle) {
    threadController_master(handle,THREAD_RESUME);
}
 
int main(int argC,char* argV[])
{
	int i,j;

    if (argC != 7) {
        fprintf(stderr, "Incorrect parameters\n");
        return 0;
    }
    bodies = atoi(argV[1]);
    timeSteps = atoi(argV[2]);
    GravConstant = atof(argV[3]);
    timescale = atof(argV[4]);
    int print_step = atoi(argV[5]);
    int number_of_threads = atoi(argV[6]);
	initiateSystem("data.ssv");

    thread_barrier_t acc_barrier = THREAD_BARRIER_INITIALIZER;
    thread_barrier_init(&acc_barrier,number_of_threads);
    thread_barrier_t vel_barrier = THREAD_BARRIER_INITIALIZER;
    thread_barrier_init(&vel_barrier,number_of_threads);
    thread_barrier_t pos_barrier = THREAD_BARRIER_INITIALIZER;
    thread_barrier_init(&pos_barrier,number_of_threads);
    thread_barrier_t col_barrier = THREAD_BARRIER_INITIALIZER;
    thread_barrier_init(&col_barrier,number_of_threads);

    pthread_mutex_t simulate_control_handle_mutex = PTHREAD_MUTEX_INITIALIZER;
    thread_barrier_t simulate_inst_ready = THREAD_BARRIER_INITIALIZER;
    thread_barrier_t simulate_inst_ack = THREAD_BARRIER_INITIALIZER;
    thread_barrier_t simulate_thread_complete = THREAD_BARRIER_INITIALIZER;
    ThreadControl* simulate_control_handle = initControlHandle(&simulate_control_handle_mutex, &simulate_inst_ready, &simulate_inst_ack, &simulate_thread_complete, number_of_threads);

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_t* simulator_thread = malloc(sizeof(pthread_t)*number_of_threads);

    SlaveArgs** simulator_arguments = malloc(sizeof(SlaveArgs*)*number_of_threads);

    for(int i=0;i<number_of_threads;i++) {
        simulator_arguments[i] = (SlaveArgs*) malloc(sizeof(SlaveArgs));
        
        (simulator_arguments[i])->mass = masses;
        (simulator_arguments[i])->acc = accelerations;
        (simulator_arguments[i])->vel = velocities;
        (simulator_arguments[i])->pos = positions;
        (simulator_arguments[i])->ts = timescale;
        (simulator_arguments[i])->id = i;
        (simulator_arguments[i])->number_of_threads = number_of_threads;
        (simulator_arguments[i])->acc_barrier = &acc_barrier;
        (simulator_arguments[i])->vel_barrier = &vel_barrier;
        (simulator_arguments[i])->pos_barrier = &pos_barrier;
        (simulator_arguments[i])->col_barrier = &col_barrier;
        (simulator_arguments[i])->handle = simulate_control_handle;

        int create_thread_error;
        create_thread_error = pthread_create(&simulator_thread[i],&attr,simulator,simulator_arguments[i]);
        if (create_thread_error) {
            printf("Error happened while creating slave threads\n");
            exit(-1);
        }
        
    }

	for(i=0;i<timeSteps;i++){
        if (i % 1000 == 0) {
            fprintf(stderr, "Iteration %d\n",i);
        }
        int p;
        for(p=0;p<(int) 1/timescale;p++) {
            //printf("DEBUG: Start to run simulation\n");
		    run_simulation(simulate_control_handle);
            //printf("DEBUG: simulation finished\n");
            
            if (p%print_step == 0) {
                //printf("%d ",i+1);
                for(j=0;j<bodies;j++) {
                    printf("%f %f %f ",positions[j].x,positions[j].y,positions[j].z);
                }
                printf("\n");
            }
        }
	}
	return 0;
}
 
