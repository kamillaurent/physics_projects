#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<ctime>

using namespace std;

# define SECTION 4 // Define the part of the code that you want to use

void dYdt(double, double *, double *);
void RK4(double, double *, void(double, double *, double *), double, double);
void Boris(double *, double *, double *, double *, double, int);

int main(){
    cout<< setiosflags(ios::scientific)<<setprecision(12);
    ofstream fdata;
    ofstream flog;

    // define the initial conditions and the constants
    double t0 = 0.0; // Initial time
    double tf = 100.0; // Final time
    double dt = 0.1; // Time step
    int Neq = 14; // Number of equations
    int i;
    double t;
    
    // Analytical formula for a known case (particle in B = (0, 0, 1))
    double x_pos, y_pos, z_pos;
    t = 0.0;
    fdata.open("analytical_trajectory.dat");
    for (i = 0; i <= 10000; i ++) {
        // Compute the current time and state to the output file
        x_pos = 0.1 * cos(t)  -0.1;
        y_pos = 0.1 * sin(t);
        z_pos = 0.1 * t;
        fdata << t << " " << x_pos << " " << y_pos << " " << z_pos << endl;
        t += 0.01;
    }
    fdata.close();

    #if SECTION == 1
    // Initial conditions
    double Y[14] = {0.0, 0.0, 0.0, 0.0, 0.1, 0.1, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    double err_x;

    // Open the output file
    fdata.open("particle_1.dat");
    // Loop to perform integration steps
    for (t = t0; t <= tf; t += dt) {
        // Print the current time and state to the output file
        fdata << t << " ";
        for (i = 0; i < Neq; i++) {
            fdata << Y[i] << " ";
        }
        fdata << endl;

        // Call the Runge-Kutta method to update the state
        RK4(t, Y, dYdt, dt, Neq);
    }
    fdata.close();

    
    //Open the "error" output file
    fdata.open("error_RK4.dat");
    double x_tf;

    // Compute the system's error using different time steps
    for (i = 128; i <= 16384; i *= 2){
        // Reset initial conditions
        Y[0] = Y[1] = Y[2] = Y[3] = 0.0;
        Y[4] = Y[5] = 0.1;
        t = t0;
        dt = fabs(tf - t0)/double(i); // Define a new time step

        for (t = t0; t < tf; t += dt) {
            // Call the Runge-Kutta method to update the state
            RK4(t, Y, dYdt, dt, Neq);
        }
        // Print the error at t = tf using i time steps
        x_tf = 0.1 * cos(t) - 0.1;
        fdata << dt << " " << fabs(Y[0] - x_tf) << " " << Y[0] << " " << x_tf << endl;
    }
    fdata.close();

    # elif SECTION == 2

    int N = 3;  // Number of dimensions
    double x[N] = {0.0};  // Initial position
    double v[N] = {0.0, 0.1, 0.1};  // Initial velocity
    double E[N] = {0.0, 0.0, 0.0};  // Electric field 
    double B[N] = {0.0, 0.0, 1.0};  // Magnetic field
    
    // Check if the algorithm works with a known case
    fdata.open("particle_2.dat");
    for (t = t0; t<tf; t+=dt){

        // Print current state in the output file
        fdata << t << " " << x[0] << " " << x[1] << " " << x[2] << " "<< v[0] << " " << v[1] << " " << v[2] << endl;

        // Call the Boris integrator
        Boris(x, v, E, B, dt, N);
    }
    fdata.close();

    //Open the "error" output file
    fdata.open("error_Boris.dat");
    double x_tf;

    // Compute the error for different time steps
    for (i = 128; i <= 16384; i *= 2){
        // Reset initial conditions
        x[0] = x[1] = x[2] = v[0] = 0.0;  
        v[1] = v[2] = 0.1;
        t = t0;
        dt = fabs(tf - t0)/double(i); // Define a new time step

        for (t = t0; t < tf; t += dt) {
            // Call the Boris method to update the state
            Boris(x, v, E, B, dt, N);
        }
        // Print the error at t = tf using i time steps
        x_tf = 0.1 * cos(t) - 0.1;
        fdata << dt << " " << fabs(x[0]-x_tf) << " " << x[0] << " " << x_tf << endl;
    }
    fdata.close();

    #elif SECTION == 3

    // Initial conditions
    double Y[14] = {0.0, 0.0, 0.0, 0.0, 0.1, 0.0, -1.0, 1.0, 0.0, 0.1, 0.0, 0.0, 0.0, 1.0};
     
    // Open the output file
    fdata.open("particle_3.dat");

    // Loop to perform integration steps
    for (double t = t0; t <= tf; t += dt) {
        // Print the current time and state to the output file
        fdata << t << " ";
        for (int i = 0; i < Neq; i++) {
            fdata << Y[i] << " ";
        }
        fdata << sqrt(Y[4]*Y[4] + Y[3]*Y[3]) << endl;

        // Call the Runge-Kutta method to update the state
        RK4(t, Y, dYdt, dt, Neq);
    }
    fdata.close();

    // Analytical solution using a new system X' in motion with v_x = 0.1
    double x, y;
    double v_p;       // Velocity of the particle in the system X'
    double B_p = Y[13] + 0.1 * Y[9]; // Magnetic field in the system X'
    double r_L = v_p/1.01; //Larmor radius in the system x'
    t = 0.0;
    v_p = 0.1 * sqrt(2.0);

    fdata.open("analytical_trajectory_2.dat");
    for (i = 0; i <= 200; i ++) {
        // Compute the current time and state to the output file
        x = r_L * cos(t + M_PI/4.0) - r_L/sqrt(2.0) + 0.1 * t;
        y = r_L * sin(t + M_PI/4.0) - r_L/sqrt(2.0);
        fdata << t << " " << x << " " << y << endl;
        t += dt;
    }
    fdata.close();

    // Compute the error for t from 0 to 20 s
    fdata.open("error_RK.dat");
    // Reset initial conditions
    Y[0] = Y[1] = Y[2] = Y[3] = 0.0;
    Y[4] = Y[5] = 0.1;
    t = t0;

    for (i = 0; i <= 200; i++){
        // Print the error at current time
        x = r_L * cos(t + M_PI / 4.0) - r_L / sqrt(2.0) + 0.1 * t;
        fdata << t << " " << fabs(Y[0] - x) << endl;

        // Call the RK4 function to update the state and update time
        RK4(t, Y, dYdt, dt, Neq);
        t += dt;
    }
    fdata.close();

    #elif SECTION == 4

    int N = 3;
    int p; 
    double x[N], v[Neq]; // Position and Velocity
    double E[N] = {0.0, 0.0, 0.5}; //Electic field
    double B[N];   // Magnetic fields
    double v0 = 0.1; 
    double energy; //total energy of the particle  
    srand48(time(NULL));
    
    fdata.open("particle_1000.dat");
    flog.open("energy.dat");

    // Loop over p=1000 particle
    for(p=0; p<1000; p++){
        // Define random angles theta in [0, Pi] and phi in [0, 2Pi]
        double theta, phi;
        theta = drand48() * M_PI;
        phi = drand48() * 2 * M_PI;

        // Define position and velocity of one particle
        x[0] = (drand48() - 0.5) * 2000.0;
        x[1] = (drand48() - 0.5) * 2000.0;
        x[2] = 0.0;

        v[0] = v0 * sin(theta) * cos(phi);
        v[1] = v0 * sin(theta) * sin(phi);
        v[2] = v0 * cos(theta);
    
        // Define initial Magnetic fields
        B[0] = 0.5 * x[1];
        B[1] = 0.5 * x[0];
        B[2] = 0.0;
    
        for (double t = t0; t<tf; t+=dt){

           // Print current state in the output file
           fdata << p << " " << t << " " << x[0] << " " << x[1] << " " << x[2];
           fdata << " " << v[0] << " " << v[1] << " " << v[2] << endl;

           // Call the Boris integrator
           Boris(x, v, E, B, dt, N);

           // Compute the magnetic field for the new position
           B[0] = x[1] / 1000.0;
           B[1] = x[0] / 1000.0;
           B[2] = 0.0;
        }
        fdata << endl;
        fdata << endl;

        // Compute the energy of the particle and write it to an output file
        energy = 0.5 * v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        flog << energy << endl;
    }

    fdata.close();
    flog.close();
    #endif

    return 0;

}

void RK4(double t, double * Y, void(*RHSFunc)(double, double *, double *), double dt, double Neq){
    int i;
    double k1[14], k2[14], k3[14], k4[14], Y1[14], Y2[14], Y3[14];

    RHSFunc(t, Y, k1); // compute k1
    for(i = 0; i < Neq; i++) Y1[i] = Y[i] + 0.5 * dt * k1[i];

    RHSFunc(t + 0.5 * dt, Y1, k2); // compute k2
    for(i = 0; i < Neq; i++) Y2[i] = Y[i] + 0.5 * dt * k2[i];

    RHSFunc(t + 0.5 * dt, Y2, k3); // compute k3
    for(i = 0; i < Neq; i++) Y3[i] = Y[i] + dt * k3[i];

    RHSFunc(t + dt, Y3, k4); // compute k4
    for(i = 0; i < Neq; i++) Y[i] += dt * (k1[i] + 2
     * k2[i] + 2 * k3[i] + k4[i]) / 6.0; //final result of the differential equation
}

void dYdt (double t, double *Y, double *R){
    double x  = Y[0]; 
    double y  = Y[1];
    double z  = Y[2]; //(x,y,z) = Position
    double vx = Y[3];
    double vy = Y[4];
    double vz = Y[5]; // (vx, vy, vz) = Velocity
    double q  = Y[6]; // Charge
    double m  = Y[7]; // Mass
    double Ex = Y[8];
    double Ey = Y[9];
    double Ez = Y[10]; // Electric Field
    double Bx = Y[11];
    double By = Y[12];
    double Bz = Y[13]; // Magnetic Field


    R[0] = vx; // Solution for each Variable dY/dt= R(Y,t)
    R[1] = vy;
    R[2] = vz;
    R[3] = (q / m) * (Ex + vy * Bz - vz * By);
    R[4] = (q / m) * (Ey + vz * Bx - vx * Bz);
    R[5] = (q / m) * (Ez + vx * By - vy * Bx);
    R[6] = 0.0;   
    R[7] = 0.0;
    R[8] = 0.0;
    R[9] = 0.0;
    R[10]= 0.0;
    R[11]= 0.0;
    R[12]= 0.0;
    R[13]= 0.0; // Each constant is put in the form of (dy/dt=0), y=constant

}

void Boris(double *x, double *v, double *E, double *B, double dt, int N){
    int i;
    double m = 1.0;   // Mass
    double q = -1.0;  // Charge
    double gamma;     // Lorentz factor
    double h = (q / m) * dt;
    double b_squared;
    double v_squared;

    double b[N]; 
    double u_minus[N];
    double u_plus[N];
    double u[N];
    double uxb[N];      
    
    //Current velocity squared and Lorentz factor
    v_squared = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    gamma = 1.0; // we define gamma = 1 for the non relativistic formula.
                 // For the relativistic case gamma = sqrt(1.0/(1.0-v_squared));

    //Modify the magnetic field vector
    for (i = 0; i < N; i++) b[i] = h * 0.5 * B[i] / gamma; 

    // Square the modified magnetic field vector
    b_squared = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];

    // Compute position using current velocity
    for (i=0; i< N; i++) x[i] += 0.5 * dt * v[i];

    // Compute velocity using acceleration generated by electric field (kick)
    for (i=0; i< N; i++) u_minus[i] = v[i] * gamma + 0.5 * h * E[i];

    // Compute the cross product (u_minus x b)
    uxb[0] = u_minus[1] * b[2] - u_minus[2] * b[1];
    uxb[1] = u_minus[2] * b[0] - u_minus[0] * b[2];
    uxb[2] = u_minus[0] * b[1] - u_minus[1] * b[0];

    // Add the rotation due to magnetic field (rotation)
    u_plus[0] = u_minus[0] + 2*(uxb[0] + uxb[1]*b[2]-uxb[2]*b[1])/(1+b_squared);
    u_plus[1] = u_minus[1] + 2*(uxb[1] + uxb[2]*b[0]-uxb[0]*b[2])/(1+b_squared);
    u_plus[2] = u_minus[2] + 2*(uxb[2] + uxb[0]*b[1]-uxb[1]*b[0])/(1+b_squared);
    
    // Compute current velocity and position using electric field and velocity
    for (i=0; i< N; i++) v[i] = (u_plus[i] + 0.5 * h * E[i]) / gamma;
    for (i=0; i< N; i++) x[i] += 0.5 * dt * v[i];
}