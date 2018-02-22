/**
 * Solar System
 *
 * This example integrates all planets of the Solar
 * System. The data comes from the NASA HORIZONS system. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

double ss_mass[10] = { // GM=1 with distance in AU and time in days
    2.9591220828559E-04,
    4.9125495718311E-11,
    7.2434531799395E-10,
    8.8876925870232E-10,
    1.0931896589621E-11,
    9.5495401748807E-11,
    2.8247604533652E-07,
    8.4576151711856E-08,
    1.2918949220207E-08,
    1.5240407045482E-08,
};

double ss_radius[10] = { // In AU
    4.6544780132355e-03,
    1.630705028477387448e-05,
    4.045512126396863864e-05,
    4.263429666582814611e-05,
    1.161179629009251657e-05,
    2.264069658312322796e-05,
    4.778945025452157572e-04,
    4.028666966848747089e-04,
    1.708513622580592106e-04,
    1.655371154958557966e-04,
};

double ss_pos[10][3] = { // In AU
    {0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00},
    {2.303471432830715981e-01, -3.641428643096443207e-01, -5.088713652760864375e-02},
    {6.777156115147005000e-01, -2.615440985251679118e-01, -4.269626038165600518e-02},
    {-7.662129455848691872e-01, 6.217897740204972878e-01, -2.570841958077960370e-05},
    {-7.668303568173396867e-01, 6.191671210294678040e-01, 1.861475764952425535e-04},
    {-1.344650640881984938e+00, -8.540472053528166407e-01, 1.510388387588550325e-02},
    {-4.071390698448524859e+00, -3.585036625965199342e+00, 1.059900141150976194e-01},
    {2.571976539468074918e-01, -1.006066100820767062e+01, 1.646307186561387359e-01},
    {1.764841201673070259e+01, 9.189142828845191957e+00, -1.943915369997064602e-01},
    {2.871500671858860443e+01, -8.476877557412295872e+00, -4.872526183251342791e-01},
};

double tesla_pos[3] = {
    -7.706431161491976711e-01, 6.173221163760826968e-01, -1.279227125742519663e-03
};

double ss_vel[10][3] = { // In AU/day
    {0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00},
    {1.815986862482058678e-02, 1.642299178819348754e-02, -3.240064549949408135e-04},
    {7.170231292269824662e-03, 1.878201481396940889e-02, -1.561281282806556360e-04},
    {-1.112688908283271383e-02, -1.342184146093371340e-02, 9.307833963532550489e-07},
    {-1.058668458179885287e-02, -1.356343231772417363e-02, -2.544213039788941757e-05},
    {8.026699294846572491e-03, -1.061652753517928344e-02, -4.194462999723903201e-04},
    {4.900315299659889490e-03, -5.313130356766265945e-03, -8.760113064667799580e-05},
    {5.277312551390259121e-03, 1.206441060316985772e-04, -2.120796839867706239e-04},
    {-1.839947612523925875e-03, 3.300033649428975400e-03, 3.608181382823715021e-05},
    {8.737578121952224979e-04, 3.025178809253320500e-03, -8.248249054337306659e-05},
};

double tesla_vel[3] = {
    -1.251899167082707286e-02, -1.488669486874137464e-02, -3.766015398112917034e-04
};

void heartbeat(struct reb_simulation* r);
double e_init;
double tmax;

double ranf(void);
double randGauss(void);
void ranVec(double vec[3]);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    double tes_earth_pos[3];
    double tes_earth_vel[3];
    double dist, speed;
    double dfact = 0.00001, vfact = 0.0001;
    double interval = 365.25 / 30.0; // 1 year per second, 30 FPS
    int doGrab = 0;
    int i;

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-grab") == 0) {
            doGrab = 1;
        }
    }
    if (!doGrab) {
        interval = 0.0;
    }

    // Setup constants
    r->dt = 1;                   // in days
    tmax = 73050 + interval*0.5; // 200 years plus slop
    r->G = 1;                    // All other units corrected to this standard
    r->N_active = 10;            // All other particles are test particles
    r->usleep = 2000000. / interval; // It takes 1.8 seconds to dump a grab
    r->ri_whfast.safe_mode = 0;  // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs. 
    r->ri_whfast.corrector = 11; // 11th order symplectic corrector
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->heartbeat = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    //r->integrator        = REB_INTEGRATOR_IAS15;        // Alternative non-symplectic integrator
    if (doGrab) {
        r->start_paused = 1;
        r->screenDumpInterval = interval;
        r->screenDumpPath = "GRABS/tesla%05d.png";
    }

    // Initial conditions
    for (int i = 0; i < 10; i++ ){
        struct reb_particle p = {0};
        p.x  = ss_pos[i][0];         p.y  = ss_pos[i][1];         p.z  = ss_pos[i][2];
        p.vx = ss_vel[i][0];         p.vy = ss_vel[i][1];         p.vz = ss_vel[i][2];
        p.m  = ss_mass[i];
        p.r = ss_radius[i];
        reb_add(r, p); 
    }

    //srandomdev();

    /* Get Earth-relative position and velocity */
    tes_earth_pos[0] = tesla_pos[0] - ss_pos[3][0];
    tes_earth_pos[1] = tesla_pos[1] - ss_pos[3][1];
    tes_earth_pos[2] = tesla_pos[2] - ss_pos[3][2];

    tes_earth_vel[0] = tesla_vel[0] - ss_vel[3][0];
    tes_earth_vel[1] = tesla_vel[1] - ss_vel[3][1];
    tes_earth_vel[2] = tesla_vel[2] - ss_vel[3][2];

    dist = tes_earth_pos[0]*tes_earth_pos[0] +
           tes_earth_pos[1]*tes_earth_pos[1] +
           tes_earth_pos[2]*tes_earth_pos[2];
    dist = sqrt(dist);
printf("dist = %g\n", dist);

    speed = tes_earth_vel[0]*tes_earth_vel[0] +
            tes_earth_vel[1]*tes_earth_vel[1] +
            tes_earth_vel[2]*tes_earth_vel[2];
    speed = sqrt(speed);
printf("speed = %g\n", speed);

    // Test particles
    for (int i = 0; i < 1000; i++) {
        struct reb_particle p = {0};
        double vec[3], fact;

        /* Position */
        ranVec(vec);
        fact = randGauss(); if (fact < 0.) fact = -fact;
        fact *= dfact * dist * 0.5;
        vec[0] *= fact; vec[1] *= fact; vec[2] *= fact;
        p.x = tesla_pos[0] + vec[0];
        p.y = tesla_pos[1] + vec[1];
        p.z = tesla_pos[2] + vec[2];

        /* Velocity */
        ranVec(vec);
        fact = randGauss(); if (fact < 0.) fact = -fact;
        fact *= vfact * speed * 0.5;
        vec[0] *= fact; vec[1] *= fact; vec[2] *= fact;
        p.vx = tesla_vel[0] + vec[0];
        p.vy = tesla_vel[1] + vec[1];
        p.vz = tesla_vel[2] + vec[2];

        reb_add(r, p); 
    }
 
    reb_move_to_com(r);
    e_init = reb_tools_energy(r);
    system("rm -f energy.txt");
    reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 10000.)){
        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);
        FILE* f = fopen("energy.txt","a");
        double e = reb_tools_energy(r);
        fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
        fclose(f);
    }
}

/* Generate a random double from 0 to 1 inclusive */
double ranf(void)
{
    long long hi, lo;

    hi = random();
    lo = random();
    return ((hi << 31LL) | lo) * (1. / (double)0x3FFFFFFFFFFFFFFFLL);
}

/* Generate a random direction vector uniformly distributed on
 * the positive unit sphere octant
 */
void ranVec(double vec[3])
{
    double dist;

    do {
        vec[0] = 2.0 * ranf() - 1.0;
        vec[1] = 2.0 * ranf() - 1.0;
        vec[2] = 2.0 * ranf() - 1.0;
        dist = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
    } while ((dist < 0.00000001) || (dist > 1.0));
    dist = 1.0 / sqrt(dist);
    vec[0] *= dist;
    vec[1] *= dist;
    vec[2] *= dist;
}

/* Return a randum number from a standard distribution with mean
 * zero and standard deviation one
 */
double randGauss(void)
{
    double x1, x2, w;
    static int parity = 1;
    static double y1, y2;

    if (parity) {
        do {
            x1 = 2.0 * ranf() - 1.0;
            x2 = 2.0 * ranf() - 1.0;
            w = x1*x1 + x2*x2;
        } while (w >= 1.0);

        w = sqrt((-2.f*log(w)) / w);
        y1 = x1*w;
        y2 = x2*w;
        parity = 0;
        return y1;
    }
    else {
        parity = 1;
        return y2;
    }
}

