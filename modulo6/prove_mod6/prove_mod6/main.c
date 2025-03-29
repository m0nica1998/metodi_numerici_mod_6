#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define nx 40
#define nt 8
#define mass 0.04
#define pigr 3.141592653
#define nvol (nx * nt)

float field[nx][nt];
int npp[nx][2], nmm[nx][2];
float  mass2, mass2p4;
long idum = 0, idum2 = 0, iy = 0, iv[NTAB];

float ran2(void);
void ranstart(void);
void ranfinish(void);

void geometry(void);
void initialize_lattice(int iflag);
void energy(float *xene_mass, float *xene_spat, float *xene_temp);
void update_heatbath(void);
void update_overrelax(void);


int main(void) {
    FILE *input_file, *output_file, *lattice_file;
    int iflag, measures, i_decorrel, iter;
    double mass2, mass2p4;
    double xene_mass, xene_spat, xene_temp;
    clock_t start_time, end_time;
    double cpu_time_used;
    int count = 0;
    double xene_mass_sum = 0.0, xene_spat_sum = 0.0, xene_temp_sum = 0.0;

    // Open input file for reading simulation parameters
    input_file = fopen("/Users/monicacesario/Desktop/input6.txt", "r");
    if (input_file == NULL) {
        printf("Error opening input file.\n");
        return 1;
    }

    // Open output file for writing measurements
    output_file = fopen("/Users/monicacesario/Desktop/modulo6/ene_cont_0.05.dat", "w");
    if (output_file == NULL) {
        printf("Error opening output file.\n");
        return 1;
    }

    // Read simulation parameters from the input file
    fscanf(input_file, "%d", &iflag);
    fscanf(input_file, "%d", &measures);
    fscanf(input_file, "%d", &i_decorrel);
    
    fclose(input_file);

    mass2 = pow(mass,2);
    mass2p4 = pow(mass,2) + 4.0;

    geometry();
    initialize_lattice(iflag);

    start_time = clock();

    for (int i = 0; i < measures; i++) {
        
        update_heatbath();
        update_overrelax();
        update_overrelax();
        update_overrelax();
        update_overrelax();
        // Update spin configuration
        if((i+1) % i_decorrel == 0){
            
            // Measure physical observables
            energy(&xene_mass, &xene_spat, &xene_temp);
            
            // Write measurements to the output file
            
            fprintf(output_file, "%d %f\n", i, xene_mass);
            count +=1;
            xene_mass_sum += xene_mass;
            xene_spat_sum += xene_spat;
            xene_temp_sum += xene_temp;
        }
        
        
    }
    printf("xene_mass = %f, xene_spat = %f, xene_temp = %f\n", xene_mass_sum/count, xene_spat_sum/count, xene_temp_sum/count);
    end_time = clock();

    fclose(output_file);

    // Save configuration and random number generator state
    lattice_file = fopen("lattice", "w");
    if (lattice_file == NULL) {
        printf("Error opening lattice file.\n");
        return 1;
    }
    // Non ho capito cosa va scritto nel file lattice !!


    fclose(lattice_file);

    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    printf("Durata del programma: %.2f secondi\n", cpu_time_used);

    

    return 0;
}

void geometry(void) {
    int i;
    for (i = 0; i < nx; i++) {
        npp[i][0] = i + 1;
        nmm[i][0] = i - 1;
    }
    npp[nx - 1][0] = 0;
    nmm[0][0] = nx - 1;

    for (i = 0; i < nt; i++) {
        npp[i][1] = i + 1;
        nmm[i][1] = i - 1;
    }
    npp[nt - 1][1] = 0;
    nmm[0][1] = nt - 1;
}

void initialize_lattice(int iflag) {
    int ix, it;
    float x;

    if (iflag == 0) {
        for (ix = 0; ix < nx; ix++) {
            for (it = 0; it < nt; it++) {
                field[ix][it] = 0.0;
            }
        }
    } else if (iflag == 1) {
        for (ix = 0; ix < nx; ix++) {
            for (it = 0; it < nt; it++) {
                x = (float)ran2();
                printf("nr random nell'initialize %f\n", x);
                field[ix][it] = 1.0 - 2.0 * x;
            }
        }
    }
}

void energy(float *xene_mass, float *xene_spat, float *xene_temp) {
    int ix, it;
    float phi, force_s, force_t;
    *xene_mass = 0.0;
    *xene_spat = 0.0;
    *xene_temp = 0.0;

    for (ix = 0; ix < nx; ix++) {
        for (it = 0; it < nt; it++) {
            phi = field[ix][it];
            force_s = field[npp[ix][0]][it];
            force_t = field[ix][npp[it][1]];

            *xene_mass += mass2 * phi * phi;
            *xene_spat -= 2.0 * phi * force_s + 2.0 * phi * phi;
            *xene_temp -= 2.0 * phi * force_t + 2.0 * phi * phi;
        }
    }

    *xene_mass /= nvol;
    printf("ene_mass = %f \n",&xene_mass);
    *xene_spat /= nvol;
    *xene_temp /= nvol;
}

void update_heatbath(void) {
    int ix, it;
    float force, phi, sigma2, aver, x, y;

    for (ix = 0; ix < nx; ix++) {
        for (it = 0; it < nt; it++) {
            force = 0.0;
            phi = field[ix][it];
            force += field[npp[ix][0]][it];
            force += field[nmm[ix][0]][it];
            force += field[ix][npp[it][1]];
            force += field[ix][nmm[it][1]];

            sigma2 = 1.0 / mass2p4;
            aver = force * sigma2;

            x = sqrt(sigma2) * sqrt(-2.0 * log(ran2()));
            y = x * cos(2.0 * pigr * ran2()) + aver;
            field[ix][it] = y;
        }
    }
}

void update_overrelax(void) {
    int ix, it;
    float force, phi, aver;

    for (ix = 0; ix < nx; ix++) {
        for (it = 0; it < nt; it++) {
            force = 0.0;
            phi = field[ix][it];
            force += field[npp[ix][0]][it];
            force += field[nmm[ix][0]][it];
            force += field[ix][npp[it][1]];
            force += field[ix][nmm[it][1]];

            aver = force / mass2p4;

            field[ix][it] = 2.0 * aver - phi;
        }
    }
}


float ran2(void)

{

    int j;

    long k;

    float temp;

   
    if (idum <= 0) {

        if (-(idum) < 1) idum=1;

        else idum = -(idum);

        idum2=(idum);

        for (j=NTAB+7;j>=0;j--) {

            k=(idum)/IQ1;

            idum=IA1*(idum-k*IQ1)-k*IR1;

            if (idum < 0) idum += IM1;

            if (j < NTAB) iv[j] = idum;

        }

        iy=iv[0];

    }

    k=(idum)/IQ1;

    idum=IA1*(idum-k*IQ1)-k*IR1;

    if (idum < 0) idum += IM1;

    k=idum2/IQ2;

    idum2=IA2*(idum2-k*IQ2)-k*IR2;

    if (idum2 < 0) idum2 += IM2;

    j=iy/NDIV;

    iy=iv[j]-idum2;

    iv[j] = idum;

    if (iy < 1) iy += IMM1;

    if ((temp=AM*iy) > RNMX) return RNMX;

    else return temp;

}


void ranstart(void) {
    FILE *fp;
    int i;

    fp = fopen("/Users/monicacesario/Desktop/modulo3/path_int/randomseed3.txt", "r");
    if (fp == NULL) {
        fprintf(stderr, "Errore nell'apertura del file randomseed\n");
        exit(1);
    }

    fscanf(fp, "%li", &idum);
    fscanf(fp, "%li", &idum2);
    for (i = 0; i < 32; i++) {
        fscanf(fp, "%li", &iv[i]);
    }
    fscanf(fp, "%li", &iy);

    fclose(fp);
}

void ranfinish(void) {
    FILE *fp;
    int i;

    fp = fopen("/Users/monicacesario/Desktop/modulo3/path_int/randomseed3.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Errore nell'apertura del file randomseed\n");
        exit(1);
    }

    fprintf(fp, "%li\n", idum);
    fprintf(fp, "%li\n", idum2);
    for (i = 0; i < 32; i++) {
        fprintf(fp, "%li\n", iv[i]);
    }
    fprintf(fp, "%li\n", iy);

    fclose(fp);
}

