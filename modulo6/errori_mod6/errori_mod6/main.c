//
//  main.c
//  errori_mod6_simone
//
//  Created by Monica  on 15/07/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SKIP_ROWS 0
#define NUM_RICAMPIONAMENTI 200

const char* INPUT_FILE = "/Users/monicacesario/Desktop/modulo6/termo_160.dat";
const char* OUTPUT_FILE = "/Users/monicacesario/Desktop/obs1_vs_k.dat";
const double M = 0.05;

int contaRighe(const char* nomeFile);
double* readData(int num_righe);
double bootstrap(double *lista, int N, int dim_blocco, int num_ricampionamenti);

int main(void)
{
    int num_righe = contaRighe(INPUT_FILE) - SKIP_ROWS;
    printf("%d\n", num_righe);
    double* data_array = readData(num_righe); //num_righe e' la dimensione di data_array

    //printf("%f %f \n", data_array[0], data_array[1]);

    FILE *file = fopen(OUTPUT_FILE, "w");

    if (file == NULL)
    {
        printf("Errore nell'apertura del file di output.\n");
        exit(1);
    }

    //per il test:
    //int k = 1000;
    //double d_ene = bootstrap(data_array, num_righe, 1000, NUM_RICAMPIONAMENTI);
    //fprintf(file, "%d\t%.6f\n", k, d_ene);

    int k;
    for (k = 3; k <=12 ; k++)
    {
        double d_ene = bootstrap(data_array, num_righe, pow(2, k), NUM_RICAMPIONAMENTI);
        fprintf(file, "%d\t%.6f\n", k, d_ene);
    }


    fclose(file);


    double sum = 0.0;
    int i;
    for (i = 0; i < num_righe; i++)
    {
        sum += data_array[i];
    }
    double average = sum / num_righe;
    printf("%.2f\t%.6f\n", M, average);



    free(data_array);
    return 0;
}

double* readData(int num_righe)
{
    double* result = (double*)malloc(sizeof(double) * num_righe);

    FILE *file = fopen(INPUT_FILE, "r");

    if (file == NULL)
    {
        printf("Errore nell'apertura del file.\n");
        exit(1);
    }

    int line_count = 0;

    while (line_count < SKIP_ROWS)
    {
        char buffer[256];
        fgets(buffer, sizeof(buffer), file);
        line_count++;
    }

    int i = 0;
    double mass, spat, temp;
    while (!feof(file))
    {
        if (fscanf(file, "%*d %lf %lf %lf", &mass, &spat, &temp) == 3)
        {
            if(i >= num_righe){
                printf("Nel file sono state trovate piu' di %d righe. Le successive saranno saltate.\n", i);
                break;
            }
            result[i] = mass + spat - temp;
            i++;
        }
    }

    fclose(file);

    return result;
}


int contaRighe(const char* nomeFile) {
    FILE* file = fopen(nomeFile, "r");
    if (file == NULL) {
        perror("Errore durante l'apertura del file");
        return -1;
    }

    int conteggio = 0;
    char carattere;

    while ((carattere = fgetc(file)) != EOF) {
        if (carattere == '\n') {
            conteggio++;
        }
    }

    fclose(file);
    return conteggio;
}

double bootstrap(double *lista, int N, int dim_blocco, int num_ricampionamenti)
{
    int num_blocchi = N / dim_blocco;

    if (N % dim_blocco != 0)
    {
        printf("Warning: %d does not divide %d. Throwing away last %d elements.\n", dim_blocco, N, N % dim_blocco);
    }

    double medie_campioni[num_ricampionamenti];
    int num_ricampionamento;
    for (num_ricampionamento = 0; num_ricampionamento < num_ricampionamenti; num_ricampionamento++)
    {
        double somma_campione = 0.0;
        int num_blocco;
        for (num_blocco = 0; num_blocco < num_blocchi; num_blocco++)
        {
            int j = rand() % N;
            int extra = j + dim_blocco - N;

            if (extra <= 0)
            {
                double somma = 0.0;
                int l;
                for (l = j; l < j + dim_blocco; l++)
                {
                    somma += lista[l];
                }
                somma_campione += somma;
            }
            else
            {
                double somma1 = 0.0;
                double somma2 = 0.0;
                int l;
                for (l = j; l < N; l++)
                {
                    somma1 += lista[l];
                }
                for (l = 0; l < extra; l++)
                {
                    somma2 += lista[l];
                }
                somma_campione += somma1 + somma2;
            }
        }
        medie_campioni[num_ricampionamento] = somma_campione / N;
    }

    double media_x = 0.0;
    double dx = 0.0;
    int i;
    for (i = 0; i < num_ricampionamenti; i++)
    {
        media_x += medie_campioni[i];
    }
    media_x /= num_ricampionamenti;

    for (i = 0; i < num_ricampionamenti; i++)
    {
        dx += pow((medie_campioni[i] - media_x), 2);
    }
    dx /= (num_ricampionamenti - 1);

    return sqrt(dx);
}
