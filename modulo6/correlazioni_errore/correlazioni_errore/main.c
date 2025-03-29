//
//  main.c
//  correlazioni_errore
//
//  Created by Monica  on 01/06/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double calcolaMedia(float* samples, int num_samples) {
    float sum = 0;
    for (int i = 0; i < num_samples; i++) {
        sum += samples[i];
    }
    double mean = (double)sum / num_samples;
    return mean;
}

double calcoloVarianza(float *samples, int num_samples){
    double sum2 = 0.0, var;
    double mean = calcolaMedia(samples,num_samples);
    for (int i = 0;i < num_samples; i++){
        sum2 += pow(samples[i],2);
    }
    sum2 /= num_samples;
    var = num_samples/(num_samples - 1)*(sum2 -  pow(mean, 2));
    return var;
}

double compute_autocorrelation(float* samples, int num_samples, int delay) {
    double mean = 0.0;
    double covariance = 0.0;
    double variance = 0.0;
    
    // Calcola la media e la varianza
    mean = calcolaMedia(samples,num_samples);
    variance = calcoloVarianza(samples, num_samples);
    
    // Calcola la covarianza
    for (int i = 0; i < num_samples - delay; i++) {
        covariance += (samples[i] - mean) * (samples[i + delay] - mean);
    }
    
    covariance /= (num_samples - delay);
    
    // Calcola l'autocorrelazione normalizzata
    return covariance / variance;
}



int countLines(FILE* file) {
    int lines = 0;
    int ch;
    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            lines++;
        }
    }
    rewind(file);
    return lines;
}

int main(void) {
    FILE* inputFile = fopen("/Users/monicacesario/Desktop/modulo3/prove3.dat", "r");
    if (inputFile == NULL) {
        fprintf(stderr, "Impossibile aprire il file di input.\n");
        return 1;
    }

    FILE* outputFile = fopen("/Users/monicacesario/Desktop/modulo3/out_prove3.dat", "w");
    if (outputFile == NULL) {
        fprintf(stderr, "Impossibile aprire il file di output.\n");
        return 1;
    }

    // Salta le prime 5000 righe del file
   // for (int i = 0; i < 5000; i++) {
    //    char line[100];
   //     fgets(line, sizeof(line), inputFile);
   // }

    int num_samples = countLines(inputFile);
    float* samples = malloc(num_samples * sizeof(float));

    // Leggi i campioni dalla seconda colonna del file di input
    for (int i = 0; i < num_samples; i++) {
        float sample1, sample2;
        fscanf(inputFile, "%f %f", &sample1, &sample2);
        samples[i] = sample2;
    }

    fclose(inputFile);
    
    int max_delay = 2000;

    // Calcola l'autocorrelazione e gli errori, scrivendo i risultati su file
    for (int delay = 0; delay <= max_delay; delay++) {
        double autocorrelation = compute_autocorrelation(samples, num_samples, delay);
       
        fprintf(outputFile, "%d %lf %lf\n", delay, autocorrelation);
    }

    fclose(outputFile);

    free(samples);

    return 0;
}
