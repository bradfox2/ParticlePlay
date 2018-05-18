//Performs FFT for a frequency value provided in the Particle.io Console with configurable variables
//Small Change
#include <math.h>

//set pins
const int mic = A0;
const int BOARD_LED_PIN = D7;

unsigned int sampling_period_us;
unsigned long microseconds;

// set the sizes for the FFT
const int m = 10; // 11 is as high as photon can handle before a failed flash
const int FFT_SIZE = pow(2, m); //2048 unique samples

// initialize the buffers for the FFT
float samples[FFT_SIZE * 2];
float imagine[FFT_SIZE * 2];
int sampleCounter = 0;

int freq = 10000;
int freq_range_high = freq/2;
double num_bins = FFT_SIZE/2;
int freq_step_per_bin = freq_range_high/num_bins;
int high_freq;

void setup() {
    Particle.variable("MaxNumSampls", FFT_SIZE);
    Particle.variable("HFreqRange", freq_range_high);
    Particle.variable("NumFreqBins", num_bins);
    Particle.variable("FrqRngPerBin", freq_step_per_bin);
    Particle.variable("HighFreqNow", high_freq);
    Particle.function("GetFFTForFrq", find_freq);
    Particle.function("SetHighFreq", set_high_freq);

	// Set up ADC and audio input
	pinMode(mic, INPUT);

	//divide this number by 2 to get the total frequency range we can observe,
	//divide that result by FFT_size to get the frequency range per 'bin',
	//number of bins will be equal to FFT_size
	sampling_period_us = round(1000000*(1.0/freq));
}

void loop() {
    run_fft(1);
    high_freq = 1;
    for(int i = 2; i < sizeof(samples)/sizeof(int); i++){
        if(samples[i] > samples[i-1]){
            high_freq  = samples[i] * (freq_step_per_bin/2);
        }
    }
    Serial.println(high_freq);
}


int run_fft(int bin){
	// if we're currently filling the buffer
    for(int i=0; i<(FFT_SIZE * 2); i++)
        {
	    microseconds = micros();

		samples[i] = analogRead(mic);
		// the buffer full of imaginary numbers for the FFT doesn't matter for our use, so fill it with 0s
		imagine[i] = 0.0;

		//delay so that the samples have even timing
		while(micros() < (microseconds + sampling_period_us)){
        }
	}

	// run FFT
 	FFT(1, m, samples, imagine);
 	return(samples[bin]);
}

//Particle linked function,  find the freq bin for given freq, return the fft value
int find_freq(String requested_freq){
    int bin_p = round(atoi(requested_freq)/freq_step_per_bin);
    if(atoi(requested_freq) > freq_range_high){
        return(-1);
    }
    return(run_fft(bin_p));
}

int set_high_freq(String new_high_freq){
    freq = atoi(new_high_freq)*2;
    freq_range_high = freq/2;
    freq_step_per_bin = freq_range_high/num_bins;
    return(1);
}

//Perform FFT
short FFT(short int dir, int m, float *rx, float *iy) {

	/*
	FFT() from Paul Bourke: http://paulbourke.net/miscellaneous/dft/
	as referenced by @phec on the Particle forums: https://community.particle.io/t/fast-fourier-transform-library/3784/4
	This computes an in-place complex-to-complex FFT
	rx and iy are the real and imaginary arrays of 2^m points.

	dir gives FFT_FORWARD or FFT_REVERSE transform
	rx is the array of real numbers on input, and the x coordinates on output
	iy is the array of imaginary numbers on input, and the y coordinates on output
	*/

	// \/ \/ \/ DO NOT EDIT THIS CODE UNLESS YOU KNOW WHAT YOU'RE DOING \/ \/ \/

	int n, i, i1, j, k, i2, l, l1, l2;
	float c1, c2, tx, ty, t1, t2, u1, u2, z;

	/* Calculate the number of points */
	n = 1;
	for (i=0;i<m;i++)
		n *= 2;

	/* Do the bit reversal */
	i2 = n >> 1;
	j = 0;
	for (i=0;i<n-1;i++) {
		if (i < j) {
			tx = rx[i];
			ty = iy[i];
			rx[i] = rx[j];
			iy[i] = iy[j];
			rx[j] = tx;
			iy[j] = ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l=0;l<m;l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for (j=0;j<l1;j++) {
			for (i=j;i<n;i+=l2) {
				i1 = i + l1;
				t1 = u1 * rx[i1] - u2 * iy[i1];
				t2 = u1 * iy[i1] + u2 * rx[i1];
				rx[i1] = rx[i] - t1;
				iy[i1] = iy[i] - t2;
				rx[i] += t1;
				iy[i] += t2;
			}
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == 1)
			c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);
	}

	/* Scaling for forward transform */
	if (dir == 1) {
		for (i=0;i<n;i++) {
			rx[i] = abs(rx[i]/n);
			iy[i] = abs(iy[i]/n);
		}
	}

	return(0);
}
