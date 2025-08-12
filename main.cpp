/*
    This code simulates the sound wave produced by a vibrating string based on physical parameters such as string material, tension, and length. Using harmonic frequency models, it calculates the fundamental frequency and overtones of the string. The resulting waveform shows how changes in string properties affect pitch, timbre, and resonance. However, it is incomplete. Acoustic enhancements like sympathetic string resonance, convolution with body impulse responces, and a detailed soundboard or instrument body vibration model are not yet implemented. These elements are essential for capturing the full richness and realism of an asoustic instrument's tone and will be added in future updates.
*/
#include <cmath>
#include <vector>
#include <random>
#include <iostream>


const double vAir  = 1.789e-5;
const double dAir  = 1.225;
const double twoPI = 2 * M_PI;

struct String {
    double rI, rO, dI, dO, Rho, E, vs, T, L, A, I, mu, angle, Ys, f0, inharmonicity;
    
    void Update() {
        dI = 2 * rI;
        dO = 2 * rO;
        
        A  = M_PI*rI*rI;
        I  = A*rI*rI/4;
        mu = Rho * A;
        Ys = 1. / sqrt(T * mu) ;
        f0 = ((1. / (2 * L)) * sqrt(T/mu));
        
        const double k = (dI / dO);
        inharmonicity = ((M_PI*M_PI*E) / (64*Rho*L*L*L*L*f0*f0)) * (k*k*k*k*dO*dO);
    }
    double getYsn(const double n) {
        return (twoPI * n * f0 * 0.025) / T;
    }
    
    double getInharmonicity() {
        return inharmonicity;
    }
    double getHarmonic(const double n) {
        return f0 * n * sqrt(1 + inharmonicity * n * n);
    }
}; 
struct Bridge {
    double Ewb, vwb, Rho_w, Yb, Lcwb;
    
    void Update() {
        double th = 0.2;
        double D = (Ewb * th*th*th) / (12 * (1 - vwb*vwb));
        Yb = 1. / (8 * sqrt(Rho_w * th * D));
    }
};


struct WavHeader {
    char riff[4] =  { 'R', 'I', 'F', 'F' };
    uint32_t riffSize = 0;
    
    char format[4] = { 'W', 'A', 'V', 'E' };
    
    //char fact[4] = { 'f', 'a', 'c', 't' };
    //uint32_t factSize = 4;
    //char factData[4] = {0};
    
    char ftm[4] = { 'f', 'm', 't', ' ' };    
    uint32_t ftmSize = 16;    
    uint16_t audioFormat;
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
    
    char data[4] = {  'd', 'a', 't', 'a' };
    uint32_t dataSize;
};

enum { FORMAT_PCM = 1, FORMAT_IEEEFLOAT = 3 };

void setHeader(WavHeader& h,
               const int sampleRate,
               const int bitsPerSample,
               const int channels,
               const double durationSeconds,
               const int format) {
    h.audioFormat   = format;
    h.numChannels   = channels;
    h.sampleRate    = sampleRate;
    h.bitsPerSample = bitsPerSample;
    
    h.byteRate      = h.sampleRate * h.numChannels * (h.bitsPerSample / 8);
    h.blockAlign    = h.numChannels * (h.bitsPerSample / 8);
    
    uint32_t dataSize = static_cast<uint32_t>(h.sampleRate * durationSeconds) * h.blockAlign;
    h.dataSize = dataSize;
    h.riffSize = 36 + h.dataSize;
}

std::random_device rd;
std::uniform_int_distribution<int> genValue(0, 0x7fffffff);

double Rand() {
    return ((double)genValue(rd) / RAND_MAX) * 2 * M_PI;   
}

void genWave(std::vector<float>& data, const double amplitude, const double phase, String& s, Bridge& b,
             const double frequency, const double duration, const int sampleRate, const double n) {
    const int numSamples = duration * sampleRate;
    const double intr = (s.E*M_PI*M_PI*s.dI*s.dI*n*n)/(64*s.Rho*s.L*s.L*s.L*s.L*s.f0*s.f0+s.E*M_PI*M_PI*s.dI*s.dI*n*n);
    const double t    = (s.dI / 4) * sqrt((twoPI * frequency) / vAir); 
    const double ext  = (dAir / s.Rho) * (2 * sqrt(2) * t + 1) / (t * t);    
    const double br   = b.Yb / (M_PI * n * (s.Ys + s.getYsn(n)));
    const double decayRate = M_PI * frequency * (intr + ext + br);
    
    std::cout << std::fixed << "decayRate: " << decayRate << "\n";

    for(int i = 0; i < numSamples; i++) {
        double timePos = static_cast<double>(i) / sampleRate;
        double s = tanh(20 * timePos) * sin(twoPI * (frequency) * timePos + phase);        
        data[i] += amplitude * pow(2.718, -decayRate * timePos) * s;
    }
}

void gen(std::vector<float>& data, String& s, Bridge& wb, const double force) {
    double point = 0.975;
    
    const int num = 13;

    
    double f[num];
    double amp[num];
    for(int i = 0; i < num; i++) 
        f[i] = s.getHarmonic(i + 1);
    
    
    for(int i = 0; i < num; i++)
        amp[i] = force*sin(  (i+1) * M_PI * (s.L * point)/s.L)/((i+1) * (i+1));
    for(int i = 0; i < num; i++) 
        genWave(data, amp[i], Rand(), s, wb, f[i],                 20, 96000, i+1);  
    for(int i = 0; i < num; i++) 
        genWave(data, amp[i], Rand(), s, wb, f[i] + 0.0001 * f[i], 20, 96000, i+1);
    for(int i = 0; i < num; i++) 
        genWave(data, amp[i], Rand(), s, wb, f[i] - 0.0001 * f[i], 20, 96000, i+1);
}

int main(int argc, char **argv) {
    String st_466_164;     //B4 flat
    String st_440;         //A4
    String st_415_305;     //A4 flat
    Bridge wb;
    

    st_466_164.rI    = 0.000475;
    st_466_164.rO    = 0.000475;
    st_466_164.vs    = 0.28;
    st_466_164.E     = 2e9;
    st_466_164.Rho   = 7840;
    st_466_164.T     = 636.508;
    st_466_164.L     = 0.363; 
    st_466_164.angle = 2;
    st_466_164.Update();
    
    st_440.rI    = 0.000475;
    st_440.rO    = 0.000475;
    st_440.vs    = 0.28;
    st_440.E     = 2e9;
    st_440.Rho   = 7840;
    st_440.T     = 627.98;
    st_440.L     = 0.382; 
    st_440.angle = 2;
    st_440.Update();

    st_415_305.rI    = 0.000475;
    st_415_305.rO    = 0.000475;
    st_415_305.vs    = 0.28;
    st_415_305.E     = 2e9;
    st_415_305.Rho   = 7840;
    st_415_305.T     = 638.216;
    st_415_305.L     = 0.408; 
    st_415_305.angle = 2;
    st_415_305.Update();


    wb.Ewb   = 1.5e7;
    wb.vwb   = 0.3;
    wb.Rho_w = 700;
    wb.Lcwb  = 0.015;
    wb.Update();
    
    std::vector<float> data(20 * 96000);
    
    gen(data, st_466_164, wb, 0.05);
    gen(data, st_440,     wb, 1.5);
    gen(data, st_415_305, wb, 0.05);
    
    WavHeader wh;
    setHeader(wh, 96000, 32, 1, 20.0, FORMAT_IEEEFLOAT);
    
    FILE* file = fopen("sound_sample_st_440.wav", "wb");
    if(file) {
        fwrite(&wh, sizeof(wh), 1, file);
        fwrite(&data[0], sizeof(float), data.size(), file);
        fclose(file);
    }
    
    
    return 0;
}

