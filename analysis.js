
// Returns 2 ** ceil(log2(x))
function pow2ceil(x) {
    x -= 1;
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);
    return x + 1;
}

// Some common windowing functions

function hann(i, N) {
    const x = Math.sin(Math.PI * i / N);
    return x * x;
}
function cos4Window(a0, a1, a2, a3, i, N) {
    const t = 2 * Math.PI * i / N;
    return a0 - a1 * Math.cos(1 * t)
              + a2 * Math.cos(2 * t)
              - a3 * Math.cos(3 * t);
}
function blackmanHarris(i, N) {
    return cos4Window(0.35875, 0.48829, 0.14128, 0.01168, i, N);
}
function blackmanNutall(i, N) {
    return cos4Window(0.3635819, 0.4891775, 0.1365995, 0.0106411, i, N);
}
function nutall(i, N) {
    return cos4Window(0.355768, 0.487396, 0.144232, 0.012604, i, N);
}

function fft(timeArray) {
    const pow2length = pow2ceil(timeArray.length);
    const lengthAdj = pow2length - timeArray.length;
    if (lengthAdj > 0) {
        timeArray = timeArray.concat(Array(lengthAdj).fill(0.0));
    }
    const f = new FFT(pow2length);
    const complexFreqArray = f.createComplexArray();
    f.realTransform(complexFreqArray, timeArray);
    const freqArray = [];
    for (let j = 0; j < complexFreqArray.length / 2; j += 1) {
        const r = complexFreqArray[2*j];
        const i = complexFreqArray[2*j+1];
        freqArray[j] = Math.sqrt(r*r + i*i);
    }
    return freqArray;
}

function sincFilter(j, M) {
    const x = 0.5 * Math.PI * (j - M/2);
    const n = x == 0 ? 1 : Math.sin(x) / x;
    return n * blackmanNutall(j, M);
}

var freqCanvas = document.getElementById("freqCanvas");
var freqCtx = freqCanvas.getContext("2d");
var amplCanvas = document.getElementById("amplCanvas");
var amplCtx = amplCanvas.getContext("2d");
var filtCanvas = document.getElementById("filtCanvas");
var filtCtx = filtCanvas.getContext("2d");

// A duplicate of Tract.init, but only including variables which change after
// each call to Tract.runStep
function TractModel() {
    this.lipOutput = 0;
    this.noseOutput = 0;
    this.transients = Tract.transients.map((transient) => ({...transient}));

    this.R = new Float64Array(Tract.n);
    this.L = new Float64Array(Tract.n);
    this.junctionOutputR = new Float64Array(Tract.n+1);
    this.junctionOutputL = new Float64Array(Tract.n+1);
    this.maxAmplitude = new Float64Array(Tract.n);
    
    this.noseR = new Float64Array(Tract.noseLength);
    this.noseL = new Float64Array(Tract.noseLength);
    this.noseJunctionOutputR = new Float64Array(Tract.noseLength+1);
    this.noseJunctionOutputL = new Float64Array(Tract.noseLength+1);
    this.noseMaxAmplitude = new Float64Array(Tract.noseLength);

    return this;
}

// A duplicate of the relevant part of AudioSystem.doScriptProcessor
TractModel.prototype.runAnalysisStep = function runAnalysisStep(glottalOutput, j, M) {
    const lambda1 = j/M;
    const lambda2 = (j+0.5)/M;

    let vocalOutput = 0;
    Tract.runStep(glottalOutput, 0, lambda1, this);
    vocalOutput += this.lipOutput + this.noseOutput;
    Tract.runStep(glottalOutput, 0, lambda2, this);
    vocalOutput += this.lipOutput + this.noseOutput;
    return vocalOutput * 0.125;
}

// Port of scypy's _local_maxima_1d (the implementation of find_peaks)
// Source: https://github.com/scipy/scipy/blob/87c46641a8b3b5b47b81de44c07b840468f7ebe7/scipy/signal/_peak_finding_utils.pyx#L20
function find_peaks(x) {
    let maxima = [];
    for (let i = 1; i < x.length - 1; i += 1) {
        // Test if previous sample is smaller
        if (x[i - 1] < x[i] && x[i] > 0.1) {
            let i_ahead = i + 1; // Index to look ahead of current sample

            // Find next sample that is unequal to x[i]
            while (i_ahead < x.length - 1 && Math.abs(x[i_ahead] - x[i]) < 1e-7) {
                i_ahead += 1;
            }

            // Maxima is found if next unequal sample is smaller than x[i]
            if (x[i_ahead] < x[i]) {
                maxima.push(Math.floor((i + i_ahead - 1) / 2));
                // Skip samples that can't be maximum
                i = i_ahead;
            }
        }
    }
    return maxima;
}

var pulseFFTArray = [];
var sincFFTArray = [];
var outFFTArray = [];
function copyData(toCopy = sincFFTArray) {
    navigator.clipboard.writeText(toCopy.slice(0, M/4 - 4));
}

var Analysis = {
    M: 4096,

    init : function()
    {
    },

    draw : function(outArray)
    {
        const fsample = Glottis.smoothFrequency / sampleRate;
        const fblock = this.M * fsample;

        let pulseModel = new TractModel();
        let sincModel = new TractModel();
        let pulseArray = new Array(this.M);
        let sincArray = new Array(this.M);
        for (let j = 0; j < this.M; j += 1) {
            const t = fsample * j;
            const pulse = t > 1 ? 0 : Glottis.normalizedLFWaveform(t);
            const sinc = sincFilter(j, this.M);
            pulseArray[j] = fblock * pulseModel.runAnalysisStep(pulse, j, this.M);
            sincArray[j] = sincModel.runAnalysisStep(sinc, j, this.M);
        }
        pulseFFTArray = fft(pulseArray);
        sincFFTArray = fft(sincArray);
        outFFTArray = fft(outArray.map((x, j) => 2 * x * hann(j, this.M)));

        amplCtx.clearRect(0, 0, amplCanvas.width, amplCanvas.height);
        amplCtx.strokeStyle = "black";
        amplCtx.beginPath()
        for (let i = 0; i < this.M; i += 1) {
          const x = i / this.M * amplCanvas.width
          const y = (0.5 - 0.5 * outArray[i]) * amplCanvas.height;
          amplCtx.lineTo(x, y)
        }
        amplCtx.stroke()

        freqCtx.clearRect(0, 0, freqCanvas.width, freqCanvas.height);
        freqCtx.strokeStyle = "blue";
        freqCtx.beginPath()
        for (let i = 0; i < this.M/4; i += 1) {
          const x = Math.log2(1 + 4 * i / this.M) * freqCanvas.width
          const y = (0.25 - Math.log2(Math.max(1e-10, 0.05 * pulseFFTArray[i])) / 20) * freqCanvas.height;
          freqCtx.lineTo(x, y)
        }
        freqCtx.stroke()
        freqCtx.strokeStyle = "black";
        freqCtx.beginPath()
        for (let i = 0; i < this.M/4; i += 1) {
          const x = Math.log2(1 + 4 * i / this.M) * freqCanvas.width
          const y = (0.25 - Math.log2(Math.max(1e-10, 0.05 * outFFTArray[i])) / 20) * freqCanvas.height;
          freqCtx.lineTo(x, y)
        }
        freqCtx.lineTo(freqCanvas.width, freqCanvas.height);
        freqCtx.lineTo(0, freqCanvas.height);
        freqCtx.fill()

        filtCtx.clearRect(0, 0, filtCanvas.width, filtCanvas.height);
        filtCtx.strokeStyle = "red";
        filtCtx.beginPath()
        const sincFFTMax = Math.max(...sincFFTArray);
        const sincFFTMin = Math.min(...sincFFTArray);
        for (let i = 0; i < this.M/4; i += 1) {
          const x = Math.log2(1 + 4 * i / this.M) * filtCanvas.width
          const y = filtCanvas.height * (1 - (sincFFTArray[i] - sincFFTMin) / (sincFFTMax - sincFFTMin));
          filtCtx.lineTo(x, y)
        }
        filtCtx.stroke()
        filtCtx.strokeStyle = "black";
        const sincFFTPeaks = find_peaks(sincFFTArray);
        for (let j = 0; j < sincFFTPeaks.length; j += 1) {
            const i = sincFFTPeaks[j];
            const x = Math.log2(1 + 4 * i / this.M) * filtCanvas.width
            filtCtx.beginPath()
            filtCtx.moveTo(x, 0)
            filtCtx.lineTo(x, filtCanvas.height)
            filtCtx.stroke()
        }
        filtCtx.stroke()
    }
}

Analysis.init();
