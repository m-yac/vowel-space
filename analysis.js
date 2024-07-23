
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
var spaceCanvas = document.getElementById("spaceCanvas");
var spaceCtx = spaceCanvas.getContext("2d");

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

function find_peaks(x, epsilon) {
    let peaks = [];
    let i = 0;
    let v1 = x[i+1];
    let d10 = x[i+1] - x[i+0];
    let d21 = x[i+2] - x[i+1];
    let d32 = x[i+3] - x[i+2];
    let v210 = (d21 + d10) / 2;
    let d210 = d21 - d10;
    let d321 = d32 - d21;
    while (i < x.length - 3) {
        if (v1 > epsilon && d10 > 0 && d21 < 0 ||
            v210 < -epsilon && d10 < 0 && d21 < 0 && d210 > 0 && d321 < 0 ||
            v210 >  epsilon && d10 > 0 && d21 > 0 && d210 < 0 && d321 > 0) {
            peaks.push(i+1);
        }
        i += 1;
        v1 = x[i+1];
        d10 = d21;
        d21 = d32;
        d32 = x[i+3] - x[i+2];
        v210 = (d10 + d21) / 2;
        d210 = d321;
        d321 = d32 - d21;
    }
    return peaks;
}

var pulseFFTArray = [];
var sincFFTArray = [];
var outFFTArray = [];
function copyData(toCopy = sincFFTArray) {
    navigator.clipboard.writeText(toCopy.slice(0, M/4 - 4));
}

var Analysis = {
    M: 4096,
    minF1: 10,
    maxF1: 75,
    minF2: 40,
    maxF2: 300,

    init : function()
    {
        const log2MinF1 = Math.log2(this.minF1);
        const log2MaxF1 = Math.log2(this.maxF1);
        const log2MinF2 = Math.log2(this.minF2);
        const log2MaxF2 = Math.log2(this.maxF2);
        this.bF1 = 1 / (log2MaxF1 / log2MinF1 - 1);
        this.bF2 = 1 / (log2MaxF2 / log2MinF2 - 1);
        this.mF1 = this.bF1 / log2MinF1;
        this.mF2 = this.bF2 / log2MinF2;
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
        const iOff = Math.round((this.M / 2) % (1 / fsample));
        for (let i = 0; i < this.M/2; i += 1) {
            const x = 2 * i / this.M * amplCanvas.width
            const y = (0.5 - 1.25 * outArray[this.M/4 + i - iOff]) * amplCanvas.height;
            amplCtx.lineTo(x, y)
        }
        amplCtx.stroke()

        freqCtx.clearRect(0, 0, freqCanvas.width, freqCanvas.height);
        freqCtx.strokeStyle = "blue";
        freqCtx.beginPath()
        for (let i = 0; i < this.M/4; i += 1) {
            const x = Math.log2(1 + 4 * i / this.M) * freqCanvas.width
            const h = Math.log2(Math.max(1e-10, pulseFFTArray[i])) / 18;
            const y = (0.5 - Math.clamp(h, -0.5, 0.5)) * freqCanvas.height;
            freqCtx.lineTo(x, y)
        }
        freqCtx.stroke()
        freqCtx.strokeStyle = "black";
        freqCtx.beginPath()
        for (let i = 0; i < this.M/4; i += 1) {
            const x = Math.log2(1 + 4 * i / this.M) * freqCanvas.width;
            const h = Math.log2(Math.max(1e-10, outFFTArray[i])) / 18;
            const y = (0.5 - Math.clamp(h, -0.5, 0.5)) * freqCanvas.height;
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
        const sincFFTPeaks = find_peaks(sincFFTArray, 0.1);
        for (let j = 0; j < sincFFTPeaks.length; j += 1) {
            const i = sincFFTPeaks[j];
            const x = Math.log2(1 + 4 * i / this.M) * filtCanvas.width
            filtCtx.beginPath()
            filtCtx.moveTo(x, 0)
            filtCtx.lineTo(x, filtCanvas.height)
            filtCtx.stroke()
        }
        filtCtx.stroke()

        document.getElementById("formantDiv").innerHTML = `F1 = ${sincFFTPeaks[0]*sampleRate/this.M}, F2 = ${sincFFTPeaks[1]*sampleRate/this.M}`;

        // spaceCtx.clearRect(0, 0, spaceCanvas.width, spaceCanvas.height);
        
        spaceCtx.beginPath()
        spaceCtx.moveTo(spaceCanvas.width, spaceCanvas.height)
        const xboundary = (1 - this.mF2 * Math.log2(this.maxF1) + this.bF2) * spaceCanvas.width;
        const yboundary = (    this.mF1 * Math.log2(this.minF2) - this.bF1) * spaceCanvas.height;
        spaceCtx.lineTo(spaceCanvas.width, yboundary);
        spaceCtx.lineTo(xboundary, spaceCanvas.height);
        spaceCtx.fill();

        if (sincFFTPeaks.length >= 2) {
            spaceCtx.beginPath();
            const x = (1 - this.mF2 * Math.log2(sincFFTPeaks[1]) + this.bF2) * spaceCanvas.width;
            const y = (    this.mF1 * Math.log2(sincFFTPeaks[0]) - this.bF1) * spaceCanvas.height;
            spaceCtx.arc(x, y, 10, 0, 2*Math.PI);
            spaceCtx.fill();
        }
    }
}

Analysis.init();
