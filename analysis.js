
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

// Inserts an element into a sorted array
function insertSorted(arr, x, compareFn) {
    let i = 0;
    while (i < arr.length && compareFn(arr[i], x) < 0) { i++; }
    arr.splice(i, 0, x);
}

function det2(m00, m01,
              m10, m11) {
  return m00 * m11 - m01 * m10;
}

function det3(m00, m01, m02,
              m10, m11, m12,
              m20, m21, m22) {
    return   m00 * det2(m11, m12,
                        m21, m22)
           - m01 * det2(m10, m12,
                        m20, m22)
           + m02 * det2(m10, m11,
                        m20, m21);
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

// Hilbert curve code
// Source: https://observablehq.com/@mourner/3d-hilbert-curves

// extract every 3rd bit from a number to decode a component of a Morton index
compact1By2 = (x) => { // https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
    x &= 0x09249249;                 
    x = (x ^ (x >> 2)) & 0x030c30c3;
    x = (x ^ (x >> 4)) & 0x0300f00f;
    x = (x ^ (x >> 8)) & 0xff0000ff;
    x = (x ^ (x >> 16)) & 0x000003ff;
    return x;
}

// lookup table for converting a Hilbert value to a Morton index
const hilbertToMortonTable = Uint8Array.of(
    48, 33, 35, 26, 30, 79, 77, 44,
    78, 68, 64, 50, 51, 25, 29, 63,
    27, 87, 86, 74, 72, 52, 53, 89,
    83, 18, 16,  1,  5, 60, 62, 15,
     0, 52, 53, 57, 59, 87, 86, 66,
    61, 95, 91, 81, 80,  2,  6, 76,
    32,  2,  6, 12, 13, 95, 91, 17,
    93, 41, 40, 36, 38, 10, 11, 31,
    14, 79, 77, 92, 88, 33, 35, 82,
    70, 10, 11, 23, 21, 41, 40,  4,
    19, 25, 29, 47, 46, 68, 64, 34,
    45, 60, 62, 71, 67, 18, 16, 49
);

// convert a Hilbert value into a Morton index using a lookup table
function transformCurve(index, bits, lookupTable) {
    let out = 0;
    for (let i = 3 * (bits - 1), transform = 0; i >= 0; i -= 3) {
      transform = lookupTable[transform | ((index >> i) & 7)];
      out = (out << 3) | (transform & 7);
      transform &= ~7;
    }
    return out;
};

// ----------------------------------------------------------------

var rawFormantData;
var boundaryComponents;
var minF1; var minF2; var minF3; var minF4;
var maxF1; var maxF2; var maxF3; var maxF4;

var oppOrchid = '#968ffa'; // OKLAB = 96.45% 0.02 284.69
var oppPalePink = '#e9e9ff';

var allDiv = document.getElementById("allDiv");
var containersDiv = document.getElementById("containersDiv");
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
    this.transients = []; //Tract.transients.map((transient) => ({...transient}));

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

function find_peaks(y) {
    let peaks = [];
    let x = 0;
    let dy_1 = y[x+2] - y[x+0];
    let dy_2 = y[x+3] - y[x+1];
    let dy_3 = y[x+4] - y[x+2];
    let dy_4 = y[x+5] - y[x+3];
    let dy_5 = y[x+6] - y[x+4];
    let ddy_2 = dy_3 - dy_1;
    let ddy_3 = dy_4 - dy_2;
    let ddy_4 = dy_5 - dy_3;
    let x_start = 0;
    let dy_start = 0;
    while (x < y.length - 6) {
        if (ddy_2 >= 0 && ddy_3 < 0) {
            const t = - ddy_2 / (ddy_3 - ddy_2);
            x_start = x + 2 + t;
            dy_start = dy_2 + (dy_3 - dy_2) * t;
        }
        else if (ddy_2 < 0 && ddy_3 >= 0) {
            const t = - ddy_2 / (ddy_3 - ddy_2);
            const x_end = x + 2 + t;
            const dy_end = dy_2 + (dy_3 - dy_2) * t;
            if (dy_start > dy_end) {
                peaks.push([(x_start + x_end) / 2, x_end - x_start, dy_start - dy_end]);
            }
        }
        x += 1;
        dy_1 = dy_2;
        dy_2 = dy_3;
        dy_3 = dy_4;
        dy_4 = dy_5;
        dy_5 = y[x+6] - y[x+4];
        ddy_2 = ddy_3;
        ddy_3 = ddy_4;
        ddy_4 = dy_5 - dy_3;
    }
    return peaks;
}

// var outArray = [];
// var middlePeakOffset = 0;
// var pulseFFTArray = [];
// var sincFFTArray = [];
// var sincFFTPeaks = [];
// var outFFTArray = [];
// var normSincFFTArray = [];

function copyData(toCopy = sincFFTArray) {
    navigator.clipboard.writeText(toCopy.slice(0, Analysis.M/4 - 4));
}

function copyPeakData() {
    let obj = { peaks: sincFFTPeaks, spectrum: normSincFFTArray };
    navigator.clipboard.writeText(JSON.stringify(obj));
}

function saveGlobalVars() {
    const js_strs = [`var outArray = ${JSON.stringify(outArray)};`,
                     `var middlePeakOffset = ${JSON.stringify(middlePeakOffset)};`,
                     `var pulseFFTArray = ${JSON.stringify(pulseFFTArray)};`,
                     `var sincFFTArray = ${JSON.stringify(sincFFTArray)};`,
                     `var sincFFTPeaks = ${JSON.stringify(sincFFTPeaks)};`,
                     `var outFFTArray = ${JSON.stringify(outFFTArray)};`,
                     `var normSincFFTArray = ${JSON.stringify(normSincFFTArray)};`]
    const jsBlob = new Blob([js_strs.join("\n")], {type: "text/javascript"});
    saveAs(jsBlob, `globalVars.js`);
}

function clickIndexDiameter(index, diameter) {
    if (isNaN(index) || isNaN(diameter)) { return; }
    if (index == TractUI.lipIndex) {
        TractUI.lipDiameter = diameter;
    }
    else {
        TractUI.tongueIndex = index;
        TractUI.tongueDiameter = diameter;
    }
    TractUI.setRestDiameter();   
    for (var i=0; i<Tract.n; i++) Tract.targetDiameter[i] = Tract.restDiameter[i];
}

// Formant data coordinate systems

function toFormantDataIndex(iLipDiam, iTongueY, iTongueX) {
    return (iLipDiam << 10) | (iTongueY << 6) | (iTongueX << 2);
}
function fromFormantDataIndex(i) {
    return [i >> 10, (i >> 6) & 0xF, (i >> 2) & 0xF];
}

function formantData(iLipDiam, iTongueY, iTongueX, iFormant) {
    const i = toFormantDataIndex(iLipDiam, iTongueY, iTongueX);
    if (iFormant !== undefined) {
        return rawFormantData[i | iFormant];
    }
    return [rawFormantData[i|0], rawFormantData[i|1], rawFormantData[i|2], rawFormantData[i|3]];
}

// Normalized lip/tongue values in the range [0,1] (or [-1,1] for normTongueIndex)
function toNormLipTongueValues(iLipDiam, iTongueY, iTongueX) {
    const normLipDiam = iLipDiam / 0xF;
    const iTongueDiam = Math.max(iTongueY, iTongueX);
    const normTongueDiam = iTongueDiam / 0xF;
    const normTongueIndex = iTongueDiam > 0 ? (iTongueY - iTongueX) / iTongueDiam : 0;
    return [normLipDiam, normTongueDiam, normTongueIndex];
}
function fromNormLipTongueValues(normLipDiam, normTongueDiam, normTongueIndex) {
    const iLipDiam = normLipDiam * 15.0;
    const iTongueDiam = normTongueDiam * 15.0;
    const iTongueY = iTongueDiam * (1 + Math.min(normTongueIndex, 0));
    const iTongueX = iTongueDiam * (1 - Math.max(normTongueIndex, 0));
    return [iLipDiam, iTongueY, iTongueX];
}

var c = 0.625;

// Actual lip/tongue values
function toLipTongueValues(iLipDiam, iTongueY, iTongueX) {
    const [normLipDiam, normTongueDiam, normTongueIndex] = toNormLipTongueValues(iLipDiam, iTongueY, iTongueX);
    // Computing lip diameter value in the range [TractUI.innerLipControlRadius, TractUI.outerLipControlRadius]
    const lipDiameter = (TractUI.innerLipControlRadius - c) * Math.pow((TractUI.outerLipControlRadius - c) / (TractUI.innerLipControlRadius - c), normLipDiam) + c;
    // Computing tongue diameter value in the range [TractUI.outerTongueControlRadius, TractUI.innerTongueControlRadius]
    const tongueDiameter = TractUI.outerTongueControlRadius * Math.pow(TractUI.innerTongueControlRadius / TractUI.outerTongueControlRadius, normTongueDiam);
    // Computing tongue index value in the range [TractUI.tongueIndexCentre - out, TractUI.tongueIndexCentre + out], where `out` is defined in terms of the following code from TractUI.handleTouches
    let fromPoint = (TractUI.outerTongueControlRadius-tongueDiameter)/(TractUI.outerTongueControlRadius-TractUI.innerTongueControlRadius);
    fromPoint = Math.pow(fromPoint, 0.58) - 0.2*(fromPoint*fromPoint-fromPoint); //horrible kludge to fit curve to straight line
    const out = fromPoint*0.5*(TractUI.tongueUpperIndexBound-TractUI.tongueLowerIndexBound);
    const tongueIndex = TractUI.tongueIndexCentre + out * normTongueIndex;
    return [lipDiameter, tongueDiameter, tongueIndex];
}
function fromLipTongueValues(lipDiameter, tongueDiameter, tongueIndex) {
    const normLipDiam = Math.clamp(Math.log((lipDiameter - c) / (TractUI.innerLipControlRadius - c)) / Math.log((TractUI.outerLipControlRadius - c) / (TractUI.innerLipControlRadius - c)), 0.0, 1.0);
    const normTongueDiam = Math.clamp(Math.log(tongueDiameter / TractUI.outerTongueControlRadius) / Math.log(TractUI.innerTongueControlRadius / TractUI.outerTongueControlRadius), 0.0, 1.0);
    // As above, the code for computing `out` is from TractUI.handleTouches
    let fromPoint = (TractUI.outerTongueControlRadius-tongueDiameter)/(TractUI.outerTongueControlRadius-TractUI.innerTongueControlRadius);
    fromPoint = Math.pow(fromPoint, 0.58) - 0.2*(fromPoint*fromPoint-fromPoint); //horrible kludge to fit curve to straight line
    const out = fromPoint*0.5*(TractUI.tongueUpperIndexBound-TractUI.tongueLowerIndexBound);
    const normTongueIndex = out > 0 ? Math.clamp((tongueIndex - TractUI.tongueIndexCentre) / out, -1.0, 1.0) : 0;
    return fromNormLipTongueValues(normLipDiam, normTongueDiam, normTongueIndex);
}

function normFormant(f, i) {
    return i != 2 ? f : Math.log2(f);
}
function invNormFormant(f, i) {
    return i != 2 ? f : Math.pow(2, f);
}

const nPadLo = [0.14, 0.10, 0.10];
const nPadHi = [0.05, 0.06, 0.15];
const nMin = [normFormant(minF1, 1),
              normFormant(minF2, 2),
              normFormant(minF3, 3)];
const nRng = [normFormant(maxF1, 1) - nMin[1-1],
              normFormant(maxF2, 2) - nMin[2-1],
              normFormant(maxF3, 3) - nMin[3-1]];
const nSlp = Array.from([0,1,2], (i) => (1 - nPadHi[i] - nPadLo[i]) / nRng[i]);
const nInt = Array.from([0,1,2], (i) => nPadLo[i] - nSlp[i] * nMin[i]);

function toNormalizedFormant(fi, i) {
    return nSlp[i-1] * normFormant(fi, i) + nInt[i-1];
}
function fromNormalizedFormant(nFi, i) {
    return invNormFormant((nFi - nInt[i-1]) / nSlp[i-1], i);
}

function toNormalizedFormantVec(f) {
    return [toNormalizedFormant(f[0],0),
            toNormalizedFormant(f[1],1),
            toNormalizedFormant(f[2],2)];
}

// ...

function FormantTet(corner_index, x, detC) {
    this.corner_index = corner_index;
    this.x = x;
    this.detC = detC;
    this.uid = ((this.x[0][0] + this.x[1][0] + this.x[2][0] + this.x[3][0]) << 12) |
               ((this.x[0][1] + this.x[1][1] + this.x[2][1] + this.x[3][1]) <<  6) |
               ((this.x[0][2] + this.x[1][2] + this.x[2][2] + this.x[3][2]) <<  0);
}

FormantTet.prototype.f = function f(i, j) {
    if (j === undefined) {
        return [formantData(this.x[i][0], this.x[i][1], this.x[i][2], 0),
                formantData(this.x[i][0], this.x[i][1], this.x[i][2], 1),
                formantData(this.x[i][0], this.x[i][1], this.x[i][2], 2), 1.0];
    }
    return j == 3 ? 1.0 : formantData(this.x[i][0], this.x[i][1], this.x[i][2], j);
}

FormantTet.prototype.compute_detC = function compute_detC(i, j) {
    const i0 = (i+1)%4; const i1 = (i+2)%4; const i2 = (i+3)%4;
    const j0 = (j+1)%4; const j1 = (j+2)%4; const j2 = (j+3)%4;
    const detC = det3(this.f(i0,j0), this.f(i0,j1), this.f(i0,j2),
                      this.f(i1,j0), this.f(i1,j1), this.f(i1,j2),
                      this.f(i2,j0), this.f(i2,j1), this.f(i2,j2));
    return i % 2 == 0 ? detC : -detC;
}

FormantTet.prototype.compute_det = function compute_det(i, row) {
    if (!this.detC[i]) {
        this.detC[i] = [this.compute_detC(i,0), this.compute_detC(i,1),
                        this.compute_detC(i,2), this.compute_detC(i,3)];
    }
    return   row[0] * this.detC[i][0] - row[1] * this.detC[i][1]
           + row[2] * this.detC[i][2] - row[3] * this.detC[i][3];
}

FormantTet.prototype.barycentric_coord = function barycentric_coord(f1, f2, f3, i) {
    if (!this.det) {
        this.det = this.compute_det(0, this.f(0));
    }
    return this.compute_det(i, [f1, f2, f3, 1]) / this.det;
}

FormantTet.prototype.barycentric_coords = function barycentric_coords(f1, f2, f3, params) {
    if (params) {
        const [m, b] = params;
        return [m[0] * f3 + b[0],
                m[1] * f3 + b[1],
                m[2] * f3 + b[2],
                m[3] * f3 + b[3]]
    }
    return [this.barycentric_coord(f1, f2, f3, 0),
            this.barycentric_coord(f1, f2, f3, 1),
            this.barycentric_coord(f1, f2, f3, 2),
            this.barycentric_coord(f1, f2, f3, 3)];
}

FormantTet.prototype.contains = function contains(f1, f2, f3, params) {
    return this.barycentric_coords(f1, f2, f3, params).every((li) => li >= -1e-7);
}

FormantTet.prototype.interp = function interp(l) {
    l[0] = Math.max(0, l[0]);
    l[1] = Math.max(0, l[1]);
    l[2] = Math.max(0, l[2]);
    l[3] = Math.max(0, l[3]);
    const sum = l[0] + l[1] + l[2] + l[3];
    l[0] = l[0] / sum;
    l[1] = l[1] / sum;
    l[2] = l[2] / sum;
    l[3] = l[3] / sum;
    return [this.x[0][0] * l[0] + this.x[1][0] * l[1] + this.x[2][0] * l[2] + this.x[3][0] * l[3],
            this.x[0][1] * l[0] + this.x[1][1] * l[1] + this.x[2][1] * l[2] + this.x[3][1] * l[3],
            this.x[0][2] * l[0] + this.x[1][2] * l[1] + this.x[2][2] * l[2] + this.x[3][2] * l[3]]
}

FormantTet.prototype.interp_coords = function interp_coords(f1, f2, f3) {
    return this.interp(this.barycentric_coords(f1, f2, f3));
}

// t.compute_det(i,[f1,f2,f3,1]) = t.det_f3_slope(i) * f3 + t.det_f3_intercept(i,f1,f2)
FormantTet.prototype.det_f3_slope = function det_f3_slope(i) {
    return this.compute_det(i, [0, 0, 1, 0]);
}
FormantTet.prototype.det_f3_intercept = function det_f3_intercept(i, f1, f2) {
    return this.compute_det(i, [f1, f2, 0, 1]);
}
FormantTet.prototype.det_f3_params = function det_f3_params(f1, f2) {
    const m0 = this.det_f3_slope(0);
    const m1 = this.det_f3_slope(1);
    const m2 = this.det_f3_slope(2);
    const m3 = this.det_f3_slope(3);
    const b0 = this.det_f3_intercept(0, f1, f2);
    const b1 = this.det_f3_intercept(1, f1, f2);
    const b2 = this.det_f3_intercept(2, f1, f2);
    const b3 = this.det_f3_intercept(3, f1, f2);
    return [[m0, m1, m2, m3], [b0, b1, b2, b3]];
}

FormantTet.prototype.closest_f3 = function closest_f3(f1, f2, params) {
    // li = this.barycentric_coord(f1, f2, f3, i)
    // li = (mi * f3 + bi) / this.det
    const [m, b] = params ? params : this.det_f3_params(f1, f2);
    // d/df3 li^2 = 2 * (mi^2 * f3 + mi * bi) / this.det
    // d/df3 (l0^2 + l1^2 + l2^2 + l3^2) = 2 * (M * f3 + B) / this.det
    const B = m[0] * b[0] + m[1] * b[1] + m[2] * b[2] + m[3] * b[3];
    const M = m[0] * m[0] + m[1] * m[1] + m[2] * m[2] + m[3] * m[3];
    // d/df3 (l0^2 + l1^2 + l2^2 + l3^2) = 0  ==>  f3 = -B/M
    return -B / M;
}

FormantTet.prototype.intersect_face_f3 = function intersect_face_f3(f1, f2, i, params) {
    const [m, b] = params ? params : this.det_f3_params(f1, f2);
    const f3 = -m[i] / b[i]
    if (this.contains(f1, f2, f3, params)) {
        return f3;
    }
    return null;
}


FormantTet.prototype.intersect_f3 = function intersect_f3(f1, f2, params) {
    if (!params) { params = this.det_f3_params(f1, f2); }
    const intersections = [];
    for (let i = 0; i < 4; i++) {
        const intersection = this.intersect_face_f3(f1, f2, i, params);
        if (intersection != null) {
            intersections.push([intersection, i]);
        }
    }
    return intersections;
}

FormantTet.prototype.normal = function normal(i) {
    const f0 = toNormalizedFormantVec(this.f((i+0)%4));
    const f1 = toNormalizedFormantVec(this.f((i+1)%4));
    const f2 = toNormalizedFormantVec(this.f((i+2)%4));
    const f3 = toNormalizedFormantVec(this.f((i+3)%4));
    const d = [[f0[0]-f3[0], f0[1]-f3[1], f0[2]-f3[2]],
               [f1[0]-f3[0], f1[1]-f3[1], f1[2]-f3[2]],
               [f2[0]-f3[0], f2[1]-f3[1], f2[2]-f3[2]]];
    const c = [det3(      1,       0,       0,
                    d[1][0], d[1][1], d[1][2],
                    d[2][0], d[2][1], d[2][2]),
               det3(      0,       1,       0,
                    d[1][0], d[1][1], d[1][2],
                    d[2][0], d[2][1], d[2][2]),
               det3(      0,       0,       1,
                    d[1][0], d[1][1], d[1][2],
                    d[2][0], d[2][1], d[2][2])];
    const b = c[0] * d[0][0] + c[1] * d[0][1] + c[2] * d[0][2];
    const m = c[0] * c[0]    + c[1] * c[1]    + c[2] * c[2];
    const k = -Math.sign(b) / Math.sqrt(m);
    return [k * c[0], k * c[1], k * c[2]];
}

FormantTet.prototype.neighbor = function neighbor(i) {
    let new_corner_index = this.corner_index;
    let new_x = this.x.slice(); // shallow copy
    if (this.corner_index == i) {
        new_corner_index = -1;
        new_x[i] = [this.x[(i+1)%4][0] + this.x[(i+2)%4][0] + this.x[(i+3)%4][0] - 2 * this.x[i][0],
                    this.x[(i+1)%4][1] + this.x[(i+2)%4][1] + this.x[(i+3)%4][1] - 2 * this.x[i][1],
                    this.x[(i+1)%4][2] + this.x[(i+2)%4][2] + this.x[(i+3)%4][2] - 2 * this.x[i][2], 1.0];
    }
    else if (this.corner_index == -1) {
        new_corner_index = i;
        new_x[i] = [(this.x[(i+1)%4][0] + this.x[(i+2)%4][0] + this.x[(i+3)%4][0] - this.x[i][0]) / 2,
                    (this.x[(i+1)%4][1] + this.x[(i+2)%4][1] + this.x[(i+3)%4][1] - this.x[i][1]) / 2,
                    (this.x[(i+1)%4][2] + this.x[(i+2)%4][2] + this.x[(i+3)%4][2] - this.x[i][2]) / 2, 1.0];
    }
    else {
        new_x[i] = [2 * this.x[this.corner_index][0] - this.x[i][0],
                    2 * this.x[this.corner_index][1] - this.x[i][1],
                    2 * this.x[this.corner_index][2] - this.x[i][2], 1.0];
        if (new_x[i][0] < 0 || new_x[i][0] > 0xF ||
            new_x[i][1] < 0 || new_x[i][1] > 0xF ||
            new_x[i][2] < 0 || new_x[i][2] > 0xF) {
            return null;
        }
    }
    let new_detC = [null, null, null, null];
    new_detC[i] = this.detC[i];
    return new FormantTet(new_corner_index, new_x, new_detC);
}

FormantTet.prototype.nearest_wrt = function nearest_wrt(coords_fn, seen) {
    if (!seen) { seen = new Set(); }
    if (seen.has(this.uid)) { return null; }
    seen.add(this.uid);

    let coords = coords_fn(this);
    coords = coords.map((x,i) => [x,i]);
    coords.sort((a,b) => a[0] - b[0]);
    if (coords[0][0] < 0) {
        const neighbor_i = this.neighbor(coords[0][1]);
        if (neighbor_i) {
            const nearest = neighbor_i.nearest_wrt(coords_fn, seen);
            if (nearest) {
                return nearest;                
            }
        }
    }
    return this;
}

FormantTet.prototype.nearest = function nearest(f1, f2, f3) {
    return this.nearest_wrt((tet) => tet.barycentric_coords(f1, f2, f3));
}

FormantTet.prototype.nearest_f1_f2 = function nearest_f1_f2(f1, f2, f3) {
    return this.nearest_wrt(function (tet) {
        const params = tet.det_f3_params(f1, f2);
        const intersections = tet.intersect_f3(f1, f2, params);
        if (tet.contains(f1, f2, f3, params)) {
            return [1/4, 1/4, 1/4, 1/4];
        }
        else if (intersections.length > 0) {
            intersections.sort((a,b) => Math.abs(a[0] - f3) - Math.abs(b[0] - f3));
            let coords = [2/3, 2/3, 2/3, 2/3];
            coords[intersections[0][1]] = -1;
            return coords;
        }
        else {
            const cls_f3 = tet.closest_f3(f1, f2, params);
            return tet.barycentric_coords(f1, f2, cls_f3, params);
        }
    });
}

FormantTet.treePoint = function treePoint(f1, f2, f3, iLipDiam, iTongueX, iTongueY) {
    const p = [toNormalizedFormant(f1, 1),
               toNormalizedFormant(f2, 2),
               toNormalizedFormant(f3, 3)];
    if (iLipDiam !== undefined) {
        p.index = [iLipDiam, iTongueX, iTongueY];
    }
    return p;
}

FormantTet.treePointFromIndex = function treePointFromIndex(i) {
    return FormantTet.treePoint(rawFormantData[i|0],
                                rawFormantData[i|1],
                                rawFormantData[i|2],
                                (i >> 10), (i >> 6) & 0xF, (i >> 2) & 0xF);
}

FormantTet.treeDistance = function treeDistance(a, b) {
    return Math.pow(a[0] - b[0], 2) + Math.pow(a[1] - b[1], 2) + Math.pow(a[2] - b[2], 2);
}

FormantTet.treeDistance_f1_f2 = function treeDistance_f1_f2(a, b) {
    return Math.pow(a[0] - b[0], 2) + Math.pow(a[1] - b[1], 2);
}

FormantTet.tree = new kdTree(Array.from({ length: 0x1000 }, ((_, i) => FormantTet.treePointFromIndex(i << 2))),
                             FormantTet.treeDistance, [0,1,2]);

FormantTet.boundaryTrees = boundaryComponents.map((a) => a.map(function (b) {
        return new kdTree(Array.from(b, FormantTet.treePointFromIndex),
                          FormantTet.treeDistance_f1_f2, [0,1]);
    }));

FormantTet.withVertexNearestTo = function withVertexNearestTo(p, index, nearest) {
    if (!nearest) {
        const indexRnd = [Math.round(index[0]), Math.round(index[1]), Math.round(index[2])];
        nearest = formantData(indexRnd[0], indexRnd[1], indexRnd[2]);
    }
    const x = [[index[0], index[1], index[2], 1.0],
               [index[0], index[1], index[2], 1.0],
               [index[0], index[1], index[2], 1.0],
               [index[0], index[1], index[2], 1.0]];
    x[0][0] += p[0] > nearest[0] && x[0][0] < 0xF || x[0][0] == 0 ? 1 : -1;
    x[1][1] += p[1] > nearest[1] && x[1][1] < 0xF || x[1][1] == 0 ? 1 : -1;
    x[2][2] += p[2] > nearest[2] && x[2][2] < 0xF || x[2][2] == 0 ? 1 : -1;
    return new FormantTet(3, x, [null, null, null, null]);
}

FormantTet.withNearestVertexTo = function withNearestVertexTo(f1, f2, f3) {
    const p = FormantTet.treePoint(f1, f2, f3);
    const nearest = FormantTet.tree.nearest(p, 1)[0][0];
    return this.withVertexNearestTo(p, nearest.index, nearest);
}

FormantTet.boundaryIntersections = function boundaryIntersections(p) {
    let nearests = [];
    for (const [i ,j] of [[0,0], [2,1], [0,1], [2,0]]) {
        const nearest = FormantTet.boundaryTrees[i][j].nearest(p, 1)[0];
        if (nearest[1] < 0.1) {
            nearests.push([nearest[0], Math.abs(nearest[2] - p[2]), i, j]);
        }
    }
    if (nearests.length == 0) {
        nearests = [[FormantTet.boundaryTrees[1][0].nearest(p, 1)[0][0], 0, 1, 0]];
    }
    return nearests;
}

FormantTet.withNearestBoundaryVertexTo = function withNearestBoundaryVertexTo(f1, f2, f3) {
    const p = FormantTet.treePoint(f1, f2, f3);
    const nearests = FormantTet.boundaryIntersections(p);
    nearests.sort((a,b) => a[1] - b[1]);
    return this.withVertexNearestTo(p, nearests[0][0].index, nearests[0][0]);
}

FormantTet.nearest = function nearest(f1, f2, f3) {
    return FormantTet.withNearestVertexTo(f1, f2, f3).nearest(f1, f2, f3);
}

FormantTet.fromPointAndIndex = function fromPointAndIndex(pt, isCorner, i) {
    const x = [[pt[0], pt[1], pt[2], 1],
               [pt[0], pt[1], pt[2], 1],
               [pt[0], pt[1], pt[2], 1],
               [pt[0], pt[1], pt[2], 1]];
    for (const j of (isCorner ? [0] : [1,2])) {
        const j0 = (j+0)%3; const j1 = (j+1)%3; const j2 = (j+2)%3;
        x[0][j0] += (i >> j0) & 1 ? 1 : -1;
        x[1][j1] += (i >> j1) & 1 ? 1 : -1;
        x[2][j2] += (i >> j2) & 1 ? 1 : -1;
        if (x[0][j0] < 0 || x[0][j0] > 0xF) { return null; }
        if (x[1][j1] < 0 || x[1][j1] > 0xF) { return null; }
        if (x[2][j2] < 0 || x[2][j2] > 0xF) { return null; }
    }
    return new FormantTet(isCorner ? 3 : -1, x, [null, null, null, null]);
}

FormantTet.forEach = function forEach(fn) {
    for (let iLipDiam = 0; iLipDiam <= 0xF; iLipDiam++) {
        for (let iTongueY = 0; iTongueY <= 0xF; iTongueY++) {
            for (let iTongueX = 0; iTongueX <= 0xF; iTongueX++) {
                const isCorner = (iLipDiam + iTongueY + iTongueX) % 2 == 0;
                const is2x2Center = iLipDiam % 2 == 1 && iTongueY % 2 == 1 && iTongueX % 2 == 1;
                if (!isCorner && !is2x2Center) { continue; }
                for (let i = 0; i < 8; i++) {
                    const tet = FormantTet.fromPointAndIndex([iLipDiam, iTongueY, iTongueX], isCorner, i);
                    if (tet) { fn(tet); }
                }
            }
        }
    }
}

function FormantDataCollection() {
    this.wait = -1;
}

FormantDataCollection.prototype.start = function start() {
    if (this.wait < 0) {
        if (!rawFormantData) { rawFormantData = new Array(0x4000); }
        this.iHilbert = 0;
        this.update();
        this.wait = 3;
    }
}

FormantDataCollection.prototype.update = function update() {
    const iMorton = transformCurve(this.iHilbert, 4, hilbertToMortonTable);
    const iLipDiam = compact1By2(iMorton >> 0);
    const iTongueY = compact1By2(iMorton >> 1);
    const iTongueX = compact1By2(iMorton >> 2);
    this.i = toFormantDataIndex(iLipDiam, iTongueX, iTongueY);
    const [lipDiameter, tongueDiameter, tongueIndex] = toLipTongueValues(iLipDiam, iTongueY, iTongueX);
    console.log(this.iHilbert, iMorton, [iLipDiam, iTongueX, iTongueY], [lipDiameter, tongueDiameter, tongueIndex]);
    clickIndexDiameter(TractUI.lipIndex, lipDiameter);
    clickIndexDiameter(tongueIndex, tongueDiameter);
}

FormantDataCollection.prototype.step = function step(peaks) {
    if (this.wait != 0) {
        if (this.wait > 0) { this.wait -= 1; }
        return;
    }
    rawFormantData[this.i + 0] = peaks[0];
    rawFormantData[this.i + 1] = peaks[1];
    rawFormantData[this.i + 2] = peaks[2];
    rawFormantData[this.i + 3] = peaks[3];
    if (this.iHilbert < 0x1000) {
        this.iHilbert += 1;
        this.update();
        this.wait = 1;
    }
    else {
        this.saveAs("formantData", "rawFormantData", Array.from(rawFormantData), []);
        this.find_boundary();
        this.wait = -1;
    }
}

FormantDataCollection.prototype.saveAs = function (filename, varname, data, addl_vars) {
    const json_str = JSON.stringify(data);
    let js_str = `var ${varname} = ${json_str};`
    for (const [addl_varname, addl_data] of addl_vars) {
        js_str += `\nvar ${addl_varname} = ${JSON.stringify(addl_data)};`;
    }
    const jsonBlob = new Blob([json_str], {type: "application/json"});
    const jsBlob = new Blob([js_str], {type: "text/javascript"});
    saveAs(jsonBlob, `${filename}.json`);
    saveAs(jsBlob, `${filename}.js`);
}

FormantDataCollection.prototype.find_boundary = function find_boundary() {
    let boundary_tet_map = new Map();
    let boundary_tet_points_map = new Map();
    FormantTet.forEach(function (tet) {
        for (let i = 0; i < 4; i++) {
            const neighbor_i = tet.neighbor(i);
            if (neighbor_i && tet.barycentric_coord(neighbor_i.f(i,0), neighbor_i.f(i,1), neighbor_i.f(i,2), i) < 1e-7) { continue; }
            if (!boundary_tet_map.has(tet.uid)) {
                boundary_tet_map.set(tet.uid, tet);
                boundary_tet_points_map.set(tet.uid, []);
            }
            const i0 = (i+1)%4; const i1 = (i+2)%4; const i2 = (i+3)%4;
            const pts = [toFormantDataIndex(tet.x[i0][0], tet.x[i0][1], tet.x[i0][2]),
                         toFormantDataIndex(tet.x[i1][0], tet.x[i1][1], tet.x[i1][2]),
                         toFormantDataIndex(tet.x[i2][0], tet.x[i2][1], tet.x[i2][2])];
            boundary_tet_points_map.get(tet.uid).push([i, tet.normal(i)[2], pts]);
        }
    });
    let boundary_components = [[new Set(), new Set()],
                               [new Set(), new Set()],
                               [new Set(), new Set()]];
    const seen_points = new Set();
    for (const uid1 of boundary_tet_map.keys()) {
        for (const [i1, n1, pts1] of boundary_tet_points_map.get(uid1)) {
            // For every point on the relevant face...
            for (const pt1 of pts1) {
                if (seen_points.has(pt1)) { continue; }
                seen_points.add(pt1);
                const f1 = rawFormantData[pt1|0];
                const f2 = rawFormantData[pt1|1];
                const f3 = rawFormantData[pt1|2];
                const f4 = rawFormantData[pt1|2];
                // Evaluate max and min
                if (!minF1 || minF1 > f1) { minF1 = f1; }
                if (!minF2 || minF2 > f2) { minF2 = f2; }
                if (!minF3 || minF3 > f3) { minF3 = f3; }
                if (!minF4 || minF4 > f4) { minF4 = f4; }
                if (!maxF1 || maxF1 < f1) { maxF1 = f1; }
                if (!maxF2 || maxF2 < f2) { maxF2 = f2; }
                if (!maxF3 || maxF3 < f3) { maxF3 = f3; }
                if (!maxF4 || maxF4 < f4) { maxF4 = f4; }
                // Collect all intersections with other faces
                let intersections = [];
                for (const uid2 of boundary_tet_map.keys()) {
                    for (const [i2, n2, pts2] of boundary_tet_points_map.get(uid2)) {
                        if (pts2.includes(pt1)) { continue; }
                        const intersection = boundary_tet_map.get(uid2).intersect_face_f3(f1, f2, i2);
                        if (intersection) {
                            intersections.push([intersection, n2]);
                        }
                    }
                }
                intersections.sort((a,b) => a[0] - b[0]);
                // Count how many net crossings and complete volumes are above and below our point
                let count_below = 0;
                let complete_below = 0;
                for (let i = 0; i < intersections.length && intersections[i][0] < f3; i++) {
                    count_below -= Math.sign(intersections[i][1]);
                    if (count_below == 0) { complete_below = 1; }
                }
                let count_above = 0;
                let complete_above = 0;
                for (let i = intersections.length-1; i >= 0 && intersections[i][0] > f3; i--) {
                    count_above += Math.sign(intersections[i][1]);
                    if (count_above == 0) { complete_above = 1; }
                }
                // Fix a weird case
                if (count_below == 0 && count_above == 0 && complete_below == 1 && complete_above == 1) {
                    count_below = 1;
                    complete_below = 0;
                }
                // Add to boundary parts
                if (count_below == 0) {
                    boundary_components[0][complete_below].add(pt1);
                }
                if (count_below == 0 && count_above == 0) {
                    boundary_components[1][Math.max(complete_below, complete_above)].add(pt1);
                }
                if (count_above == 0) {
                    boundary_components[2][complete_above].add(pt1);
                }
            }
        }
    }
    boundary_components = boundary_components.map((a) => a.map((s) => Array.from(s.values())));
    boundary_components[1][0].sort(function (a,b) {
        const a_f1 = rawFormantData[a|0];
        const b_f1 = rawFormantData[b|0];
        const a_f2 = rawFormantData[a|1];
        const b_f2 = rawFormantData[b|1];
        const a_t = Math.atan2(a_f2 - (maxF2 + minF2) / 2, a_f1 - (maxF1 + minF1) / 2);
        const b_t = Math.atan2(b_f2 - (maxF2 + minF2) / 2, b_f1 - (maxF1 + minF1) / 2);
        return a_t - b_t;
    })
    this.saveAs("formantDataBoundary", "boundaryComponents", boundary_components,
                                       [["minF1", minF1], ["minF2", minF2], ["minF3", minF3], ["minF4", minF4],
                                        ["maxF1", maxF1], ["maxF2", maxF2], ["maxF3", maxF3], ["maxF4", maxF4]]);
}

var formantDataCollection = new FormantDataCollection();

const f3_region_width = 320;
const f3_region_start = 0;
const f3_region_end = f3_region_width - 170;

var Analysis = {
    M: 4096,

    init : function()
    {
        spaceCanvas.addEventListener('mousedown', Analysis.mouseMove.bind(this, true));
        spaceCanvas.addEventListener('mousemove', Analysis.mouseMove.bind(this, false));
        spaceCanvas.addEventListener('mouseup', Analysis.mouseLeave.bind(this));
        spaceCanvas.addEventListener('mouseleave',Analysis.mouseLeave.bind(this));

        allDiv.addEventListener("mousedown", function () {
            if (!AudioSystem.started)
            {
                AudioSystem.started = true;
                AudioSystem.startSound();
                // remove class: dimmed
                allDiv.className = "";
                // remove class: click text
                containersDiv.className = "containers";
            }
        });
    },


    mouseMove: function(mouseDown, event)
    {
        this.mouse_pageX = event.pageX;
        this.mouse_pageY = event.pageY;
        if (mouseDown) {
            this.mouseDown = true;
            this.orig_f1 = sincFFTPeaks[0];
            this.orig_f2 = sincFFTPeaks[1];
            this.orig_f3 = sincFFTPeaks[2];
            this.mouseStep();
        }
    },

    mouseStep: function() {
        if (this.mouseDown) {
            const rect = spaceCanvas.getBoundingClientRect();
            const x = (this.mouse_pageX-window.scrollX-rect.left)/rect.width*spaceCanvas.width;
            const y = (this.mouse_pageY-window.scrollY-rect.top)/rect.height*spaceCanvas.height;
            const f1 = x <= f3_region_width ? this.orig_f1 : fromNormalizedFormant(     y                   /spaceCanvas.height, 1);
            const f2 = x <= f3_region_width ? this.orig_f2 : fromNormalizedFormant(1 - (x - f3_region_width)/spaceCanvas.height, 2);
            const f3 = x >  f3_region_width ? this.orig_f3 : fromNormalizedFormant(1 -  y                   /spaceCanvas.height, 3);
            if (!this.mouseTet) {
                const index = fromLipTongueValues(TractUI.lipDiameter, TractUI.tongueDiameter, TractUI.tongueIndex);
                this.mouseTet = FormantTet.withVertexNearestTo([f1, f2, f3], index);
            }
            this.mouseTet = this.mouseTet.nearest(f1, f2, f3);
            if (!this.mouseTet.contains(f1, f2, f3)) {
                this.mouseTet = FormantTet.withNearestBoundaryVertexTo(f1, f2, f3);
            }
            // this.mouseTet = this.mouseTet.nearest_f1_f2(f1, f2, f3);
            const is = this.mouseTet.interp_coords(f1, f2, f3);
            const [lipDiameter, tongueDiameter, tongueIndex] = toLipTongueValues(is[0], is[2], is[1]);
            clickIndexDiameter(TractUI.lipIndex, lipDiameter);
            clickIndexDiameter(tongueIndex, tongueDiameter);
        }
    },

    mouseLeave: function(event)
    {
        this.mouseDown = false;
        this.mouse_pageX = undefined;
        this.mouse_pageY = undefined;
        this.mouseTet = undefined;
    },

    analyze : function(outArrayIn, middlePeakOffsetIn)
    {
        outArray = outArrayIn;
        middlePeakOffset = middlePeakOffsetIn;

        // console.log(TractUI.lipDiameter, TractUI.tongueDiameter, TractUI.tongueIndex)
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

        this.draw();
    },

    draw : function()
    {
        amplCtx.clearRect(0, 0, amplCanvas.width, amplCanvas.height);
        amplCtx.lineWidth = 3;
        amplCtx.strokeStyle = "orchid";
        amplCtx.beginPath()
        const iOff = this.M/4 - middlePeakOffset;
        for (let i = 0; i < this.M/2; i += 1) {
            const x = 2 * i / this.M * amplCanvas.width
            const y = (0.5 - 1.25 * outArray[iOff + i]) * amplCanvas.height;
            amplCtx.lineTo(x, y)
        }
        amplCtx.stroke()

        freqCtx.clearRect(0, 0, freqCanvas.width, freqCanvas.height);
        freqCtx.lineWidth = 3;
        freqCtx.strokeStyle = oppOrchid;
        freqCtx.beginPath()
        for (let i = 0; i < this.M/4; i += 1) {
            const x = Math.log2(1 + 4 * i / this.M) * freqCanvas.width
            const h = Math.log2(Math.max(1e-10, pulseFFTArray[i])) / 18;
            const y = (0.5 - Math.clamp(h, -0.5, 0.5)) * freqCanvas.height;
            freqCtx.lineTo(x, y)
        }
        freqCtx.stroke()
        freqCtx.lineWidth = 3;
        freqCtx.strokeStyle = "orchid";
        freqCtx.fillStyle = palePink;
        freqCtx.beginPath()
        for (let i = 0; i < this.M/4; i += 1) {
            const x = Math.log2(1 + 4 * i / this.M) * freqCanvas.width;
            const h = Math.log2(Math.max(1e-10, outFFTArray[i])) / 18;
            const y = (0.5 - Math.clamp(h, -0.5, 0.5)) * freqCanvas.height;
            freqCtx.lineTo(x, y)
        }
        freqCtx.lineTo(freqCanvas.width+1, freqCanvas.height+1);
        freqCtx.lineTo(-1, freqCanvas.height+1);
        freqCtx.fill()
        freqCtx.stroke()

        filtCtx.clearRect(0, 0, filtCanvas.width, filtCanvas.height);
        filtCtx.lineWidth = 3;
        filtCtx.strokeStyle = oppOrchid;
        filtCtx.fillStyle = oppOrchid;
        filtCtx.beginPath()
        filtCtx.lineTo(-1, filtCanvas.height+1);
        const truncSincFFTArray = sincFFTArray.slice(0, this.M/4)
        const truncSincFFTMax = Math.max(...truncSincFFTArray);
        normSincFFTArray = truncSincFFTArray.map((x) => x / truncSincFFTMax);
        for (let i = 0; i < this.M/4; i += 1) {
            const x = Math.log2(1 + 4 * i / this.M) * filtCanvas.width
            const y = filtCanvas.height * (1 - 0.8 * normSincFFTArray[i]);
            filtCtx.lineTo(x, y)
        }
        filtCtx.lineTo(filtCanvas.width+1, filtCanvas.height+1);
        filtCtx.stroke()
        filtCtx.globalAlpha = 0.15;
        filtCtx.fill()
        filtCtx.globalAlpha = 1;
        filtCtx.lineWidth = 2;
        filtCtx.strokeStyle = "black";
        sincFFTPeaks = find_peaks(normSincFFTArray).filter(function ([_, x_range, dy_range]) {
            // Note: This numbers were found experimentally with the goal of not having F3 split into two for rounded front vowels
            return x_range >= 6.6857944001708915 || dy_range >= 0.005642906942830599;
        });
        const r0 = 36.436593939856664;
        const r1 = 22.14870729719607;
        if (sincFFTPeaks.length >= 4 && sincFFTPeaks[4-1][0] > 400) {
            if (sincFFTPeaks[2-1][0] < 150) {
                sincFFTPeaks.splice(3-1, 0, sincFFTPeaks[3-1]);
                let [x, x_range, dy_range] = sincFFTPeaks[3-1];
                const r = Math.pow(r1, x_range / r0);
                sincFFTPeaks[3-1] = [x - r/2, r, dy_range];
                sincFFTPeaks[4-1] = [x + r/2, r, dy_range];
            }
            else {
                sincFFTPeaks.splice(2-1, 0, sincFFTPeaks[2-1]);
                let [x, x_range, dy_range] = sincFFTPeaks[2-1];
                const r = Math.pow(r1, x_range / r0);
                didIt = x_range;
                sincFFTPeaks[2-1] = [x - r/2, r, dy_range];
                sincFFTPeaks[3-1] = [x + r/2, r, dy_range];
            }
        } 
        sincFFTPeaks = sincFFTPeaks.map((x) => x[0]);
        // console.log(`${sincFFTPeaks[0]} ${sincFFTPeaks[1]} ${sincFFTPeaks[2]} ${sincFFTPeaks[3]}`)

        filtCtx.font = "24px Arial";
        filtCtx.fillStyle = "#AAA";
        filtCtx.textAlign = "center";
        filtCtx.textBaseline = 'top';
        for (let j = 0; j < sincFFTPeaks.length; j += 1) {
            const i = sincFFTPeaks[j];
            const x = Math.log2(1 + 4 * i / this.M) * filtCanvas.width
            filtCtx.beginPath()
            filtCtx.moveTo(x, filtCanvas.height * 0.2)
            filtCtx.lineTo(x, filtCanvas.height)
            filtCtx.stroke()
            filtCtx.fillText(`F${j+1}`, x, filtCanvas.height * 0.1);
        }

        for (let i = 0; i < 4; i += 1) {
            let f = sincFFTPeaks[i] * sampleRate / this.M;
            const str = `F${i+1} = ${f.toFixed(3)} Hz`;
            if (str != document.getElementById(`f${i+1}Div`).innerHTML) {
                document.getElementById(`f${i+1}Div`).innerHTML = str;
            }
        }

        spaceCtx.clearRect(0, 0, spaceCanvas.width, spaceCanvas.height);

        spaceCtx.lineWidth = 45*2.5;
        spaceCtx.fillStyle = oppPalePink;
        spaceCtx.strokeStyle = oppPalePink;
        spaceCtx.lineCap = "round";
        spaceCtx.lineJoin = "round";
        spaceCtx.beginPath()
        for (let i = 0; i <= boundaryComponents[1][0].length; i++) {
            const j = boundaryComponents[1][0][i % boundaryComponents[1][0].length];
            const x = (1 - toNormalizedFormant(rawFormantData[j + 1], 2)) * spaceCanvas.height;
            const y = (    toNormalizedFormant(rawFormantData[j + 0], 1)) * spaceCanvas.height;
            spaceCtx.lineTo(f3_region_width + x, y)
        }
        spaceCtx.stroke()
        spaceCtx.fill()

        if (sincFFTPeaks.length > 2) {
            const p = FormantTet.treePoint(sincFFTPeaks[0], sincFFTPeaks[1], sincFFTPeaks[2]);
            const f3_bounds = FormantTet.boundaryIntersections(p);
            let pts = [[f3_bounds[0][0], f3_bounds[f3_bounds.length-1][0]]];
            if (f3_bounds.length > 2 &&
                f3_bounds[1][0][2] > f3_bounds[0][0][2] &&
                f3_bounds[f3_bounds.length-2][0][2] < f3_bounds[f3_bounds.length-1][0][2]) {
                pts = [[f3_bounds[0][0], f3_bounds[1][0]],
                       [f3_bounds[f3_bounds.length-2][0], f3_bounds[f3_bounds.length-1][0]]];
            }
            for (const [minPoint, maxPoint] of pts) {
                const minTet = FormantTet.withVertexNearestTo(p, minPoint.index, minPoint);
                const maxTet = FormantTet.withVertexNearestTo(p, maxPoint.index, maxPoint);
                const minIntersections = minTet.intersect_f3(sincFFTPeaks[0], sincFFTPeaks[1]);
                const maxIntersections = maxTet.intersect_f3(sincFFTPeaks[0], sincFFTPeaks[1]);
                minIntersections.sort((a,b) => a[0] - b[0]);
                maxIntersections.sort((a,b) => b[0] - a[0]);
                const minBound = minIntersections.length > 0 ? toNormalizedFormant(minIntersections[0][0],3) : minPoint[2];
                const maxBound = maxIntersections.length > 0 ? toNormalizedFormant(maxIntersections[0][0],3) : maxPoint[2];
                spaceCtx.beginPath()
                let x = (f3_region_start + f3_region_end) / 2;
                let y = (1 - minBound) * spaceCanvas.height;
                spaceCtx.lineTo(x, y)
                y = (1 - maxBound) * spaceCanvas.height;
                spaceCtx.lineTo(x, y)
                spaceCtx.stroke()
            }
        }

        if (TractUI.tongueTouch != 0 || TractUI.lipTouch != 0) {
            // palePink = (255, 238, 245)
            // a*bi +  = ai - 255*(1-a)
            const alpha = 1;
            const r = (255 - 255*(1-alpha)) / alpha;
            const g = (238 - 255*(1-alpha)) / alpha;
            const b = (245 - 255*(1-alpha)) / alpha;
            spaceCtx.fillStyle = `rgb(${r}, ${g}, ${b})`;
            spaceCtx.strokeStyle = `rgb(${r}, ${g}, ${b})`;
            spaceCtx.globalAlpha = alpha;
            const [iLipDiam, iTongueY, iTongueX] = fromLipTongueValues(TractUI.lipDiameter, TractUI.tongueDiameter, TractUI.tongueIndex);
            const iLipDiam0 = Math.floor(iLipDiam);
            const iLipDiam1 = Math.ceil(iLipDiam);
            const t = iLipDiam % 1;
            for (let iTongueX = 0; iTongueX < 0xF; iTongueX++) {
                for (let iTongueY = 0; iTongueY < 0xF; iTongueY++) {
                    spaceCtx.beginPath();
                    for (let k = 0; k < 5; k++) {
                        const iTongueX_k = iTongueX + (((k % 4) >> 1) ^ ((k % 4) &  1));
                        const iTongueY_k = iTongueY + ((k % 4) >> 1);
                        const f1 = formantData(iLipDiam0, iTongueX_k, iTongueY_k, 0) * (1 - t) +
                                   formantData(iLipDiam1, iTongueX_k, iTongueY_k, 0) * t;
                        const f2 = formantData(iLipDiam0, iTongueX_k, iTongueY_k, 1) * (1 - t) +
                                   formantData(iLipDiam1, iTongueX_k, iTongueY_k, 1) * t;
                        const x = (1 - toNormalizedFormant(f2, 2)) * spaceCanvas.height;
                        const y = (    toNormalizedFormant(f1, 1)) * spaceCanvas.height;
                        spaceCtx.lineTo(f3_region_width + x, y, 10, 0, 2*Math.PI);
                    }
                    spaceCtx.fill();
                    spaceCtx.stroke();
                }
            }
            spaceCtx.beginPath()
            let x = (f3_region_start + f3_region_end) / 2;
            let y = (1 - toNormalizedFormant(sincFFTPeaks[2], 3)) * spaceCanvas.height;
            spaceCtx.arc(x, y, 45/2*2.5, 0, 2*Math.PI);
            spaceCtx.fill();
            spaceCtx.globalAlpha = 1.0;
        }

        // draw dashed lines
        spaceCtx.lineWidth = 2;
        spaceCtx.setLineDash([10,10]);
        spaceCtx.fillStyle = "#AAA";
        spaceCtx.strokeStyle = "#AAA";
        spaceCtx.font = "32px Arial";
        spaceCtx.globalCompositeOperation = "color-burn";
        const f1_tick = 100;
        const f2_tick = 500;
        const f3_tick = 100;
        const f1_factor = f1_tick * this.M / sampleRate;
        const f2_factor = f2_tick * this.M / sampleRate;
        const f3_factor = f3_tick * this.M / sampleRate;
        const f1_axis_start = 80;
        const f2_axis_start = (1 - toNormalizedFormant(f2_factor * (Math.ceil(minF2 / f2_factor) - 1), 2)) * spaceCanvas.height;
        spaceCtx.textAlign = "center";
        spaceCtx.textBaseline = 'bottom';
        spaceCtx.save();
        spaceCtx.translate(f3_region_width + f2_axis_start + 80, (f1_axis_start + spaceCanvas.height) / 2);
        spaceCtx.rotate(Math.PI/2);
        spaceCtx.fillText("F1 (Hz)", 0, 0);
        spaceCtx.restore();
        spaceCtx.save();
        spaceCtx.translate(f3_region_end + 100, (f1_axis_start + spaceCanvas.height) / 2);
        spaceCtx.rotate(Math.PI/2);
        spaceCtx.fillText("F3 (Hz)", 0, 0);
        spaceCtx.restore();
        spaceCtx.textAlign = "left";
        spaceCtx.textBaseline = 'middle';
        for (let f1 = Math.floor(minF1 / f1_factor); f1 < Math.ceil(maxF1 / f1_factor) + 1; f1++) {
            const y = toNormalizedFormant(f1_factor * f1, 1) * spaceCanvas.height;
            if (y > f1_axis_start) {
                spaceCtx.beginPath();
                spaceCtx.lineTo(f3_region_width, y);
                spaceCtx.lineTo(f3_region_width + f2_axis_start, y);
                spaceCtx.stroke();
            }
            if (y >= f1_axis_start) {
                spaceCtx.fillText(`${f1_tick * f1}`, f3_region_width + f2_axis_start + 8, y + 2);
            }
        }
        for (let f3 = Math.ceil(minF3 / f3_factor) - 1; f3 < Math.ceil(maxF3 / f3_factor) + 1; f3++) {
            const y = (1 - toNormalizedFormant(f3_factor * f3, 3)) * spaceCanvas.height;
            if (y > f1_axis_start) {
                spaceCtx.beginPath();
                spaceCtx.lineTo(f3_region_start, y);
                spaceCtx.lineTo(f3_region_end, y);
                spaceCtx.stroke();
            }
            if (y >= f1_axis_start) {
                spaceCtx.fillText(`${f3_tick * f3}`, f3_region_end + 8, y + 2);
            }
        }
        spaceCtx.textAlign = "center";
        spaceCtx.textBaseline = 'top';
        spaceCtx.fillText("F2 (Hz)", f3_region_width + f2_axis_start / 2, 0);
        spaceCtx.textBaseline = 'bottom';
        for (let f2 = Math.ceil(minF2 / f2_factor) - 1; f2 < Math.ceil(maxF2 / f2_factor) + 1; f2++) {
            const x = (1 - toNormalizedFormant(f2_factor * f2, 2)) * spaceCanvas.height;
            if (x < f2_axis_start) {
                spaceCtx.beginPath();
                spaceCtx.lineTo(f3_region_width + x, f1_axis_start);
                spaceCtx.lineTo(f3_region_width + x, spaceCanvas.height);
                spaceCtx.stroke();
            }
            if (x <= f2_axis_start) {
                spaceCtx.fillText(`${f2_tick * f2}`, f3_region_width + x, f1_axis_start - 8);
            }
        }
        // draw boundary arrows
        spaceCtx.setLineDash([]);
        spaceCtx.lineWidth = 4;
        spaceCtx.beginPath();
        spaceCtx.lineTo(f3_region_width + 10, f1_axis_start-10);
        spaceCtx.lineTo(f3_region_width + 0, f1_axis_start);
        spaceCtx.lineTo(f3_region_width + 10, f1_axis_start+10);
        spaceCtx.lineTo(f3_region_width + 0, f1_axis_start);
        spaceCtx.lineTo(f3_region_width + f2_axis_start, f1_axis_start);
        spaceCtx.lineTo(f3_region_width + f2_axis_start, spaceCanvas.height - 2);
        spaceCtx.lineTo(f3_region_width + f2_axis_start - 10, spaceCanvas.height - 2 - 10);
        spaceCtx.lineTo(f3_region_width + f2_axis_start, spaceCanvas.height - 2);
        spaceCtx.lineTo(f3_region_width + f2_axis_start + 10, spaceCanvas.height - 2 - 10);
        spaceCtx.stroke();
        spaceCtx.beginPath();
        spaceCtx.lineTo(f3_region_end, spaceCanvas.height);
        spaceCtx.lineTo(f3_region_end, f1_axis_start - 2);
        spaceCtx.lineTo(f3_region_end - 10, f1_axis_start - 2 + 10);
        spaceCtx.lineTo(f3_region_end, f1_axis_start - 2);
        spaceCtx.lineTo(f3_region_end + 10, f1_axis_start - 2 + 10);
        spaceCtx.stroke();

        spaceCtx.globalCompositeOperation = "source-over";

        // draw IPA symbols
        spaceCtx.fillStyle = "#474747";
        spaceCtx.strokeStyle = "#888";
        spaceCtx.font = "64px Arial";
        spaceCtx.textAlign = "center";
        spaceCtx.textBaseline = 'middle';
        spaceCtx.lineWidth = 3;
        const ipa = [["i", 293.189, 2752.617, [1,2]], ["y", 290.559, 1939.181, [3]],
                     ["e", 471.685, 2217.909, [3,4]], ["ø", 458.650, 1684.690 , [5]],
                     ["ɛ", 627.898, 1899.734, [5,6]], ["œ", 617.047, 1474.398, []],
                     ["æ", 755.531, 1574.874, [14]], ["ɐ", 55.63831, 111.10644, []],

                     ["u", 280.814, 587.882, [9,10]], ["ɯ", 289.417, 836.786, [11]],
                     ["o", 427.646, 674.643, [11,12]], ["ɤ", 439.720, 942.954, [13]],
                     ["ɔ", 572.700, 770.287, [13,15]], ["ʌ", 600.275, 1029.261, []],
                     ["a", 811.613, 1228.786, []], ["ɒ", 722.026, 949.818, []],

                     ["ɪ", 383.169, 2116.273, []],
                     ["ɵ", 361.728, 1286.303 , []],
                     ["ʊ", 377.588, 751.881, []],
                     ["ə", 552.198, 1242.380, []],
                     ["ɐ", 705.832, 1232.888, []]];
        for (const [str, f1, f2, lines] of ipa) {
            const x = (1 - toNormalizedFormant(f2 * this.M / sampleRate, 2)) * spaceCanvas.height;
            const y = (    toNormalizedFormant(f1 * this.M / sampleRate, 1)) * spaceCanvas.height;
            spaceCtx.fillText(str, f3_region_width + x, y);
            for (const i of lines) {
                spaceCtx.beginPath();
                const xi = (1 - toNormalizedFormant(ipa[i][2] * this.M / sampleRate, 2)) * spaceCanvas.height;
                const yi = (    toNormalizedFormant(ipa[i][1] * this.M / sampleRate, 1)) * spaceCanvas.height;
                const r = 0.2;
                spaceCtx.lineTo(f3_region_width +      r  * x + (1 - r) * xi,      r  * y + (1 - r) * yi);
                spaceCtx.lineTo(f3_region_width + (1 - r) * x +      r  * xi, (1 - r) * y +      r  * yi);
                spaceCtx.stroke();
            }
        }
        
        if (sincFFTPeaks.length >= 3) {
            spaceCtx.fillStyle = oppOrchid;
            x = (1 - toNormalizedFormant(sincFFTPeaks[1], 2)) * spaceCanvas.height;
            y = (    toNormalizedFormant(sincFFTPeaks[0], 1)) * spaceCanvas.height;

            spaceCtx.lineWidth = 4*2.5;
            spaceCtx.strokeStyle = oppOrchid;

            spaceCtx.globalAlpha = 0.7;
            spaceCtx.beginPath();
            spaceCtx.arc(f3_region_width + x, y, 18*2.5, 0, 2*Math.PI);
            spaceCtx.stroke();        
            spaceCtx.globalAlpha = 0.15;
            spaceCtx.fill();
            spaceCtx.globalAlpha = 1.0;

            x = (f3_region_start + f3_region_end) / 2;
            y = (1 - toNormalizedFormant(sincFFTPeaks[2], 3)) * spaceCanvas.height;

            spaceCtx.globalAlpha = 0.7;
            spaceCtx.beginPath();
            spaceCtx.arc(x, y, 18*2.5, 0, 2*Math.PI);
            spaceCtx.stroke();        
            spaceCtx.globalAlpha = 0.15;
            spaceCtx.fill();
            spaceCtx.globalAlpha = 1.0;
        }

        if (this.mouseDown) {
            this.mouseStep();
        }

        formantDataCollection.step(sincFFTPeaks);
    }
}

Analysis.init();
requestAnimationFrame(Analysis.draw.bind(Analysis));
