
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

// Source: https://github.com/stonkpunk/my-npm-modules/blob/main/barycentric-coordinates/index.js
function determinant(m00,m01,m02,m03,
                     m10,m11,m12,m13,
                     m20,m21,m22,m23,
                     m30,m31,m32,m33){ //http://www.euclideanspace.com/maths/algebra/matrix/functions/determinant/fourD/index.htm
    return  m03 * m12 * m21 * m30 - m02 * m13 * m21 * m30-
            m03 * m11 * m22 * m30 +m01 * m13 * m22 * m30+
            m02 * m11 * m23 * m30 -m01 * m12 * m23 * m30-
            m03 * m12 * m20 * m31 +m02 * m13 * m20 * m31+
            m03 * m10 * m22 * m31 -m00 * m13 * m22 * m31-
            m02 * m10 * m23 * m31 +m00 * m12 * m23 * m31+
            m03 * m11 * m20 * m32 -m01 * m13 * m20 * m32-
            m03 * m10 * m21 * m32 +m00 * m13 * m21 * m32+
            m01 * m10 * m23 * m32 -m00 * m11 * m23 * m32-
            m02 * m11 * m20 * m33 +m01 * m12 * m20 * m33+
            m02 * m10 * m21 * m33 -m00 * m12 * m21 * m33-
            m01 * m10 * m22 * m33 +m00 * m11 * m22 * m33;
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

var allCanvasContainers = document.querySelectorAll(".canvas-container");
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

var pulseFFTArray = [];
var sincFFTArray = [];
var normSincFFTArray = [];
var sincFFTPeaks = [];
var outFFTArray = [];
function copyData(toCopy = sincFFTArray) {
    navigator.clipboard.writeText(toCopy.slice(0, Analysis.M/4 - 4));
}

function copyPeakData() {
    let obj = { peaks: sincFFTPeaks, spectrum: normSincFFTArray };
    navigator.clipboard.writeText(JSON.stringify(obj));
}

function clickIndexDiameter(index, diameter) {
    // getIndex : function(x,y)
    // {
    //     var xx = this.originX-x; var yy = y-this.originY;
    //     var angle = Math.atan2(yy, xx);
    //     while (angle> 0) angle -= 2*Math.PI;
    //     return (Math.PI + angle - this.angleOffset)*(Tract.lipStart-1) / (this.angleScale*Math.PI);
    // },
    // getDiameter : function(x,y)
    // {
    //     var xx = this.originX-x; var yy = y-this.originY;
    //     return (this.radius-Math.sqrt(xx*xx + yy*yy))/this.scale;
    // },
    const r = TractUI.radius - diameter * TractUI.scale;
    const t = - TractUI.angleScale * Math.PI * index / (Tract.lipStart - 1) - TractUI.angleOffset - Math.PI;
    const xx = r * Math.cos(t);
    const yy = r * Math.sin(t);
    const x = TractUI.originX-xx;
    const y = TractUI.originY-yy;
    // touch.x = (event.pageX-window.scrollX-UI.left_margin)/UI.width*600;
    // touch.y = (event.pageY-window.scrollY-UI.top_margin)/UI.width*600;
    const pageX = (x + window.scrollX + UI.left_margin) * UI.width / 600;
    const pageY = (y + window.scrollY + UI.top_margin) * UI.width / 600;
    const event = { pageX: pageX, pageY: pageY };
    UI.mouseDown = true;
    UI.startMouse(event);
    UI.mouseDown = false;
    UI.endMouse(event);
}

// Formant data

function formantData(iLipDiam, iTongueY, iTongueX, iFormant) {
    const i = (iLipDiam << 10) | (iTongueY << 6) | (iTongueX << 2);
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
    const detC =  this.f(i0,j0) * (this.f(i1,j1) * this.f(i2,j2) - this.f(i2,j1) * this.f(i1,j2))
                - this.f(i0,j1) * (this.f(i1,j0) * this.f(i2,j2) - this.f(i1,j2) * this.f(i2,j0))
                + this.f(i0,j2) * (this.f(i1,j0) * this.f(i2,j1) - this.f(i1,j1) * this.f(i2,j0));
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

FormantTet.prototype.barycentric_coords = function barycentric_coords(f1, f2, f3) {
    return [this.barycentric_coord(f1, f2, f3, 0),
            this.barycentric_coord(f1, f2, f3, 1),
            this.barycentric_coord(f1, f2, f3, 2),
            this.barycentric_coord(f1, f2, f3, 3)];
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

FormantTet.prototype.nearest = function nearest(f1, f2, f3, seen) {
    if (!seen) { seen = new Set(); }
    if (seen.has(this.uid)) { return null; }
    seen.add(this.uid);

    let coords = this.barycentric_coords(f1, f2, f3);
    coords = coords.map((x,i) => [x,i]);
    coords.sort((a,b) => a[0] - b[0]);
    if (coords[0][0] < 0) {
        const neighbor_i = this.neighbor(coords[0][1]);
        if (neighbor_i) {
            const nearest = neighbor_i.nearest(f1, f2, f3, seen);
            if (nearest) {
                return nearest;                
            }
        }
    }
    return this;
}

FormantTet.treePoint = function treePoint(f1, f2, f3, iLipDiam, iTongueX, iTongueY) {
    const p = [(f1 - minF1) / (maxF1 - minF1),
               (f2 - minF2) / (maxF2 - minF2),
               (f3 - minF3) / (maxF3 - minF3)];
    if (iLipDiam !== undefined) {
        p.index = [iLipDiam, iTongueX, iTongueY];
    }
    return p;
}

FormantTet.treePointFromIndex = function treePointFromIndex(_, i) {
    return FormantTet.treePoint(rawFormantData[(i << 2) | 0],
                                rawFormantData[(i << 2) | 1],
                                rawFormantData[(i << 2) | 2],
                                (i >> 8), (i >> 4) & 0xF, i & 0xF);
}

FormantTet.treeDistance = function treeDistance(a, b) {
    return Math.pow(a[0] - b[0], 2) + Math.pow(a[1] - b[1], 2) + Math.pow(a[2] - b[2], 2);
}

FormantTet.tree = new kdTree(Array.from({ length: 0x1000 }, FormantTet.treePointFromIndex),
                             FormantTet.treeDistance, [0,1,2]);

FormantTet.withNearestVertexTo = function withNearestVertexTo(f1, f2, f3) {
    const p = FormantTet.treePoint(f1, f2, f3);
    const nearest = FormantTet.tree.nearest(p, 1)[0][0];
    const x = [[nearest.index[0], nearest.index[1], nearest.index[2], 1.0],
               [nearest.index[0], nearest.index[1], nearest.index[2], 1.0],
               [nearest.index[0], nearest.index[1], nearest.index[2], 1.0],
               [nearest.index[0], nearest.index[1], nearest.index[2], 1.0]];
    x[0][0] += p[0] > nearest[0] && x[0][0] < 0xF || x[0][0] == 0 ? 1 : -1;
    x[1][1] += p[1] > nearest[1] && x[1][1] < 0xF || x[1][1] == 0 ? 1 : -1;
    x[2][2] += p[2] > nearest[2] && x[2][2] < 0xF || x[2][2] == 0 ? 1 : -1;
    return new FormantTet(3, x, [null, null, null, null]);
}

FormantTet.nearest = function nearest(f1, f2, f3) {
    return FormantTet.withNearestVertexTo(f1, f2, f3).nearest(f1, f2, f3);
}

function FormantDataCollection() {
    this.wait = -1;
}

FormantDataCollection.prototype.start = function start() {
    if (this.wait < 0) {
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
    this.i = toFormantDataIndex(iLipDiam, iTongueX, iTongueY, 0);
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
    formantData[this.i + 0] = peaks[0];
    formantData[this.i + 1] = peaks[1];
    formantData[this.i + 2] = peaks[2];
    formantData[this.i + 3] = peaks[3];
    if (this.iHilbert < 0x1000) {
        this.iHilbert += 1;
        this.update();
        this.wait = 1;
    }
    else {
        const obj = { formantData: Array.from(formantData), 
                      innerLipControlRadius: TractUI.innerLipControlRadius,
                      outerLipControlRadius: TractUI.outerLipControlRadius,
                      outerTongueControlRadius: TractUI.outerTongueControlRadius,
                      innerTongueControlRadius: TractUI.innerTongueControlRadius,
                      tongueUpperIndexBound: TractUI.tongueUpperIndexBound,
                      tongueLowerIndexBound: TractUI.tongueLowerIndexBound };
        const jsonBlob = new Blob([JSON.stringify(obj)], {type: "application/json"});
        const jsBlob = new Blob([`var rawFormantData = [${obj.formantData}];`], {type: "text/javascript"});
        saveAs(jsonBlob, `formantData.json`);
        saveAs(jsBlob, `formantData.js`);
        this.wait = -1;
    }
}

var formantDataCollection = new FormantDataCollection();

var Analysis = {
    M: 4096,

    init : function()
    {
        const log2MinF1 = Math.log2(minF1) - 0.125;
        const log2MinF2 = Math.log2(minF2) - 0.125;
        const log2MinF3 = Math.log2(minF3) - 0.125;
        const log2MaxF1 = Math.log2(maxF1) + 0.125;
        const log2MaxF2 = Math.log2(maxF2) + 0.125;
        const log2MaxF3 = Math.log2(maxF3) + 0.125;
        this.bF1 = 1 / (log2MaxF1 / log2MinF1 - 1);
        this.bF2 = 1 / (log2MaxF2 / log2MinF2 - 1);
        this.bF3 = 1 / (log2MaxF3 / log2MinF3 - 1);
        this.mF1 = this.bF1 / log2MinF1;
        this.mF2 = this.bF2 / log2MinF2;
        this.mF3 = this.bF3 / log2MinF3;

        tractCanvas.addEventListener('mousedown', Analysis.handleMouse);
        tractCanvas.addEventListener('mousemove', Analysis.handleMouse);

        allCanvasContainers.forEach((canvasContainer) => canvasContainer.addEventListener("mousedown", function () {
            if (!AudioSystem.started)
            {
                AudioSystem.started = true;
                AudioSystem.startSound();
                allCanvasContainers.forEach((iCanvasContainer) => iCanvasContainer.className = "canvas-container");
            }
        }));
    },

    handleMouse: function(event)
    {

    },

    draw : function(outArray, middlePeakOffset)
    {
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

        amplCtx.clearRect(0, 0, amplCanvas.width, amplCanvas.height);
        amplCtx.strokeStyle = "black";
        amplCtx.beginPath()
        const iOff = this.M/4 - middlePeakOffset;
        for (let i = 0; i < this.M/2; i += 1) {
            const x = 2 * i / this.M * amplCanvas.width
            const y = (0.5 - 1.25 * outArray[iOff + i]) * amplCanvas.height;
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
        const truncSincFFTArray = sincFFTArray.slice(0, this.M/4)
        const truncSincFFTMax = Math.max(...truncSincFFTArray);
        normSincFFTArray = truncSincFFTArray.map((x) => x / truncSincFFTMax);
        for (let i = 0; i < this.M/4; i += 1) {
            const x = Math.log2(1 + 4 * i / this.M) * filtCanvas.width
            const y = filtCanvas.height * (1 - normSincFFTArray[i]);
            filtCtx.lineTo(x, y)
        }
        filtCtx.stroke()
        filtCtx.strokeStyle = "black";
        sincFFTPeaks = find_peaks(normSincFFTArray).filter(function ([_, x_range, dy_range]) {
            // Note: This numbers were found experimentally with the goal of not having F3 split into two for rounded front vowels
            return x_range >= 6.6857944001708915 || dy_range >= 0.005642906942830599;
        });
        const r0 = 36.436593939856664;
        const r1 = 22.14870729719607;
        if (sincFFTPeaks[4-1][0] > 400) {
            if (sincFFTPeaks[2-1][0] < 150) {
                sincFFTPeaks.splice(3-1, 0, sincFFTPeaks[3-1]);
                let [x, x_range, dy_range] = sincFFTPeaks[3-1];
                const r = Math.pow(r1, x_range / r0);
                sincFFTPeaks[3-1] = [x - r/2, r, dy_range];
                sincFFTPeaks[4-1] = [x + r/2, r, dy_range];
            }
            else {
                sincFFTPeaks.splice(2-1, 0, sincFFTPeaks[2-1]);
                [x, x_range, dy_range] = sincFFTPeaks[2-1];
                const r = Math.pow(r1, x_range / r0);
                didIt = x_range;
                sincFFTPeaks[2-1] = [x - r/2, r, dy_range];
                sincFFTPeaks[3-1] = [x + r/2, r, dy_range];
            }
        } 
        sincFFTPeaks = sincFFTPeaks.map((x) => x[0]);
        // console.log(`${sincFFTPeaks[0]} ${sincFFTPeaks[1]} ${sincFFTPeaks[2]} ${sincFFTPeaks[3]}`)

        for (let j = 0; j < sincFFTPeaks.length; j += 1) {
            const i = sincFFTPeaks[j];
            const x = Math.log2(1 + 4 * i / this.M) * filtCanvas.width
            filtCtx.beginPath()
            filtCtx.moveTo(x, 0)
            filtCtx.lineTo(x, filtCanvas.height)
            filtCtx.stroke()
        }

        let formant_strings = [];
        for (let i = 0; i < 4; i += 1) {
            let f = sincFFTPeaks[i] * sampleRate / this.M;
            formant_strings.push(`F${i+1} = ${f.toFixed(2)}`);
        }
        document.getElementById("formantDiv").innerHTML = formant_strings.join(", ");

        spaceCtx.clearRect(0, 0, spaceCanvas.width, spaceCanvas.height);
        
        spaceCtx.beginPath()
        spaceCtx.moveTo(spaceCanvas.height, spaceCanvas.height)
        const xboundary = (1 - this.mF2 * Math.log2(this.maxF1) + this.bF2) * spaceCanvas.height;
        const yboundary = (    this.mF1 * Math.log2(this.minF2) - this.bF1) * spaceCanvas.height;
        spaceCtx.lineTo(spaceCanvas.height, yboundary);
        spaceCtx.lineTo(xboundary, spaceCanvas.height);
        spaceCtx.fill();

        spaceCtx.fillStyle = palePink;
        spaceCtx.beginPath()
        for (let i = 0; i < formantDataBoundary.length; i++) {
            const j = formantDataBoundary[i] << 2;
            const x = (1 - this.mF2 * Math.log2(rawFormantData[j + 1]) + this.bF2) * spaceCanvas.height;
            const y = (    this.mF1 * Math.log2(rawFormantData[j + 0]) - this.bF1) * spaceCanvas.height;
            spaceCtx.lineTo(x, y)
        }
        spaceCtx.fill()

        if (TractUI.tongueTouch != 0 || TractUI.lipTouch != 0) {
            spaceCtx.fillStyle = "orchid";
            spaceCtx.strokeStyle = "orchid";
            spaceCtx.globalAlpha = 0.15;
            const [iLipDiam, iTongueY, iTongueX] = fromLipTongueValues(TractUI.lipDiameter, TractUI.  tongueDiameter, TractUI.tongueIndex);
            const iLipDiam0 = Math.floor(iLipDiam);
            const iLipDiam1 = Math.ceil(iLipDiam);
            const t = iLipDiam % 1;
            for (let iTongueX = 0; iTongueX < 0xF; iTongueX++) {
                for (let iTongueY = 0; iTongueY < 0xF; iTongueY++) {
                    spaceCtx.beginPath();
                    for (let k = 0; k < 4; k++) {
                        const iTongueX_k = iTongueX + ((k >> 1) ^ (k &  1));
                        const iTongueY_k = iTongueY + (k >> 1);
                        const x = (1 - this.mF2 * Math.log2(formantData(iLipDiam0, iTongueX_k, iTongueY_k, 1) * (1 - t) +
                                                            formantData(iLipDiam1, iTongueX_k, iTongueY_k, 1) * t       )    + this.bF2) * spaceCanvas.height;
                        const y = (    this.mF1 * Math.log2(formantData(iLipDiam0, iTongueX_k, iTongueY_k, 0) * (1 - t) +
                                                            formantData(iLipDiam1, iTongueX_k, iTongueY_k, 0) * t       )    - this.bF1) * spaceCanvas.height;
                        spaceCtx.lineTo(x, y, 10, 0, 2*Math.PI);
                    }
                    spaceCtx.fill();
                }
            }
            /*
            spaceCtx.lineCap = "round";
            spaceCtx.lineWidth = 10;
            spaceCtx.beginPath();
            const iTongueX0 = Math.floor(iTongueX);
            const iTongueY0 = Math.floor(iTongueY);
            const iTongueX1 = Math.ceil(iTongueX);
            const iTongueY1 = Math.ceil(iTongueY);
            const tX = iTongueX % 1;
            const tY = iTongueY % 1;
            for (let iLipDiam = 0; iLipDiam <= 0xF; iLipDiam++) {
                const j00 = toFormantDataIndex(iLipDiam, iTongueX0, iTongueY0, 0);
                const j01 = toFormantDataIndex(iLipDiam, iTongueX0, iTongueY1, 0);
                const j10 = toFormantDataIndex(iLipDiam, iTongueX1, iTongueY0, 0);
                const j11 = toFormantDataIndex(iLipDiam, iTongueX1, iTongueY1, 0);
                const x = (1 - this.mF2 * Math.log2(formantData[j00 + 1] * (1 - tX) * (1 - tY) + formantData[j01 + 1] * (1 - tX) * tY + formantData[j10 + 1] * tX * (1 - tY) + formantData[j11 + 1] * tX * tY)    + this.bF2) * spaceCanvas.height;
                const y = (    this.mF1 * Math.log2(formantData[j00 + 0] * (1 - tX) * (1 - tY) + formantData[j01 + 0] * (1 - tX) * tY + formantData[j10 + 0] * tX * (1 - tY) + formantData[j11 + 0] * tX * tY)    - this.bF1) * spaceCanvas.height;
                if (iLipDiam == 0) {
                    spaceCtx.moveTo(x, y);
                }
                else {
                    spaceCtx.lineTo(x, y);
                }
            }
            spaceCtx.stroke();
            */
            spaceCtx.globalAlpha = 1.0;
        }
        
        if (sincFFTPeaks.length >= 3) {
            spaceCtx.fillStyle = 'black';
            spaceCtx.beginPath();
            let x = (1 - this.mF2 * Math.log2(sincFFTPeaks[1]) + this.bF2) * spaceCanvas.height;
            let y = (    this.mF1 * Math.log2(sincFFTPeaks[0]) - this.bF1) * spaceCanvas.height;
            spaceCtx.arc(x, y, 10, 0, 2*Math.PI);
            spaceCtx.fill();
            spaceCtx.beginPath();
            x = spaceCanvas.height + (spaceCanvas.width - spaceCanvas.height) / 2;
            y = (this.mF3 * Math.log2(sincFFTPeaks[2]) - this.bF3) * spaceCanvas.height;
            const x_spread = (spaceCanvas.width - spaceCanvas.height) / 3;
            spaceCtx.arc(x, y, 10, 0, 2*Math.PI);
            spaceCtx.fill();
            /*
            const [iLipDiam, iTongueY, iTongueX] = fromFormants(sincFFTPeaks[0], sincFFTPeaks[1], sincFFTPeaks[2]);
            const iLipDiam0 = Math.floor(iLipDiam);
            const iTongueY0 = Math.floor(iTongueY);
            const iTongueX0 = Math.floor(iTongueX);
            const iLipDiam1 = Math.ceil(iLipDiam);
            const iTongueY1 = Math.ceil(iTongueY);
            const iTongueX1 = Math.ceil(iTongueX);
            const tL = iLipDiam % 1;
            const tY = iTongueY % 1;
            const tX = iTongueX % 1;
            const j000 = toFormantDataIndex(iLipDiam0, iTongueY0, iTongueX0, 0);
            const j001 = toFormantDataIndex(iLipDiam0, iTongueY0, iTongueX1, 0);
            const j010 = toFormantDataIndex(iLipDiam0, iTongueY1, iTongueX0, 0);
            const j011 = toFormantDataIndex(iLipDiam0, iTongueY1, iTongueX1, 0);
            const j100 = toFormantDataIndex(iLipDiam1, iTongueY0, iTongueX0, 0);
            const j101 = toFormantDataIndex(iLipDiam1, iTongueY0, iTongueX1, 0);
            const j110 = toFormantDataIndex(iLipDiam1, iTongueY1, iTongueX0, 0);
            const j111 = toFormantDataIndex(iLipDiam1, iTongueY1, iTongueX1, 0);
            const f0 = formantData[j000|0] * (1 - tL) * (1 - tY) * (1 - tX) +
                       formantData[j001|0] * (1 - tL) * (1 - tY) *      tX  +
                       formantData[j010|0] * (1 - tL) *      tY  * (1 - tX) +
                       formantData[j011|0] * (1 - tL) *      tY  *      tX  +
                       formantData[j100|0] *      tL  * (1 - tY) * (1 - tX) +
                       formantData[j101|0] *      tL  * (1 - tY) *      tX  +
                       formantData[j110|0] *      tL  *      tY  * (1 - tX) +
                       formantData[j111|0] *      tL  *      tY  *      tX;
            const f1 = formantData[j000|1] * (1 - tL) * (1 - tY) * (1 - tX) +
                       formantData[j001|1] * (1 - tL) * (1 - tY) *      tX  +
                       formantData[j010|1] * (1 - tL) *      tY  * (1 - tX) +
                       formantData[j011|1] * (1 - tL) *      tY  *      tX  +
                       formantData[j100|1] *      tL  * (1 - tY) * (1 - tX) +
                       formantData[j101|1] *      tL  * (1 - tY) *      tX  +
                       formantData[j110|1] *      tL  *      tY  * (1 - tX) +
                       formantData[j111|1] *      tL  *      tY  *      tX;*/
            const tet = FormantTet.nearest(sincFFTPeaks[0], sincFFTPeaks[1], sincFFTPeaks[2]);
            const coords = tet.barycentric_coords(sincFFTPeaks[0], sincFFTPeaks[1], sincFFTPeaks[2]);
            spaceCtx.fillStyle = coords[0] >= 0 && coords[1] >= 0 && coords[2] >= 0 && coords[3] >= 0 ? 'green' : 'red';
            for (let i = 0; i < 4; i++) {
                spaceCtx.beginPath();
                x = (1 - this.mF2 * Math.log2(formantData(tet.x[i][0], tet.x[i][1], tet.x[i][2], 1)) + this.bF2) * spaceCanvas.height;
                y = (    this.mF1 * Math.log2(formantData(tet.x[i][0], tet.x[i][1], tet.x[i][2], 0)) - this.bF1) * spaceCanvas.height;
                spaceCtx.arc(x, y, 10, 0, 2*Math.PI);
                spaceCtx.fill();
            }/*
            spaceCtx.fillStyle = projected ? 'blue' : 'green';
            spaceCtx.beginPath();
            x = (1 - this.mF2 * Math.log2(f1) + this.bF2) * spaceCanvas.height;
            y = (    this.mF1 * Math.log2(f0) - this.bF1) * spaceCanvas.height;
            spaceCtx.arc(x, y, 10, 0, 2*Math.PI);
            spaceCtx.fill();*/
        }

        formantDataCollection.step(sincFFTPeaks);
    }
}

Analysis.init();
