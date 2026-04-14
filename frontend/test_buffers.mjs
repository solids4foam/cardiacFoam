import fs from 'fs';
import * as THREE from 'three';

const geomBuf = fs.readFileSync('./public/base_geometry.ply');
console.log("PLY Size (bytes):", geomBuf.length);

const dvBuf = fs.readFileSync('./public/scalars/disassemble_vectors.bin');
const dv = new Float32Array(dvBuf.buffer, dvBuf.byteOffset, dvBuf.length / 4);
console.log("Disassemble Vectors Array Length:", dv.length);
console.log("Expected Length (225139 * 3):", 225139 * 3);

let nans = 0;
for(let i=0; i<dv.length; i++) {
    if(Number.isNaN(dv[i])) nans++;
}
console.log("Disassemble Vectors NaNs:", nans);

const colorBuf = fs.readFileSync('./public/scalars/uvc_rotational.bin');
const colorBytes = new Uint8Array(colorBuf.buffer, colorBuf.byteOffset, colorBuf.length);
console.log("Color Bytes Length:", colorBytes.length);
console.log("Expected Length (225139 * 3):", 225139 * 3);
