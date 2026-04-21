"use client";

import { useRef, Suspense, useMemo, useState, useEffect } from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls, Text } from "@react-three/drei";
import * as THREE from "three";

// Cellular zoom placeholder
function CellularAnimation({ activeModel }: { activeModel: string }) {
  const meshRef = useRef<THREE.Mesh>(null);
  useFrame((_, delta) => {
    if (meshRef.current) {
      meshRef.current.rotation.x += delta * 0.2;
      meshRef.current.rotation.y += delta * 0.5;
    }
  });
  return (
    <group>
      <mesh ref={meshRef}>
        <icosahedronGeometry args={[1, 1]} />
        <meshPhysicalMaterial color="#ff4757" roughness={0.2} transmission={0.9} thickness={1} />
      </mesh>
      <Text position={[0, -1.5, 0]} fontSize={0.2} color="white">
        Simulating {activeModel} Ion Exchange
      </Text>
    </group>
  );
}

function CameraRig() {
  return (
    <OrbitControls
      makeDefault
      enablePan={false}
      target={[0, 0, 0]}
      minPolarAngle={Math.PI / 2 - Math.PI / 12}  // ~75° — slight tilt up
      maxPolarAngle={Math.PI / 2 + Math.PI / 12}  // ~105° — slight tilt down
      minDistance={50}
      maxDistance={400}
      enableDamping
      dampingFactor={0.08}
    />
  );
}

// Full-featured heart mesh built from JSON arrays
function HeartMesh({ activeField, disassembleLevel }: {
  activeField: string | null,
  disassembleLevel: number
}) {
  const [geometry, setGeometry] = useState<THREE.BufferGeometry | null>(null);
  const [originalPositions, setOriginalPositions] = useState<Float32Array | null>(null);
  const [macroVectors, setMacroVectors] = useState<Float32Array | null>(null);
  const [microVectors, setMicroVectors] = useState<Float32Array | null>(null);

  // Load geometry from JSON (runs once)
  useEffect(() => {
    fetch('/heart_geometry.json')
      .then(res => res.json())
      .then(data => {
        const geom = new THREE.BufferGeometry();
        const positions = new Float32Array(data.positions);
        geom.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geom.computeVertexNormals();
        geom.computeBoundingSphere();
        setGeometry(geom);
        setOriginalPositions(new Float32Array(positions));
        console.log(`[Heart] Loaded: ${data.vertexCount} verts, ${data.faceCount} faces`);
      })
      .catch(console.error);
  }, []);

  // Load hierarchical explosion vectors (runs once)
  useEffect(() => {
    fetch('/scalars/disassemble_macro.bin')
      .then(res => res.arrayBuffer())
      .then(buf => setMacroVectors(new Float32Array(buf)))
      .catch(console.error);
    fetch('/scalars/disassemble_micro.bin')
      .then(res => res.arrayBuffer())
      .then(buf => setMicroVectors(new Float32Array(buf)))
      .catch(console.error);
  }, []);

  // Exploded view: macro (chambers) + micro (sub-parts) applied simultaneously
  useEffect(() => {
    if (!geometry || !originalPositions || !macroVectors || !microVectors) return;
    const pos = geometry.attributes.position.array as Float32Array;
    const t = disassembleLevel;
    for (let i = 0; i < pos.length; i++) {
      pos[i] = originalPositions[i] + macroVectors[i] * t + microVectors[i] * t * 0.3;
    }
    geometry.attributes.position.needsUpdate = true;
    geometry.computeVertexNormals();
  }, [geometry, disassembleLevel, originalPositions, macroVectors, microVectors]);

  // Load and apply scalar color field
  useEffect(() => {
    if (!geometry || !activeField) return;
    let active = true;
    fetch(`/scalars/${activeField}.bin`)
      .then(res => res.arrayBuffer())
      .then(buffer => {
        if (!active) return;
        const rgb = new Uint8Array(buffer);
        geometry.setAttribute('color', new THREE.BufferAttribute(rgb, 3, true));
        geometry.attributes.color.needsUpdate = true;
      })
      .catch(console.error);
    return () => { active = false; };
  }, [geometry, activeField]);

  const material = useMemo(() => new THREE.MeshStandardMaterial({
    vertexColors: true,
    roughness: 0.7,
    metalness: 0.1,
    side: THREE.DoubleSide,
    color: 0xffffff, // white: vertex colors show through unmodified
  }), []);

  if (!geometry) {
    return <Text position={[0, 0, 0]} fontSize={0.3} color="white">Building mesh...</Text>;
  }

  return (
    <mesh
      geometry={geometry}
      material={material}
      frustumCulled={false}
      onPointerOver={() => (document.body.style.cursor = 'crosshair')}
      onPointerOut={() => (document.body.style.cursor = 'auto')}
    />
  );
}

export default function MeshViewer({ activeModel, isZoomed, activeField, disassembleLevel }: {
  activeModel: string,
  isZoomed: boolean,
  activeField: string | null,
  disassembleLevel: number
}) {


  return (
    <div className="w-full h-full bg-slate-900 absolute inset-0">
      <Canvas camera={{ position: [0, 0, 200], fov: 45 }}>
        <ambientLight intensity={0.6} />
        <directionalLight position={[100, 100, 50]} intensity={1} />
        <directionalLight position={[-50, -50, -50]} intensity={0.3} />

        <Suspense fallback={<Text position={[0, 0, 0]} fontSize={0.5} color="white">Loading...</Text>}>
          {isZoomed ? (
            <CellularAnimation activeModel={activeModel} />
          ) : (
            <HeartMesh activeField={activeField} disassembleLevel={disassembleLevel} />
          )}
        </Suspense>

        <CameraRig />
      </Canvas>
    </div>
  );
}
