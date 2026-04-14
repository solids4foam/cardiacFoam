"use client";

import { useState } from "react";
import dynamic from 'next/dynamic';
import Sidebar from "@/components/Sidebar";

const MeshViewer = dynamic(() => import('@/components/MeshViewer'), { 
  ssr: false, 
  loading: () => <div className="absolute inset-0 flex items-center justify-center text-slate-400 bg-slate-900">Loading 3D Engine...</div> 
});

export default function Home() {
  const [isZoomed, setIsZoomed] = useState(false);
  const [activeModel, setActiveModel] = useState("TNNP");
  
  const [activeField, setActiveField] = useState<string | null>(null);
  const [disassembleLevel, setDisassembleLevel] = useState(0);

  const handleFieldSelect = (field: string) => {
    setActiveField(field);
  };

  return (
    <main className="flex h-screen w-screen overflow-hidden bg-slate-900 relative">
      <Sidebar 
        onZoomToggle={() => setIsZoomed(!isZoomed)} 
        isZoomed={isZoomed}
        activeModel={activeModel}
        setActiveModel={setActiveModel}
        activeField={activeField}
        onSelectField={handleFieldSelect}
        disassembleLevel={disassembleLevel}
        onDisassembleChange={setDisassembleLevel}
      />
      <div className="flex-1 relative">
        <MeshViewer 
          activeModel={activeModel} 
          isZoomed={isZoomed} 
          activeField={activeField}
          disassembleLevel={disassembleLevel}
        />
        
        {/* Educational Overlay */}
        {isZoomed && (
          <div className="absolute bottom-8 left-1/2 -translate-x-1/2 max-w-lg bg-slate-800/80 backdrop-blur-md p-6 rounded-2xl border border-slate-700 shadow-2xl animate-in slide-in-from-bottom-8 fade-in text-center">
            <h3 className="text-xl font-bold text-white mb-2">Cellular Level Interaction</h3>
            <p className="text-slate-300 text-sm leading-relaxed">
              You are now looking at the cellular simulation governed by the <strong className="text-blue-400">{activeModel}</strong> model. 
              This model calculates the action potential based on exchanging `{activeModel === "TNNP" ? "Sodium, Calcium, and Potassium" : "simplified fast and slow"}` variables.
            </p>
          </div>
        )}
      </div>
    </main>
  );
}
