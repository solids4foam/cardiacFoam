"use client";

import { useState, useEffect } from "react";
import { Play, Settings2, Loader2 } from "lucide-react";

interface SidebarProps {
  onZoomToggle: () => void;
  activeModel: string;
  setActiveModel: (m: string) => void;
  isZoomed: boolean;
  activeField: string | null;
  onSelectField: (f: string) => void;
  disassembleLevel: number;
  onDisassembleChange: (val: number) => void;
}

interface PatchInfo {
  name: string;
  label: string;
  color: number[];
  group: string;
}

export default function Sidebar({ onZoomToggle, activeModel, setActiveModel, isZoomed, activeField, onSelectField, disassembleLevel, onDisassembleChange }: SidebarProps) {
  const [meshVars, setMeshVars] = useState<string[]>([]);
  const [loadingVars, setLoadingVars] = useState(true);
  const [patchGroups, setPatchGroups] = useState<Record<string, PatchInfo[]>>({});

  useEffect(() => {
    fetch('/scalars/manifest.json')
      .then(res => res.json())
      .then(data => {
        if (data.fields) {
          setMeshVars(Object.keys(data.fields));
        }
      })
      .catch(console.error)
      .finally(() => setLoadingVars(false));
  }, []);

  // Load anatomical patch groupings
  useEffect(() => {
    fetch('/scalars/patch_colors.json')
      .then(res => res.json())
      .then((data: { chambers?: Record<string, { patches: string[] }>; patches?: Record<string, { label: string; color: number[]; group: string }> }) => {
        if (!data.chambers || !data.patches) return;
        const groups: Record<string, PatchInfo[]> = {};
        for (const [chamberName, chamberData] of Object.entries(data.chambers)) {
          groups[chamberName] = chamberData.patches
            .filter((pname: string) => data.patches![pname])
            .map((pname: string) => ({
              name: pname,
              label: data.patches![pname].label,
              color: data.patches![pname].color,
              group: data.patches![pname].group,
            }));
        }
        setPatchGroups(groups);
      })
      .catch(console.error);
  }, []);

  return (
    <div className="w-80 h-full bg-slate-800 text-slate-200 flex flex-col shadow-2xl z-10 p-6 border-r border-slate-700 overflow-y-auto">
      <div className="flex items-center gap-3 mb-8">
        <div className="w-10 h-10 rounded-lg bg-blue-500 flex items-center justify-center">
          <Settings2 size={20} className="text-white" />
        </div>
        <div>
          <h1 className="text-xl font-bold tracking-tight text-white">driverFoam</h1>
          <p className="text-xs text-slate-400">Educational Viewer</p>
        </div>
      </div>

      <div className="flex-1 space-y-6">
        <div>
          <h2 className="text-sm font-semibold text-slate-400 uppercase tracking-wider mb-3">Simulation Config</h2>
          
          <div className="space-y-6">
            <div>
              <label className="block text-sm font-medium text-slate-300 mb-1">Mesh</label>
              <div className="px-3 py-2 bg-slate-900 rounded-md text-sm border border-slate-700 text-slate-400">
                ASCIIlegacy.vtk
              </div>
            </div>

            {/* Render Mesh Variables dynamically fetched from the API */}
            <div className="pt-2 border-t border-slate-700">
               <label className="block text-sm font-medium text-slate-300 mb-2">Available Mesh Variables</label>
               {loadingVars ? (
                 <div className="flex items-center gap-2 text-slate-400 text-sm">
                   <Loader2 size={14} className="animate-spin" /> Extracting fields...
                 </div>
               ) : (
                 <div className="flex flex-wrap gap-2">
                   {meshVars.length > 0 ? meshVars.map(v => (
                     <button 
                       key={v} 
                       onClick={() => onSelectField(v)}
                       className={`px-2 py-1 text-xs rounded border transition-colors ${
                         activeField === v 
                           ? "bg-blue-600 border-blue-500 text-white" 
                           : "bg-slate-700 text-blue-200 border-slate-600 hover:bg-slate-600"
                       }`}
                     >
                       {v}
                     </button>
                   )) : (
                     <span className="text-xs text-slate-500 italic">No point/cell data found</span>
                   )}
                 </div>
               )}
            </div>

            {/* Anatomical Tag Legend */}
            <div className="pt-2 border-t border-slate-700">
              <label className="block text-sm font-medium text-slate-300 mb-2">Anatomical Regions</label>
              {Object.keys(patchGroups).length > 0 ? (
                <div className="space-y-1">
                  {Object.entries(patchGroups).map(([groupName, patches]) => (
                    <details key={groupName} className="group">
                      <summary className="flex items-center gap-2 cursor-pointer px-2 py-1.5 rounded hover:bg-slate-700/50 text-sm text-slate-300 select-none">
                        <span className="text-[10px] text-slate-500 group-open:rotate-90 transition-transform">▶</span>
                        <span className="font-medium">{groupName}</span>
                        <span className="text-[10px] text-slate-500 ml-auto">{patches.length}</span>
                      </summary>
                      <div className="ml-5 mt-1 space-y-0.5">
                        {patches.map(p => (
                          <div key={p.name} className="flex items-center gap-2 px-2 py-0.5 text-xs text-slate-400">
                            <span
                              className="inline-block w-3 h-3 rounded-sm flex-shrink-0"
                              style={{ backgroundColor: `rgb(${p.color.join(',')})` }}
                            />
                            {p.label}
                          </div>
                        ))}
                      </div>
                    </details>
                  ))}
                </div>
              ) : (
                <span className="text-xs text-slate-500 italic">Loading regions...</span>
              )}
            </div>

            <div className="pt-2 border-t border-slate-700">
              <label className="block text-sm font-medium text-slate-300 mb-1">Ionic Model</label>
              <select 
                value={activeModel}
                onChange={(e) => setActiveModel(e.target.value)}
                className="w-full px-3 py-2 bg-slate-900 rounded-md text-sm border border-slate-700 focus:outline-none focus:ring-2 focus:ring-blue-500"
              >
                <option value="TNNP">TNNP (Detailed)</option>
                <option value="FentonKarma">Fenton-Karma (Phenomenological)</option>
                <option value="FitzHughNagumo">FitzHugh-Nagumo</option>
              </select>
            </div>
            
            <div className="pt-4 border-t border-slate-700">
              <div className="flex justify-between items-center mb-2">
                <label className="text-sm font-medium text-slate-300">Exploded View</label>
                <span className="text-xs text-blue-400 font-mono">{disassembleLevel.toFixed(2)}x</span>
              </div>
              <input 
                type="range" 
                min="0" 
                max="2.5" 
                step="0.05"
                value={disassembleLevel}
                onChange={(e) => onDisassembleChange(parseFloat(e.target.value))}
                className="w-full accent-blue-500 cursor-pointer"
              />
              <p className="text-[10px] text-slate-500 mt-1 leading-tight">
                0–1: Ventricles &amp; atria separate. 1–2.5: Valves &amp; annulus rings separate within each chamber.
              </p>
            </div>
            
            <button 
              onClick={onZoomToggle}
              className="mt-4 w-full py-2 px-4 rounded-md border border-blue-500/50 text-blue-400 text-sm hover:bg-blue-500/10 transition-colors flex items-center justify-center gap-2"
            >
              {isZoomed ? "Exit Microscope View" : "Explore Tissue (Microscope)"}
            </button>
          </div>
        </div>
      </div>

      <div className="mt-8">
        <button className="w-full py-3 px-4 bg-blue-600 hover:bg-blue-500 text-white rounded-lg font-medium transition-colors flex items-center justify-center gap-2 shadow-lg shadow-blue-500/20">
          <Play size={18} />
          Run Simulation
        </button>
      </div>
    </div>
  );
}
