import { NextResponse } from 'next/server';
import { exec } from 'child_process';
import { promisify } from 'util';

const execAsync = promisify(exec);

export async function GET() {
  try {
    // The Next.js app runs in cardiacFoamv2/frontend, so the driver workspace is one level up
    const backendCwd = process.cwd().replace('/frontend', '');
    
    // We execute the backend introspection layer
    const { stdout } = await execAsync('python3 -m openfoam_driver describe', {
      cwd: backendCwd
    });
    
    return NextResponse.json(JSON.parse(stdout));
  } catch (error) {
    return NextResponse.json({ error: String(error) }, { status: 500 });
  }
}
