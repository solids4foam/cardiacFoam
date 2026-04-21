import { NextResponse } from 'next/server';
import { exec } from 'child_process';
import { promisify } from 'util';

const execAsync = promisify(exec);

export async function GET(request: Request) {
  try {
    // The Next.js app runs in cardiacFoamv2/frontend, so the driver workspace is one level up
    const backendCwd = process.cwd().replace('/frontend', '');
    const url = new URL(request.url);
    const entry = url.searchParams.get('entry') ?? url.searchParams.get('tutorial') ?? 'singleCell';
    const entryKind = url.searchParams.get('entryKind');

    // We execute the backend introspection layer
    const command = [
      'python3 -m openfoam_driver describe',
      `--entry ${JSON.stringify(entry)}`,
      entryKind ? `--entry-kind ${JSON.stringify(entryKind)}` : null,
    ].filter(Boolean).join(' ');

    const { stdout } = await execAsync(command, {
      cwd: backendCwd
    });
    
    return NextResponse.json(JSON.parse(stdout));
  } catch (error) {
    return NextResponse.json({ error: String(error) }, { status: 500 });
  }
}
