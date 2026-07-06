from __future__ import annotations

import webbrowser
from pathlib import Path


def serve(root: Path, store_path: Path, host: str = "127.0.0.1", port: int = 8765,
          enable_3d: bool = True, blender: str | None = None,
          open_browser: bool = True) -> int:
    import uvicorn

    from .app import create_app

    app = create_app(root=Path(root), store_path=Path(store_path),
                     enable_3d=enable_3d, blender=blender)
    if open_browser:
        webbrowser.open(f"http://{host}:{port}/")
    uvicorn.run(app, host=host, port=port, log_level="info")
    return 0
