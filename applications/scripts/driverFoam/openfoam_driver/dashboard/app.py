from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path

from fastapi import FastAPI, Form, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .assets3d import export_glb, find_primary_surface
from .catalog import build_catalog
from .store import AnnotationStore

_PKG = Path(__file__).parent


def create_app(root: Path, store_path: Path, enable_3d: bool = True,
               blender: str | None = None) -> FastAPI:
    root = Path(root)
    store = AnnotationStore(store_path)
    assets_dir = Path(store_path).parent / "assets"
    templates = Jinja2Templates(directory=str(_PKG / "templates"))
    app = FastAPI()
    app.mount("/static", StaticFiles(directory=str(_PKG / "static")), name="static")

    def _catalog() -> list:
        return list(build_catalog(root))

    def _card(case_id: str):
        for c in _catalog():
            if c.case_id == case_id:
                return c
        return None

    def _glb_url(card) -> str | None:
        if not (enable_3d and card.is_3d):
            return None
        surface = find_primary_surface(root / card.case_id)
        if surface is None:
            return None
        out = assets_dir / card.case_id / "model.glb"
        result = export_glb(surface, out)
        if result is None:
            return None
        return f"/asset3d/{card.case_id}/model.glb"

    @app.get("/", response_class=HTMLResponse)
    def inventory(request: Request):
        cards = _catalog()
        families: "OrderedDict[str, list]" = OrderedDict()
        for c in sorted(cards, key=lambda c: (c.family, c.name)):
            families.setdefault(c.family, []).append(c)
        ctx = {
            "families": families,
            "total": len(cards),
            "ran": sum(1 for c in cards if c.status == "ran"),
            "threed": sum(1 for c in cards if c.is_3d),
        }
        return templates.TemplateResponse(request, "inventory.html", ctx)

    @app.get("/case/{case_id:path}", response_class=HTMLResponse)
    def case_detail(request: Request, case_id: str):
        card = _card(case_id)
        if card is None:
            return HTMLResponse("case not found", status_code=404)
        ann = store.get_annotation(case_id)
        ctx = {
            "card": card,
            "notes": ann.notes, "tags": ann.tags,
            "captions": store.get_captions(case_id),
            "glb_url": _glb_url(card),
        }
        return templates.TemplateResponse(request, "case.html", ctx)

    @app.post("/case/{case_id:path}/annotation")
    def save_annotation(case_id: str, notes: str = Form(""), tags: str = Form("")):
        store.upsert_annotation(case_id, notes=notes,
                                tags=[t for t in tags.split(",") if t.strip()])
        return RedirectResponse(f"/case/{case_id}", status_code=303)

    @app.post("/case/{case_id:path}/caption")
    def save_caption(case_id: str, figure_id: str = Form(...), caption: str = Form("")):
        store.upsert_caption(case_id, figure_id, caption)
        return RedirectResponse(f"/case/{case_id}", status_code=303)

    @app.post("/rescan")
    def rescan():
        return RedirectResponse("/", status_code=303)

    @app.get("/api/catalog.json")
    def catalog_json():
        cards = _catalog()
        payload = {"root": str(root),
                   "cases": [c.to_dict() for c in cards]}
        return JSONResponse(json.loads(json.dumps(payload)))

    @app.get("/asset/{relpath:path}")
    def serve_asset(relpath: str):
        from fastapi.responses import FileResponse
        target = (root / relpath).resolve()
        if root.resolve() not in target.parents and target != root.resolve():
            return HTMLResponse("forbidden", status_code=403)
        if not target.is_file():
            return HTMLResponse("not found", status_code=404)
        return FileResponse(target)

    @app.get("/asset3d/{relpath:path}")
    def serve_asset3d(relpath: str):
        from fastapi.responses import FileResponse
        target = (assets_dir / relpath).resolve()
        if assets_dir.resolve() not in target.parents:
            return HTMLResponse("forbidden", status_code=403)
        if not target.is_file():
            return HTMLResponse("not found", status_code=404)
        return FileResponse(target)

    return app
