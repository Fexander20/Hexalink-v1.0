# Hexalink v2.0.0

Hexalink v2.0.0 provides a hexagonal finite-volume surface-routing implementation and a separate volume-based terrain-priority screening workflow.

Included are the authoritative Python routing engine, reproducible numerical checks against an analytical nonlinear-diffusion solution and an independently implemented Cartesian solver, a controlled synthetic St. John routing setup, event-based screening evaluation scripts, and selected summary results.

The numerics and screening are research methods. The St. John routing case uses synthetic forcing and is not a calibrated flood hindcast; the GFD comparison concerns screening ranking only. Source DEM, derived mesh and observed flood raster are not redistributed. See the README and data-availability statement for dependencies and limitations.

**Licence:** MIT for the original software. **Creator:** Ifeanyi Chidozie Oli.
