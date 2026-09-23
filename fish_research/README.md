> Current API migration: see [API_MIGRATION.md](API_MIGRATION.md). The dry case now requires Reload.xml; the wet case is not yet adapted.

# Fixed-head fish research (fish_2)

## Cases and status
- `tests/2d_examples/test_2d_fish_spine_10bone`: latest dry model. Cached equal/opposite attachment forces, objective attachment damping and consistent reaction torque; head-to-tail active strain wave; configurable hinge damping (default 0.6).
- `tests/2d_examples/test_2d_fish_fixed_head_water`: existing experimental wet model, including frozen-body diagnostics. It has NOT yet received the dry attachment/drive corrections. Do not treat it as the validated wet counterpart of the dry case. Zero-drive fluid startup disturbances remain under investigation.

The two CMake targets have the same names as their directories. The parent examples CMake file discovers both automatically. Configure a NEW build directory for this checkout; do not reuse the old checkout's CMake cache.

## Reload particles and execution
The shared reference file `reload/MuscleBody_rld.xml` contains 10,078 relaxed muscle particles. Copy the `reload` directory into a separate run directory and run the chosen executable from that directory with `--relax=false --reload=true --state_recording=true`. Use separate run directories for different cases. Runtime DLLs must be available on PATH.

Dry reproduction settings (PowerShell):
```powershell
$env:FISH_END_TIME = '2'
$env:FISH_DRIVE_SCALE = '1'
$env:FISH_DT_SCALE = '1'
$env:FISH_HINGE_DAMPING = '0.6'
$env:FISH_WALL_SECONDS = '3600'
```
The default wall-time guard is 1800 seconds. The 2-second run needed 1827 seconds on the current machine; slower machines may need a larger guard.

## Latest dry result
The 2-second run completed with exit code 0. Sampled tail peak-to-peak displacement in cycles 5–8 was 23.913, 23.157, 22.065 and 18.648 mm. Cycle 8 differs from cycle 7 by -15.48% in amplitude and 3.658 mm RMS in the matched-phase displacement trace. Thus a stable periodic response has NOT been established. Later deformation/attachment behaviour needs inspection before accepting a baseline for water coupling.

The rear curvature proxy remains tail-lagging, but its fitted phase difference changes from about -28 to -31 degrees during 1.0–1.5 s to -42 to -48 degrees during 1.5–2.0 s. This is not evidence of a settled periodic state. See analysis.json and long_run_analysis.png. Statistics use 25 matched samples per 0.25-second cycle; last saved frame is 1.990002 s.

The attachment unit verification checks instantaneous force/torque/power identities and objectivity. This does not prove whole-system time-discrete energy conservation or biological validity.

Only source, reference particles, and a compact result summary are migrated. Compiled executables and full simulation output remain in their original run directories.
