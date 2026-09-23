# Dry fish migration to current SPHinXsys API

Only the dry case was adapted. No library source or physical parameters were changed.

Changes: BoundingBoxd; MultiPolygon.addPolygon and GeometricOps; defineMatterMaterial; reference-returning defineBodyLevelSetShape and parameterless writeLevelSet; body.Name(); public Simbody mass-property/center accessors; obtain system time through the owning body.

## Reload format
The current library reads reload/Reload.xml with a reload_data root and a body element named MuscleBody. The legacy MuscleBody_rld.xml remains intact for the older executable. Reload.xml contains the same 10,078 particle attribute dictionaries (positions, volumes and IDs), verified after conversion. Copy BOTH files when preparing runs for old/new comparisons. The VS Code runs/dry_0p5/reload folder has already been prepared.

## Validation
The dry target built successfully against the new checkout. Both old and new executables completed 0.03 s with drive scale 1, timestep scale 1 and hinge damping 0.6, using identical reference particles. Exit codes were zero. Matched saved frames had at most 5.023e-7 m position discrepancy; hinge9 differed by at most 1.162e-10 rad. Position differences include VTP serialization precision and should not be interpreted as exact internal solver differences. See api_migration_comparison.json.

This is a startup compatibility check, not a multi-cycle convergence test. The earlier 2 s run showed late-cycle drift; that issue is not resolved by API migration. The wet case remains unadapted and must not be selected for compilation yet. Build the dry target explicitly rather than all.

## VS Code
Select Windows ck Release; configure; build test_2d_fish_spine_10bone. Select Dry fish - 0.5 seconds in Run and Debug and press Ctrl+F5. The launch profile does not build automatically. Results go to runs/dry_0p5; preserve previous output before reusing the same directory.
