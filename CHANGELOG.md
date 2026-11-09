# Changelog

Features, fixes, and performance changes in each MODFLOW-USG release.

## Version 1.5.2 (2026-11-09)

### Fixes

* build: Link the Intel runtime libraries statically (#80) ([07f15a4](https://github.com/MODFLOW-ORG/mfusg/commit/07f15a41de90e53633fb634402c4aae5c27c5f9c))
* lpf: Avoid division by zero in anisotropic conductance (#78) ([1ac60bb](https://github.com/MODFLOW-ORG/mfusg/commit/1ac60bb85f450328f5b8bb5d6e0e229814aa4ab7))
* gncn: Allocate and deallocate the optional arrays (#81) ([7534514](https://github.com/MODFLOW-ORG/mfusg/commit/75345140a8df78e921f7961996592590030ed7af))
* cln: Compute CLN drawdown from the drawdown reference (#84) ([79e3641](https://github.com/MODFLOW-ORG/mfusg/commit/79e36412d141019bbfb77acc25cdc868a3005424))

## Version 1.5.1 (2026-11-09)

Changes since version 1.5.0.

### Features

* build: Build with meson and pixi, and build and release executables for Linux, macOS (Apple silicon), and Windows with GitHub Actions (#54, #57) ([6263424](https://github.com/MODFLOW-ORG/mfusg/commit/6263424), [694e07d](https://github.com/MODFLOW-ORG/mfusg/commit/694e07d))

### Fixes

* ets: Add NODES to GLOBAL import in GWF2ETS8U1RP (#47) ([b48bac7](https://github.com/MODFLOW-ORG/mfusg/commit/b48bac7))
* openspec: Use unformatted stream access so gfortran builds can write binary files (#65) ([2c69bb1](https://github.com/MODFLOW-ORG/mfusg/commit/2c69bb1))
* Initialize variables read before they are set in WEL, LAK, and SMS (#69) ([0069671](https://github.com/MODFLOW-ORG/mfusg/commit/0069671))
* drt: Compute layer, row, and column for printed budget rates (#71) ([89b30fb](https://github.com/MODFLOW-ORG/mfusg/commit/89b30fb))
