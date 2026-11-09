# unit tests use test-drive, which conda-forge builds with gfortran
test_drive = dependency('test-drive', required : false)
if test_drive.found() and not fc_id.contains('intel')
  tests = [
    'Cbck12',
    'Gncn',
    'Smooth',
  ]

  test_srcs = files('tester.f90')
  foreach t : tests
    test_srcs += files('Test@0@.f90'.format(t))
  endforeach

  # compile the mfusg sources here too, so ninja orders their module files
  # before the tests that use them
  tester = executable(
    'tester',
    sources: test_srcs + mfusg_sources,
    dependencies: [test_drive] + extra_deps,
  )

  foreach t : tests
    test(t, tester, args: t)
  endforeach
endif
