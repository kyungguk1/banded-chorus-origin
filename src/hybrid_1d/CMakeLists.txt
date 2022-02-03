
add_executable(hybrid_1d)
target_sources(hybrid_1d PRIVATE
    ParamSet.cc
    Core/BField.cc
    Core/EField.cc
    Core/Charge.cc
    Core/Current.cc
    Core/Species.cc
    Core/ExternalSource.cc
    Core/PartSpecies.cc
    Core/ColdSpecies.cc
    Core/Domain.cc
    Core/Domain_PC.cc
    Core/Domain_CAMCL.cc
    Boundary/Delegate.cc
    Boundary/SubdomainDelegate.cc
    Boundary/DistributedParticleDelegate.cc
    Boundary/WorkerDelegate.cc
    Boundary/MasterDelegate.cc
    Recorder/Snapshot.cc
    Recorder/Recorder.cc
    Recorder/EnergyRecorder.cc
    Recorder/FieldRecorder.cc
    Recorder/MomentRecorder.cc
    Recorder/ParticleRecorder.cc
    Recorder/VHistogramRecorder.cc
    Driver.cc
    main.cc
    )

set_project_warnings(hybrid_1d)
enable_sanitizers(hybrid_1d)
