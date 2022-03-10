// Reading time averaging interval
label avgInterval
(
    runTime.controlDict().lookupOrDefault<label>
        ("avgInterval", runTime.timeIndex())
);

// Reading time restart averaging time
scalar restartAvgTime
(
    runTime.controlDict().lookupOrDefault<scalar>
        ("restartAvgTime", runTime.endTime().value()+1)
);

// Reading maximum number of Explicit Cycles
label expItr(0);
label expItrMax(runTime.controlDict().lookupOrDefault<label>("expItrMax", P+1));

// Reading truncated P (Ptrunc)
label Ptrunc(runTime.controlDict().lookupOrDefault<label>("Ptrunc", P));
