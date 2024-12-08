settings {
    logfile    = "/var/log/lsyncd/lsyncd.log",
    statusFile = "/var/log/lsyncd/lsyncd.status",
    nodaemon   = false,
}

sync {
    default.rsync,
    source = "/home/simulateurspatial/Documents",
    target = "simulateurspatial@192.168.1.126:/home/simulateurspatial/jsatorb",
    rsync = {
        archive = true,
        compress = true,
        verbose = true,
        _extra = {"--delete"}
    }
}



