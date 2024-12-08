settings {
    logfile    = "/var/log/lsyncd/lsyncd.log",
    statusFile = "/var/log/lsyncd/lsyncd.status",
    inotifyMode = "CloseWrite or Modify"
}

sync {
    default.rsync,
    source = "/home/simulateurspatial/Documents/",
    target = "blue1@192.168.1.128:/home/blue1/Bureau/",
    rsync = {
        compress = true,
        archive = true,
        verbose = true,
        _extra = {"--delete", "--rsync-path=/usr/bin/rsync"}
    }
}
