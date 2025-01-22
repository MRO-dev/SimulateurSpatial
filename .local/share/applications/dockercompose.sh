#!/usr/bin/env bash

JSATORB_PATH="/home/simulateurspatial/jsatorb"
LOG_FILE="/tmp/docker-script.log"

# Function for logging
log_message() {
    local message="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] $message" | tee -a "$LOG_FILE"
}

# Clear previous log
> "$LOG_FILE"

# Check if the directory exists
if [ ! -d "$JSATORB_PATH" ]; then
    log_message "Error: Directory $JSATORB_PATH does not exist"
    exit 1
fi

log_message "Starting script execution..."

# 1. Launch Docker Compose in a new terminal
log_message "Launching Docker Compose in new terminal..."
gnome-terminal -- bash -c "
    cd \"$JSATORB_PATH\"
    mkdir -p jsatorb-rest-api/files/default_Sat/Data/Backup/
    docker compose up
    exec bash
"

# 2 & 3. Wait for containers to be ready then remove backup
(
    cd "$JSATORB_PATH"
    log_message "Background process started"
    log_message "Waiting for containers to be ready..."
    
    attempt=1
    max_attempts=60  # 30 seconds maximum wait time
    
    # First wait for containers to be created
    sleep 15
    log_message "Checking container status..."
    
    # Check if containers exist
    if ! docker compose ps > /dev/null 2>&1; then
        log_message "ERROR: Docker Compose not running"
        exit 1
    fi
    
    log_message "Docker Compose is well running. Checking all containers..."
    docker compose ps >> "$LOG_FILE"
    
    log_message "Removing backup folder..."
    if rm -rf "$JSATORB_PATH/jsatorb-rest-api/files/default_Sat/Data/Backup/"; then
        log_message "Backup folder successfully removed"
    else
        log_message "ERROR: Failed to remove backup folder"
    fi
    
    log_message "Background process completed"
) &

log_message "Main script completed. Check $LOG_FILE for detailed logs"
